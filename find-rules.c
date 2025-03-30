#include "Apriori.h"
#include "dynamic_hash_table.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include <time.h>

#include "io-processing.c"

// Structure to hold performance metrics
typedef struct {
  double file_split_time;
  double frequent_items_time;
  double triangle_matrix_build_time;
  double pair_generation_time;
  double triangle_prune_time;
  double large_itemset_time;
  double rule_generation_time;
  double total_time;
  int total_transactions;
  int total_frequent_items;
  int total_frequent_pairs;
  int total_frequent_itemsets;
  int total_rules;
} PerformanceMetrics;

// Create a new frequent itemset of size k+1 from two itemsets of size k
ItemSet *generate_candidate(ItemSet *itemset1, ItemSet *itemset2, int k) {
  // Ensure the first k-1 elements match
  for (int i = 0; i < k - 1; i++) {
    if (itemset1->elements[i] != itemset2->elements[i]) {
      return NULL;
    }
  }
  
  // The last elements must be different, and in ascending order
  if (itemset1->elements[k - 1] >= itemset2->elements[k - 1]) {
    return NULL;
  }
  
  // Create a new candidate itemset
  ItemSet *candidate = malloc(sizeof(ItemSet));
  candidate->size = k + 1;
  candidate->elements = malloc((k + 1) * sizeof(int));
  candidate->support = 0.0;
  
  // Copy the k-1 common elements
  for (int i = 0; i < k - 1; i++) {
    candidate->elements[i] = itemset1->elements[i];
  }
  
  // Add the two different elements at the end
  candidate->elements[k - 1] = itemset1->elements[k - 1];
  candidate->elements[k] = itemset2->elements[k - 1];
  
  return candidate;
}

// Check if all subsets of the candidate itemset are frequent
int has_frequent_subsets(ItemSet *candidate, ItemSet *frequent_itemsets, int count, int k) {
  int *subset = malloc((k) * sizeof(int));
  
  // Generate all subsets of size k by excluding one element at a time
  for (int i = 0; i < k + 1; i++) {
    int pos = 0;
    for (int j = 0; j < k + 1; j++) {
      if (j != i) {
        subset[pos++] = candidate->elements[j];
      }
    }
    
    // Check if this subset is in the frequent itemsets
    int found = 0;
    for (int j = 0; j < count; j++) {
      if (frequent_itemsets[j].size == k) {
        int match = 1;
        for (int l = 0; l < k; l++) {
          int elem_match = 0;
          for (int m = 0; m < k; m++) {
            if (subset[l] == frequent_itemsets[j].elements[m]) {
              elem_match = 1;
              break;
            }
          }
          if (!elem_match) {
            match = 0;
            break;
          }
        }
        if (match) {
          found = 1;
          break;
        }
      }
    }
    
    if (!found) {
      free(subset);
      return 0;
    }
  }
  
  free(subset);
  return 1;
}

// Check if a candidate itemset is present in a transaction
int is_subset(int *candidate, int candidate_size, int *transaction, int transaction_size) {
  int matches = 0;
  
  for (int i = 0; i < candidate_size; i++) {
    for (int j = 0; j < transaction_size; j++) {
      if (candidate[i] == transaction[j]) {
        matches++;
        break;
      }
    }
  }
  
  return (matches == candidate_size);
}

// Print a single association rule
void print_rule(AssociationRule rule) {
  printf("{");
  for (int i = 0; i < rule.antecedent_size; i++) {
    printf("%d", rule.antecedent[i]);
    if (i < rule.antecedent_size - 1) {
      printf(", ");
    }
  }
  printf("} => {");
  for (int i = 0; i < rule.consequent_size; i++) {
    printf("%d", rule.consequent[i]);
    if (i < rule.consequent_size - 1) {
      printf(", ");
    }
  }
  printf("} (Confidence: %.2f%%, Support: %.2f%%)\n", 
         rule.confidence * 100, rule.support * 100);
}

int main(int argc, char **argv) {
  int size, rank;
  double start_time, end_time;
  PerformanceMetrics metrics = {0};
  
  // Record total execution start time
  start_time = MPI_Wtime();
  
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (argc < 4) {
    if (rank == 0)
      fprintf(stderr, "Usage: %s <inputfile> <support_threshold> <confidence_threshold>\n", argv[0]);
    MPI_Finalize();
    return 1;
  }

  float global_support = strtof(argv[2], NULL);
  float confidence_threshold = strtof(argv[3], NULL);

  // Split file into roughly equal sized portions for processing
  double phase_start = MPI_Wtime();
  long *split_points = NULL;
  split_points = (long *)malloc(size * sizeof(long));
  if (rank == 0) {
    printf("#################### \n");
    printf("Begin input processing \n");
    split_file(argv[1], split_points, size);
    printf("Split points generated \n");
    for (int i = 0; i < size; i++) {
      printf("split %i is %li \n", i, split_points[i]);
    }
    printf("#################### \n");
    metrics.file_split_time = MPI_Wtime() - phase_start;
  }

  // Send split points to all processors
  MPI_Bcast(split_points, size, MPI_LONG, 0, MPI_COMM_WORLD);
  
  // Broadcast file split timing to all processes
  MPI_Bcast(&metrics.file_split_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  // Each processor finds frequent ints, support and transaction counts
  phase_start = MPI_Wtime();
  HashTable local_table;
  init_table(&local_table, INT_TYPE);
  int local_transaction_count;
  local_transaction_count =
      process_chunk(argv[1], split_points, rank, &local_table, global_support);
  printf("Proc %i has finished processing %i lines in its file chunk \n", rank,
         local_transaction_count);
  metrics.frequent_items_time = MPI_Wtime() - phase_start;
  metrics.total_frequent_items = local_table.count;

  // Get total number of transactions
  int total_transactions;
  MPI_Reduce(&local_transaction_count, &total_transactions, 1, MPI_INT, MPI_SUM,
             0, MPI_COMM_WORLD);
  metrics.total_transactions = total_transactions;

  // Gather frequent items metrics
  int total_frequent_items;
  MPI_Reduce(&local_table.count, &total_frequent_items, 1, MPI_INT, MPI_SUM,
             0, MPI_COMM_WORLD);
  
  if (rank == 0) {
    printf("#################### \n");
    printf("Total of %i transactions found \n", total_transactions);
    printf("Total of %i frequent items found\n", total_frequent_items);
    printf("#################### \n");
    printf("Input processing complete \n");
    printf("#################### \n");
    printf("Beginning SON by performing Apriori on each processor \n");
  }
  
  MPI_Barrier(MPI_COMM_WORLD);

  // Build triangular matrix for storing pairs
  phase_start = MPI_Wtime();
  TriangularMatrix *local_tri = build_tri_matrix(&local_table);
  metrics.triangle_matrix_build_time = MPI_Wtime() - phase_start;
  printf("Proc %i: Triangular matrix built \n", rank);

  // Set file buffer for reading the assigned chunk
  FILE *fp = fopen(argv[1], "r");
  long read_start = (rank == 0) ? 0 : split_points[rank - 1];
  long read_end = split_points[rank];
  
  if (read_start > 0) {
    fseek(fp, read_start, SEEK_SET);
  } else {
    rewind(fp);
  }

  // Read each line in the chunk, entering candidate pairs into triangle matrix
  phase_start = MPI_Wtime();
  char buffer[2048];
  while (ftell(fp) < read_end && fgets(buffer, sizeof(buffer), fp) != NULL) {
    int outcount;
    int *frequent_items = extract_frequent(&local_table, buffer, &outcount);
    
    if (outcount >= 2) {
      check_pairs(local_tri, frequent_items, outcount);
    }
    
    free(frequent_items);
  }
  metrics.pair_generation_time = MPI_Wtime() - phase_start;
  printf("Proc %i: Pairs generated\n", rank);

  // Create list of supported pairs from triangle matrix
  phase_start = MPI_Wtime();
  int valid_pairs;
  float local_support = global_support * ((float)local_transaction_count / total_transactions);
  ItemSet *frequent_pairs = prune_triangle(local_tri, local_support, local_transaction_count, &valid_pairs);
  metrics.triangle_prune_time = MPI_Wtime() - phase_start;
  metrics.total_frequent_pairs = valid_pairs;
  
  printf("Proc %i: Triangle pruned, found %d valid pairs\n", rank, valid_pairs);

  // Initialize array for storing all frequent itemsets (pairs, triples, etc.)
  int max_itemsets = valid_pairs * 10; // Initial guess
  ItemSet *all_frequent_itemsets = malloc(max_itemsets * sizeof(ItemSet));
  int total_itemsets = 0;
  
  // Copy the frequent pairs to the all_frequent_itemsets array
  for (int i = 0; i < valid_pairs; i++) {
    all_frequent_itemsets[total_itemsets++] = frequent_pairs[i];
  }
  
  // Find larger itemsets (k > 2)
  phase_start = MPI_Wtime();
  int k = 2; // Current itemset size
  int prev_start_idx = 0;
  int prev_count = valid_pairs;
  
  while (prev_count > 0) {
    int new_start_idx = total_itemsets;
    int candidates_count = 0;
    ItemSet *candidates = NULL;
    
    // 1. Generate candidate itemsets of size k+1
    for (int i = prev_start_idx; i < prev_start_idx + prev_count - 1; i++) {
      for (int j = i + 1; j < prev_start_idx + prev_count; j++) {
        ItemSet *candidate = generate_candidate(&all_frequent_itemsets[i], 
                                               &all_frequent_itemsets[j], k);
        
        if (candidate != NULL && has_frequent_subsets(candidate, all_frequent_itemsets, 
                                                    total_itemsets, k)) {
          // Add the candidate to our list
          candidates = realloc(candidates, (candidates_count + 1) * sizeof(ItemSet));
          candidates[candidates_count++] = *candidate;
          free(candidate); // Free the container, but not the elements (copied above)
        }
      }
    }
    
    if (candidates_count == 0) {
      break; // No new candidates, we're done
    }
    
    // 2. Count support for each candidate
    // Reset file pointer to read the chunk again
    fseek(fp, read_start, SEEK_SET);
    
    while (ftell(fp) < read_end && fgets(buffer, sizeof(buffer), fp) != NULL) {
      int outcount;
      int *frequent_items = extract_frequent(&local_table, buffer, &outcount);
      
      if (outcount >= k + 1) {
        for (int i = 0; i < candidates_count; i++) {
          if (is_subset(candidates[i].elements, candidates[i].size, 
                       frequent_items, outcount)) {
            candidates[i].support += 1.0 / local_transaction_count;
          }
        }
      }
      
      free(frequent_items);
    }
    
    // 3. Prune candidates based on support
    int new_count = 0;
    for (int i = 0; i < candidates_count; i++) {
      if (candidates[i].support >= local_support) {
        // Check if we need to resize all_frequent_itemsets
        if (total_itemsets >= max_itemsets) {
          max_itemsets *= 2;
          all_frequent_itemsets = realloc(all_frequent_itemsets, 
                                        max_itemsets * sizeof(ItemSet));
        }
        
        // Add to frequent itemsets
        all_frequent_itemsets[total_itemsets++] = candidates[i];
        new_count++;
      } else {
        // Free memory for non-frequent candidates
        free(candidates[i].elements);
      }
    }
    
    // Update for next iteration
    prev_start_idx = new_start_idx;
    prev_count = new_count;
    k++;
    
    free(candidates);
    
    printf("Proc %i: Found %d frequent itemsets of size %d\n", rank, new_count, k);
  }
  
  metrics.large_itemset_time = MPI_Wtime() - phase_start;
  metrics.total_frequent_itemsets = total_itemsets;
  
  // Clean up
  fclose(fp);
  
  // Share frequent itemsets with master process
  int *local_itemset_counts = NULL;
  int *itemset_displacements = NULL;
  int total_all_itemsets = 0;
  ItemSet *global_itemsets = NULL;
  
  // First, gather the counts from all processes
  if (rank == 0) {
    local_itemset_counts = malloc(size * sizeof(int));
  }
  
  MPI_Gather(&total_itemsets, 1, MPI_INT, local_itemset_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  // Master process prepares to receive all itemsets
  if (rank == 0) {
    // Calculate total count and displacements for MPI_Gatherv
    itemset_displacements = malloc(size * sizeof(int));
    int disp = 0;
    
    for (int i = 0; i < size; i++) {
      itemset_displacements[i] = disp;
      total_all_itemsets += local_itemset_counts[i];
      disp += local_itemset_counts[i];
    }
    
    global_itemsets = malloc(total_all_itemsets * sizeof(ItemSet));
    printf("Master will receive %d total itemsets\n", total_all_itemsets);
  }
  
  // Note: We can't easily create an MPI datatype for ItemSet because it contains a pointer
  // For simplicity, we'll send each itemset's components individually
  
  if (rank == 0) {
    // Master process receives itemsets one by one
    int current_idx = 0;
    
    for (int src = 0; src < size; src++) {
      for (int i = 0; i < local_itemset_counts[src]; i++) {
        int size_buf;
        float support_buf;
        
        // Receive size and support
        MPI_Recv(&size_buf, 1, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&support_buf, 1, MPI_FLOAT, src, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        // Allocate elements array and receive elements
        global_itemsets[current_idx].size = size_buf;
        global_itemsets[current_idx].support = support_buf;
        global_itemsets[current_idx].elements = malloc(size_buf * sizeof(int));
        
        MPI_Recv(global_itemsets[current_idx].elements, size_buf, MPI_INT, 
                src, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        current_idx++;
      }
    }
    
    printf("Master has received all %d itemsets\n", total_all_itemsets);
  } else {
    // Worker processes send their itemsets
    for (int i = 0; i < total_itemsets; i++) {
      MPI_Send(&all_frequent_itemsets[i].size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&all_frequent_itemsets[i].support, 1, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
      MPI_Send(all_frequent_itemsets[i].elements, all_frequent_itemsets[i].size, 
              MPI_INT, 0, 2, MPI_COMM_WORLD);
    }
  }
  
  // Master process generates association rules
  if (rank == 0) {
    phase_start = MPI_Wtime();
    int rule_count;
    AssociationRule *rules = generate_rules(global_itemsets, total_all_itemsets, 
                                          confidence_threshold, &rule_count);
    metrics.rule_generation_time = MPI_Wtime() - phase_start;
    metrics.total_rules = rule_count;
    
    printf("#################### \n");
    printf("Generated %d association rules with confidence >= %f\n", 
          rule_count, confidence_threshold);
    printf("#################### \n");
    
    // Print top rules (limited to 10 for brevity)
    int rules_to_print = (rule_count < 10) ? rule_count : 10;
    printf("Top %d association rules:\n", rules_to_print);
    
    for (int i = 0; i < rules_to_print; i++) {
      printf("Rule %d: ", i + 1);
      print_rule(rules[i]);
    }
    
    // Clean up rules
    free_rules(rules, rule_count);
    
    // Clean up global itemsets
    for (int i = 0; i < total_all_itemsets; i++) {
      free(global_itemsets[i].elements);
    }
    free(global_itemsets);
    free(local_itemset_counts);
    free(itemset_displacements);
  }
  
  // Clean up local resources
  for (int i = 0; i < total_itemsets; i++) {
    free(all_frequent_itemsets[i].elements);
  }
  free(all_frequent_itemsets);
  free(frequent_pairs);
  free(local_tri->matrix);
  free(local_tri->item_to_index_map);
  free(local_tri->index_to_item_map);
  free(local_tri);
  free_table(&local_table);
  free(split_points);
  
  // Calculate total time
  end_time = MPI_Wtime();
  metrics.total_time = end_time - start_time;
  
  // Output performance metrics
  if (rank == 0) {
    printf("\n#################### \n");
    printf("Performance Metrics:\n");
    printf("Total execution time: %.3f seconds\n", metrics.total_time);
    printf("File splitting time: %.3f seconds (%.2f%%)\n", 
           metrics.file_split_time, 
           (metrics.file_split_time / metrics.total_time) * 100);
    printf("Frequent items mining time: %.3f seconds (%.2f%%)\n", 
           metrics.frequent_items_time, 
           (metrics.frequent_items_time / metrics.total_time) * 100);
    printf("Triangle matrix build time: %.3f seconds (%.2f%%)\n", 
           metrics.triangle_matrix_build_time, 
           (metrics.triangle_matrix_build_time / metrics.total_time) * 100);
    printf("Pair generation time: %.3f seconds (%.2f%%)\n", 
           metrics.pair_generation_time, 
           (metrics.pair_generation_time / metrics.total_time) * 100);
    printf("Triangle pruning time: %.3f seconds (%.2f%%)\n", 
           metrics.triangle_prune_time, 
           (metrics.triangle_prune_time / metrics.total_time) * 100);
    printf("Large itemset mining time: %.3f seconds (%.2f%%)\n", 
           metrics.large_itemset_time, 
           (metrics.large_itemset_time / metrics.total_time) * 100);
    printf("Rule generation time: %.3f seconds (%.2f%%)\n", 
           metrics.rule_generation_time, 
           (metrics.rule_generation_time / metrics.total_time) * 100);
    printf("\n");
    printf("Dataset statistics:\n");
    printf("Total transactions: %d\n", metrics.total_transactions);
    printf("Total frequent items: %d\n", total_frequent_items);
    printf("Total frequent pairs: %d\n", metrics.total_frequent_pairs);
    printf("Total frequent itemsets: %d\n", metrics.total_frequent_itemsets);
    printf("Total association rules: %d\n", metrics.total_rules);
    printf("#################### \n");
    
    // Write metrics to CSV file for further analysis
    FILE *metrics_file = fopen("performance_metrics.csv", "w");
    fprintf(metrics_file, "Metric,Value\n");
    fprintf(metrics_file, "Total execution time (s),%.3f\n", metrics.total_time);
    fprintf(metrics_file, "File splitting time (s),%.3f\n", metrics.file_split_time);
    fprintf(metrics_file, "Frequent items mining time (s),%.3f\n", metrics.frequent_items_time);
    fprintf(metrics_file, "Triangle matrix build time (s),%.3f\n", metrics.triangle_matrix_build_time);
    fprintf(metrics_file, "Pair generation time (s),%.3f\n", metrics.pair_generation_time);
    fprintf(metrics_file, "Triangle pruning time (s),%.3f\n", metrics.triangle_prune_time);
    fprintf(metrics_file, "Large itemset mining time (s),%.3f\n", metrics.large_itemset_time);
    fprintf(metrics_file, "Rule generation time (s),%.3f\n", metrics.rule_generation_time);
    fprintf(metrics_file, "Total transactions,%d\n", metrics.total_transactions);
    fprintf(metrics_file, "Total frequent items,%d\n", total_frequent_items);
    fprintf(metrics_file, "Total frequent pairs,%d\n", metrics.total_frequent_pairs);
    fprintf(metrics_file, "Total frequent itemsets,%d\n", metrics.total_frequent_itemsets);
    fprintf(metrics_file, "Total association rules,%d\n", metrics.total_rules);
    fclose(metrics_file);
  }
  
  MPI_Finalize();
  return 0;
}