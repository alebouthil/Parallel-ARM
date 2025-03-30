#include "Apriori.h"
#include "dynamic_hash_table.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define MAX_ITEM_ID 20000

TriangularMatrix *build_tri_matrix(HashTable *frequent_items) {
  int n = frequent_items->count;

  // Step 1: Build item-to-index map
  int* item_ids = malloc(frequent_items->count * sizeof(int));
  sort_ids(frequent_items, item_ids);
  int *item_to_index = malloc(MAX_ITEM_ID * sizeof(int));
  int *index_to_item = malloc(n * sizeof(int));
  for (int i = 0; i < MAX_ITEM_ID; i++)
    item_to_index[i] = -1;
  for (int i = 0; i < n; i++) {
    item_to_index[item_ids[i]] = i;
    index_to_item[i] = item_ids[i];
  }

  // Step 2: Allocate matrix
  int total_pairs = n * (n - 1) / 2;
  int *matrix = calloc(total_pairs, sizeof(int));

  // Step 3: Return wrapped matrix object
  TriangularMatrix *tm = malloc(sizeof(TriangularMatrix));
  tm->matrix = matrix;
  tm->num_items = n;
  tm->item_to_index_map = item_to_index;
  tm->index_to_item_map = index_to_item;

  return tm;
}

int *extract_frequent(HashTable *frequent_items, char *basket, int *outcount) {
  // Extract an array of all the frequent items from a basket
  char *token;
  int *items = malloc(100 * sizeof(int)); // Assume max basket size of 100 items
  int x = 0;

  token = strtok(basket, " ");
  while (token != NULL) {
    if (get_count(frequent_items, strtol(token, NULL, 10)) > 0) {
      items[x] = strtol(token, NULL, 10);
      x++;
    }
    token = strtok(NULL, " ");
  }
  *outcount = x;
  return items;
}

ItemSet *prune_triangle(TriangularMatrix *matrix, float support, int baskets,
                        int *out_count) {
  // Prune a triangular matrix, returning only supported pairs
  int *index_to_item = matrix->index_to_item_map;
  int n = matrix->num_items;
  int total_pairs = n * (n - 1) / 2;

  ItemSet *results = malloc(total_pairs * sizeof(ItemSet));
  int count = 0;

  int *triangle = matrix->matrix;

  int support_threshold_count = (int)ceil(support * baskets);

  for (int i = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {
      int idx = TRIANGULAR_INDEX(i, j, n);
      float support = (float)triangle[idx] / baskets;
      if (support >= support_threshold_count) {
        ItemSet is;
        is.elements = malloc(2 * sizeof(int));
        is.elements[0] = index_to_item[i];
        is.elements[1] = index_to_item[j];
        is.size = 2;
        is.support = support;

        results[count++] = is;
      }
    }
  }

  *out_count = count;
  return results;
}

void check_pairs(TriangularMatrix *matrix, int *items, int count) {
  // Generate all pairs from an array of ints
  // Increment index in TriangularMatrix
  // items must be frequent
  int *item_to_index = matrix->item_to_index_map;
  int num_items = matrix->num_items;
  int *triangle = matrix->matrix;

  for (int i = 0; i < count - 1; i++) {
    for (int j = i + 1; j < count; j++) {
      int item_i = items[i];
      int item_j = items[j];

      // Map items to dense indices
      int idx_i = item_to_index[item_i];
      int idx_j = item_to_index[item_j];

      if (idx_i == -1 || idx_j == -1)
        continue; // should never happen, but safe check

      // Ensure i < j for triangular matrix
      if (idx_i < idx_j) {
        int index = TRIANGULAR_INDEX(idx_i, idx_j, num_items);
        triangle[index]++;
      } else {
        int index = TRIANGULAR_INDEX(idx_j, idx_i, num_items);
        triangle[index]++;
      }
    }
  }
}

// Helper function to check if item is in an array
int contains_item(int item, int *array, int size) {
    for (int i = 0; i < size; i++) {
        if (array[i] == item) {
            return 1;
        }
    }
    return 0;
}

// Find a frequent itemset in a list of itemsets
float find_itemset_support(ItemSet *itemsets, int count, int *items, int size) {
    for (int i = 0; i < count; i++) {
        if (itemsets[i].size == size) {
            // Check if all elements match
            int match_count = 0;
            for (int j = 0; j < size; j++) {
                if (contains_item(items[j], itemsets[i].elements, itemsets[i].size)) {
                    match_count++;
                }
            }
            
            if (match_count == size && size == itemsets[i].size) {
                return itemsets[i].support;
            }
        }
    }
    return 0.0;  // Not found
}

// Generate all possible k-sized subsets of an n-element array
void generate_subsets(int *set, int set_size, int subset_size, int *temp, 
                     int temp_index, int start_pos, int ***results, int *result_count) {
    // Base case: subset is complete
    if (temp_index == subset_size) {
        // Add current subset to results
        (*results)[*result_count] = (int*)malloc(subset_size * sizeof(int));
        memcpy((*results)[*result_count], temp, subset_size * sizeof(int));
        (*result_count)++;
        return;
    }
    
    // Try all possible elements for current position
    for (int i = start_pos; i <= set_size - (subset_size - temp_index); i++) {
        temp[temp_index] = set[i];
        generate_subsets(set, set_size, subset_size, temp, temp_index + 1, i + 1, results, result_count);
    }
}

// Generate all possible subsets of an array
int** get_all_subsets(int *set, int set_size, int subset_size, int *count) {
    // Maximum number of subsets: C(n,k) = n!/(k!(n-k)!)
    // This is an upper bound, we'll reallocate later
    int max_count = 1;
    for (int i = 0; i < subset_size; i++) {
        max_count *= (set_size - i);
        max_count /= (i + 1);
    }
    
    int **subsets = (int**)malloc(max_count * sizeof(int*));
    int *temp = (int*)malloc(subset_size * sizeof(int));
    *count = 0;
    
    generate_subsets(set, set_size, subset_size, temp, 0, 0, &subsets, count);
    
    free(temp);
    return subsets;
}

// Create the complement of a subset within the original set
int* get_complement(int *set, int set_size, int *subset, int subset_size, int *result_size) {
    *result_size = set_size - subset_size;
    int *complement = (int*)malloc(*result_size * sizeof(int));
    int index = 0;
    
    for (int i = 0; i < set_size; i++) {
        int found = 0;
        for (int j = 0; j < subset_size; j++) {
            if (set[i] == subset[j]) {
                found = 1;
                break;
            }
        }
        
        if (!found) {
            complement[index++] = set[i];
        }
    }
    
    return complement;
}

/**
 * Generate association rules from an array of frequent itemsets
 * 
 * @param itemsets Array of frequent itemsets
 * @param itemset_count Number of itemsets in the array
 * @param min_confidence Minimum confidence threshold
 * @param rule_count Pointer to store the number of rules generated
 * @return Array of AssociationRule structures
 */
AssociationRule* generate_rules(ItemSet *itemsets, int itemset_count, 
                              float min_confidence, int *rule_count) {
    // Initial allocation for rules
    int max_rules = 1000;  // We'll reallocate as needed
    AssociationRule *rules = (AssociationRule*)malloc(max_rules * sizeof(AssociationRule));
    *rule_count = 0;
    
    // Process each itemset with size >= 2
    for (int i = 0; i < itemset_count; i++) {
        ItemSet itemset = itemsets[i];
        
        // Skip itemsets of size 1 (can't form rules)
        if (itemset.size < 2) {
            continue;
        }
        
        // Get the support of the full itemset
        float full_support = itemset.support;
        
        // Generate all possible antecedent sizes (from 1 to size-1)
        for (int antecedent_size = 1; antecedent_size < itemset.size; antecedent_size++) {
            // Generate all possible antecedents of this size
            int subset_count;
            int **antecedents = get_all_subsets(itemset.elements, itemset.size, antecedent_size, &subset_count);
            
            // For each antecedent, create a rule
            for (int j = 0; j < subset_count; j++) {
                int *antecedent = antecedents[j];
                
                // Find the complement of the antecedent (the consequent)
                int consequent_size;
                int *consequent = get_complement(itemset.elements, itemset.size, 
                                               antecedent, antecedent_size, &consequent_size);
                
                // Find the support of the antecedent
                float antecedent_support = find_itemset_support(itemsets, itemset_count, 
                                                             antecedent, antecedent_size);
                
                // Calculate confidence = support(X ∪ Y) / support(X)
                float confidence = full_support / antecedent_support;
                
                // If confidence meets or exceeds threshold, add to rules
                if (confidence >= min_confidence) {
                    // Resize if needed
                    if (*rule_count >= max_rules) {
                        max_rules *= 2;
                        rules = (AssociationRule*)realloc(rules, max_rules * sizeof(AssociationRule));
                    }
                    
                    // Add the rule
                    rules[*rule_count].antecedent = (int*)malloc(antecedent_size * sizeof(int));
                    memcpy(rules[*rule_count].antecedent, antecedent, antecedent_size * sizeof(int));
                    rules[*rule_count].antecedent_size = antecedent_size;
                    
                    rules[*rule_count].consequent = (int*)malloc(consequent_size * sizeof(int));
                    memcpy(rules[*rule_count].consequent, consequent, consequent_size * sizeof(int));
                    rules[*rule_count].consequent_size = consequent_size;
                    
                    rules[*rule_count].confidence = confidence;
                    rules[*rule_count].support = full_support;
                    
                    (*rule_count)++;
                }
                
                free(consequent);
            }
            
            // Free memory for antecedents
            for (int j = 0; j < subset_count; j++) {
                free(antecedents[j]);
            }
            free(antecedents);
        }
    }
    
    // Trim the rules array to the actual size
    if (*rule_count > 0) {
        rules = (AssociationRule*)realloc(rules, (*rule_count) * sizeof(AssociationRule));
    }
    
    return rules;
}

// Free memory allocated for association rules
void free_rules(AssociationRule *rules, int rule_count) {
    for (int i = 0; i < rule_count; i++) {
        free(rules[i].antecedent);
        free(rules[i].consequent);
    }
    free(rules);
}