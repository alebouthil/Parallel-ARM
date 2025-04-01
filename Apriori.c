#include "Apriori.h"
#include "dynamic_hash_table.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_ITEM_ID 20000

static inline int get_triangle_index(int i, int j, int n) {
    if (i == j) return -1; // Diagonal isn't stored
    if (i > j) { int tmp = i; i = j; j = tmp; }
//   (((n) * (n - 1) - (n - (i)) * (n - (i) + 1)) / 2 + ((j) - (i)))
    return (n * (n - 1) - (n - i) * (n - i + 1)) / 2 + (j - i - 1);
}

TriangularMatrix *build_tri_matrix(HashTable *frequent_items) {
  if (frequent_items == NULL || frequent_items->count <= 0) {
    return NULL;
  }

  int n = frequent_items->count;

  // Step 1: Build item-to-index map
  int *item_ids = malloc(frequent_items->count * sizeof(int));
  if (item_ids == NULL) {
    printf("Memory allocation failed for item_ids\n");
    return NULL;
  }

  sort_ids(frequent_items, item_ids);
  int *item_to_index = malloc(MAX_ITEM_ID * sizeof(int));
  if (item_to_index == NULL) {
    printf("Memory allocation failed for item_to_index\n");
    free(item_ids);
    return NULL;
  }

  int *index_to_item = malloc(n * sizeof(int));
  if (index_to_item == NULL) {
    printf("Memory allocation failed for index_to_item\n");
    free(item_ids);
    free(item_to_index);
    return NULL;
  }

  for (int i = 0; i < MAX_ITEM_ID; i++)
    item_to_index[i] = -1;
  for (int i = 0; i < n; i++) {
    item_to_index[item_ids[i]] = i;
    index_to_item[i] = item_ids[i];
  }

  // Step 2: Allocate matrix
  int total_pairs = (n * (n - 1)) / 2;
  int *matrix = calloc(total_pairs, sizeof(int));
  if (matrix == NULL) {
    printf("Memory allocation failed for triangular matrix\n");
    free(item_ids);
    free(item_to_index);
    free(index_to_item);
    return NULL;
  }

  // Step 3: Return wrapped matrix object
  TriangularMatrix *tm = malloc(sizeof(TriangularMatrix));
  if (tm == NULL) {
    printf("Memory allocation failed for triangular matrix structure\n");
    free(item_ids);
    free(item_to_index);
    free(index_to_item);
    free(matrix);
    return NULL;
  }

  tm->matrix = matrix;
  tm->num_items = n;
  tm->item_to_index_map = item_to_index;
  tm->index_to_item_map = index_to_item;

  free(item_ids); // Free temp array after we're done with it
  return tm;
}

int *extract_frequent(HashTable *frequent_items, char *basket, int *outcount) {
  // Extract an array of all the frequent items from a basket
  if (frequent_items == NULL || basket == NULL || outcount == NULL) {
    *outcount = 0;
    return NULL;
  }

  char *token;
  int *items = malloc(100 * sizeof(int)); // Assume max basket size of 100 items
  if (items == NULL) {
    printf("Memory allocation failed for items\n");
    *outcount = 0;
    return NULL;
  }

  int x = 0;

  token = strtok(basket, " \t\r\n");
  while (token != NULL && x < 100) {
    char *endptr;
    long num = strtol(token, &endptr, 10); // Convert string read into a number
    if (*endptr == '\0') {
      if (get_count(frequent_items, (int)num) > 0) {
        items[x] = (int)num;
        x++;
      }
    }
    token = strtok(NULL, " \t\r\n");
  }
  *outcount = x;

  // If no items were found, free the array and return NULL
  if (x == 0) {
    free(items);
    return NULL;
  }

  return items;
}

ItemSet *prune_triangle(TriangularMatrix *matrix, float support, int baskets,
                        int *out_count) {
  // Prune a triangular matrix, returning only supported pairs
  if (matrix == NULL || out_count == NULL || baskets <= 0) {
    *out_count = 0;
    return NULL;
  }

  int *index_to_item = matrix->index_to_item_map;
  int n = matrix->num_items;
  int total_pairs = (n * (n - 1)) / 2;

  ItemSet *results = malloc(total_pairs * sizeof(ItemSet));
  if (results == NULL) {
    printf("Memory allocation failed for results\n");
    *out_count = 0;
    return NULL;
  }

  int count = 0;
  int *triangle = matrix->matrix;

  for (int i = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {
      int idx = get_triangle_index(i, j, n);

      assert(idx >= 0 && idx < (n * (n - 1)) / 2);

      float item_support = (float)triangle[idx] / baskets;
      if (item_support >= support) {
        ItemSet is;
        is.elements = malloc(2 * sizeof(int));
        if (is.elements == NULL) {
          printf("Memory allocation failed for itemset elements\n");
          // Free any previously allocated elements
          for (int k = 0; k < count; k++) {
            free(results[k].elements);
          }
          free(results);
          *out_count = 0;
          return NULL;
        }
        is.elements[0] = index_to_item[i];
        is.elements[1] = index_to_item[j];
        is.size = 2;
        is.support = item_support;
        is.count = triangle[idx];

        results[count++] = is;
      }
    }
  }

  *out_count = count;

  // If no pairs met the support threshold, free the array and return NULL
  if (count == 0) {
    free(results);
    return NULL;
  }

  // Resize the results array to save memory
  ItemSet *resized_results = realloc(results, count * sizeof(ItemSet));
  if (resized_results != NULL) {
    return resized_results;
  }

  // If realloc fails, just return the original array
  return results;
}

/**
 * Generate pairs of frequent items
 * @param TriangularMatrix *matrix a triangular matrix to store candidate pairs
 * in
 * @param int *items array of frequent ints from a basket
 * @param int count number of frequent ints in the basket
 * @return number of pairs found in basket. pairs stored in TriangularMatrix
 */
int check_pairs(TriangularMatrix *matrix, int *items, int count) {
  // Generate all pairs from an array of ints
  // Increment index in TriangularMatrix
  // items must be frequent
  if (matrix == NULL || items == NULL || count < 2) {
    return 0;
  }

  int pairs = 0;
  int *item_to_index = matrix->item_to_index_map;
  int num_items = matrix->num_items;
  int *triangle = matrix->matrix;

  for (int i = 0; i < count - 1; i++) {
    for (int j = i + 1; j < count; j++) {
      int item_i = items[i];
      int item_j = items[j];

      // Map items to dense indices
      int idx_i =
          (item_i >= 0 && item_i < MAX_ITEM_ID) ? item_to_index[item_i] : -1;
      int idx_j =
          (item_j >= 0 && item_j < MAX_ITEM_ID) ? item_to_index[item_j] : -1;

      if (idx_i == -1 || idx_j == -1)
        continue; // should never happen, but safe check

      // Ensure i < j for triangular matrix
      if (idx_i < idx_j) {
        int index = get_triangle_index(idx_i, idx_j, num_items);

        assert(index >= 0 && index < (num_items * (num_items - 1)) / 2);

        if (index >= 0 && index < (num_items * (num_items - 1) / 2)) {
          triangle[index]++;
          pairs++;
        }
      }
    }
  }
  return pairs;
}

// Helper function to check if item is in an array
int contains_item(int item, int *array, int size) {
  if (array == NULL || size <= 0)
    return 0;

  for (int i = 0; i < size; i++) {
    if (array[i] == item) {
      return 1;
    }
  }
  return 0;
}

// Find a frequent itemset in a list of itemsets
float find_itemset_support(ItemSet *itemsets, int count, int *items, int size) {
  if (itemsets == NULL || count <= 0 || items == NULL || size <= 0) {
    return 0.0;
  }

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
  return 0.0; // Not found
}

// Generate all possible k-sized subsets of an n-element array
void generate_subsets(int *set, int set_size, int subset_size, int *temp,
                      int temp_index, int start_pos, int ***results,
                      int *result_count) {
  if (set == NULL || temp == NULL || results == NULL || result_count == NULL) {
    return;
  }

  // Base case: subset is complete
  if (temp_index == subset_size) {
    // Add current subset to results
    (*results)[*result_count] = (int *)malloc(subset_size * sizeof(int));
    if ((*results)[*result_count] == NULL) {
      printf("Memory allocation failed in generate_subsets\n");
      return;
    }
    memcpy((*results)[*result_count], temp, subset_size * sizeof(int));
    (*result_count)++;
    return;
  }

  // Try all possible elements for current position
  for (int i = start_pos; i <= set_size - (subset_size - temp_index); i++) {
    temp[temp_index] = set[i];
    generate_subsets(set, set_size, subset_size, temp, temp_index + 1, i + 1,
                     results, result_count);
  }
}

// Generate all possible subsets of an array
int **get_all_subsets(int *set, int set_size, int subset_size, int *count) {
  if (set == NULL || count == NULL || subset_size <= 0 || set_size <= 0) {
    *count = 0;
    return NULL;
  }

  // Maximum number of subsets: C(n,k) = n!/(k!(n-k)!)
  // This is an upper bound, we'll reallocate later
  int max_count = 1;
  for (int i = 0; i < subset_size; i++) {
    max_count *= (set_size - i);
    max_count /= (i + 1);
  }

  int **subsets = (int **)malloc(max_count * sizeof(int *));
  if (subsets == NULL) {
    printf("Memory allocation failed in get_all_subsets\n");
    *count = 0;
    return NULL;
  }

  int *temp = (int *)malloc(subset_size * sizeof(int));
  if (temp == NULL) {
    printf("Memory allocation failed for temp in get_all_subsets\n");
    free(subsets);
    *count = 0;
    return NULL;
  }

  *count = 0;

  generate_subsets(set, set_size, subset_size, temp, 0, 0, &subsets, count);

  free(temp);

  // If no subsets were generated, free the array and return NULL
  if (*count == 0) {
    free(subsets);
    return NULL;
  }

  // Resize the array to save memory
  int **resized_subsets = realloc(subsets, (*count) * sizeof(int *));
  if (resized_subsets != NULL) {
    return resized_subsets;
  }

  // If realloc fails, just return the original array
  return subsets;
}

// Create the complement of a subset within the original set
int *get_complement(int *set, int set_size, int *subset, int subset_size,
                    int *result_size) {
  if (set == NULL || subset == NULL || result_size == NULL || set_size <= 0 ||
      subset_size <= 0 || subset_size >= set_size) {
    *result_size = 0;
    return NULL;
  }

  *result_size = set_size - subset_size;
  int *complement = (int *)malloc(*result_size * sizeof(int));
  if (complement == NULL) {
    printf("Memory allocation failed in get_complement\n");
    *result_size = 0;
    return NULL;
  }

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
AssociationRule *generate_rules(ItemSet *itemsets, int itemset_count,
                                float min_confidence, int *rule_count) {
  // Initial allocation for rules
  int max_rules = 1000; // We'll reallocate as needed
  AssociationRule *rules =
      (AssociationRule *)malloc(max_rules * sizeof(AssociationRule));
  if (rules == NULL) {
    printf("Memory allocation failed in generate_rules\n");
    *rule_count = 0;
    return NULL;
  }

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
    for (int antecedent_size = 1; antecedent_size < itemset.size;
         antecedent_size++) {
      // Generate all possible antecedents of this size
      int subset_count;
      int **antecedents = get_all_subsets(itemset.elements, itemset.size,
                                          antecedent_size, &subset_count);
      if (antecedents == NULL)
        continue;

      // For each antecedent, create a rule
      for (int j = 0; j < subset_count; j++) {
        int *antecedent = antecedents[j];
        if (antecedent == NULL)
          continue;

        // Find the complement of the antecedent (the consequent)
        int consequent_size;
        int *consequent =
            get_complement(itemset.elements, itemset.size, antecedent,
                           antecedent_size, &consequent_size);
        if (consequent == NULL)
          continue;

        // Find the support of the antecedent
        float antecedent_support = find_itemset_support(
            itemsets, itemset_count, antecedent, antecedent_size);

        // Calculate confidence = support(X ∪ Y) / support(X)
        // Protection against division by zero
        float confidence = 0.0;
        if (antecedent_support > 0.0001) {
          confidence = full_support / antecedent_support;

          // Clamp confidence to reasonable values
          if (confidence > 1.0)
            confidence = 1.0;
        }

        // If confidence meets or exceeds threshold, add to rules
        if (confidence >= min_confidence) {
          // Resize if needed
          if (*rule_count >= max_rules) {
            max_rules *= 2;
            AssociationRule *temp = (AssociationRule *)realloc(
                rules, max_rules * sizeof(AssociationRule));
            if (temp == NULL) {
              printf("Memory reallocation failed in generate_rules\n");
              // Clean up before returning
              for (int k = 0; k < *rule_count; k++) {
                if (rules[k].antecedent != NULL)
                  free(rules[k].antecedent);
                if (rules[k].consequent != NULL)
                  free(rules[k].consequent);
              }
              free(rules);
              free(consequent);
              for (int k = j; k < subset_count; k++) {
                if (antecedents[k] != NULL)
                  free(antecedents[k]);
              }
              free(antecedents);
              *rule_count = 0;
              return NULL;
            }
            rules = temp;
          }

          // Add the rule
          rules[*rule_count].antecedent =
              (int *)malloc(antecedent_size * sizeof(int));
          if (rules[*rule_count].antecedent == NULL) {
            printf(
                "Memory allocation failed for antecedent in generate_rules\n");
            free(consequent);
            continue;
          }
          memcpy(rules[*rule_count].antecedent, antecedent,
                 antecedent_size * sizeof(int));
          rules[*rule_count].antecedent_size = antecedent_size;

          rules[*rule_count].consequent =
              (int *)malloc(consequent_size * sizeof(int));
          if (rules[*rule_count].consequent == NULL) {
            printf(
                "Memory allocation failed for consequent in generate_rules\n");
            free(rules[*rule_count].antecedent);
            free(consequent);
            continue;
          }
          memcpy(rules[*rule_count].consequent, consequent,
                 consequent_size * sizeof(int));
          rules[*rule_count].consequent_size = consequent_size;

          rules[*rule_count].confidence = confidence;
          rules[*rule_count].support = full_support;

          (*rule_count)++;
        }

        free(consequent);
      }

      // Free memory for antecedents
      for (int j = 0; j < subset_count; j++) {
        if (antecedents[j] != NULL)
          free(antecedents[j]);
      }
      free(antecedents);
    }
  }

  // Trim the rules array to the actual size
  if (*rule_count > 0) {
    AssociationRule *temp = (AssociationRule *)realloc(
        rules, (*rule_count) * sizeof(AssociationRule));
    if (temp != NULL) {
      rules = temp;
    }
  } else if (rules != NULL) {
    free(rules);
    rules = NULL;
  }

  return rules;
}

// Free memory allocated for association rules
void free_rules(AssociationRule *rules, int rule_count) {
  if (rules == NULL)
    return;

  for (int i = 0; i < rule_count; i++) {
    if (rules[i].antecedent != NULL) {
      free(rules[i].antecedent);
      rules[i].antecedent = NULL;
    }
    if (rules[i].consequent != NULL) {
      free(rules[i].consequent);
      rules[i].consequent = NULL;
    }
  }
  free(rules);
}
