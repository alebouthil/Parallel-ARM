#include "Apriori.h"
#include "dynamic_hash_table.h"
#include <math.h>
#include <string.h>

#define MAX_ITEM_ID 20000

TriangularMatrix *build_triangular_matrix(HashTable *frequent_items) {
  int n = frequent_items->count;

  // Step 1: Build item-to-index map
  int *item_ids;
  sort_ids(frequent_items, item_ids);
  int *item_to_index = malloc(MAX_ITEM_ID * sizeof(int));
  int *index_to_item = malloc(MAX_ITEM_ID * sizeof(int));
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
