#include "dynamic_hash_table.h"
#define TRIANGULAR_INDEX(i, j, n)                                              \
  (((n) * (n - 1) - (n - (i)) * (n - (i) + 1)) / 2 + ((j) - (i)))

typedef struct {
  int *matrix;            // counts for size-2 itemsets
  int num_items;          // number of frequent items
  int *item_to_index_map; // maps item IDs to dense matrix indices
  int *index_to_item_map; // (optional) inverse lookup
} TriangularMatrix;

typedef struct {
  int *elements;
  int size;
  float support;
} ItemSet;

ItemSet *prune_triangle(TriangularMatrix *matrix, float support, int baskets,
                        int *out_count);
int *extract_frequent(HashTable *frequent_items, char *basket, int *outcount);
TriangularMatrix *build_tri_matrix(HashTable *frequent_items);
void check_pairs(TriangularMatrix *matrix, int *items, int count);
