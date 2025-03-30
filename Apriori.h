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

// Structure to represent an association rule
typedef struct {
  int *antecedent;      // Items in the "if" part of the rule
  int antecedent_size;  // Number of items in the antecedent
  int *consequent;      // Items in the "then" part of the rule
  int consequent_size;  // Number of items in the consequent
  float confidence;     // Confidence of the rule (conditional probability)
  float support;        // Support of the rule (joint probability)
} AssociationRule;

// Core functions for Apriori algorithm
ItemSet *prune_triangle(TriangularMatrix *matrix, float support, int baskets,
                        int *out_count);
int *extract_frequent(HashTable *frequent_items, char *basket, int *outcount);
TriangularMatrix *build_tri_matrix(HashTable *frequent_items);
void check_pairs(TriangularMatrix *matrix, int *items, int count);

// Helper functions
int contains_item(int item, int *array, int size);
float find_itemset_support(ItemSet *itemsets, int count, int *items, int size);
float find_itemset_support_fixed(ItemSet *itemsets, int count, int *items, int size);
int** get_all_subsets(int *set, int set_size, int subset_size, int *count);
int* get_complement(int *set, int set_size, int *subset, int subset_size, int *result_size);

// Rule generation functions
AssociationRule *generate_rules(ItemSet *itemsets, int itemset_count, 
                              float min_confidence, int *rule_count);
AssociationRule *generate_rules_fixed(ItemSet *itemsets, int itemset_count, 
                                   float min_confidence, int *rule_count);

// Function to free memory allocated for association rules
void free_rules(AssociationRule *rules, int rule_count);