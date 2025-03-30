#include "dynamic_hash_table.h"
#define TRIANGULAR_INDEX(i, j, n) ((i-n) * (n - (i/2))) + j - i

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

ItemSet *prune_triangle(TriangularMatrix *matrix, float support, int baskets,
                        int *out_count);
int *extract_frequent(HashTable *frequent_items, char *basket, int *outcount);
TriangularMatrix *build_tri_matrix(HashTable *frequent_items);
void check_pairs(TriangularMatrix *matrix, int *items, int count);

// Function to generate association rules from frequent itemsets
AssociationRule *generate_rules(ItemSet *itemsets, int itemset_count, 
                              float min_confidence, int *rule_count);

// Function to free memory allocated for association rules
void free_rules(AssociationRule *rules, int rule_count);
