#ifndef DYNAMIC_HASH_TABLE_H
#define DYNAMIC_HASH_TABLE_H

#include <stdio.h>
#include <stdlib.h>

#define INITIAL_TABLE_SIZE 65536
#define LOAD_FACTOR 0.75

typedef enum { INT_TYPE, FLOAT_TYPE } ValueType;

typedef union {
  int i;
  float f;
} HashValue;

typedef struct {
  int key;
  HashValue value;
  int occupied; // 1 if used, 0 otherwise
} HashEntry;

typedef struct {
  HashEntry *entries;
  int size;  // Current table size
  int count; // Number of used slots
  ValueType type;
} HashTable;

typedef struct {
  int key;
  int count;
} Pair;

void print_sample(const HashTable *table, int max_items, int rank);
void init_table(HashTable *table, ValueType type);
void insert(HashTable *table, int key);
int lookup(HashTable *table, int key);
void resize_table(HashTable *table);
void free_table(HashTable *table);
void merge_exact(HashTable *dest, int key, float input);
Pair *convert_table(HashTable *src, int *count);
void clone_table(HashTable *dest, const HashTable *src);

#endif
