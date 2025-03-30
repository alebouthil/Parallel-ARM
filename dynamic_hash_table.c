#include "dynamic_hash_table.h"

unsigned int hash(int key, int size) { return (unsigned int)key % size; }

void init_table(HashTable *table) {
  table->size = INITIAL_TABLE_SIZE;
  table->count = 0;
  table->entries = calloc(table->size, sizeof(HashEntry));
}

void resize_table(HashTable *table) {
  int old_size = table->size;
  HashEntry *old_entries = table->entries;

  table->size = old_size * 2;
  table->entries = calloc(table->size, sizeof(HashEntry));
  table->count = 0;

  for (int i = 0; i < old_size; i++) {
    if (old_entries[i].occupied) {
      insert(table, old_entries[i].key);
      table->entries[hash(old_entries[i].key, table->size)].count =
          old_entries[i].count;
    }
  }

  free(old_entries);
}

void print_sample(const HashTable *table, int max_items, int rank) {
  printf("Rank %d: printing up to %d hash table entries:\n", rank, max_items);
  int printed = 0;

  for (int i = 0; i < table->size && printed < max_items; i++) {
    if (table->entries[i].occupied) {
      printf("  key = %d, count = %d\n", table->entries[i].key,
             table->entries[i].count);
      printed++;
    }
  }

  if (printed == 0) {
    printf("  (No entries)\n");
  }

  fflush(stdout); // Ensure output appears in MPI jobs
}

void insert(HashTable *table, int key) {
  if ((float)table->count / table->size >= LOAD_FACTOR) {
    printf("resizing local table");
    resize_table(table);
  }

  unsigned int index = hash(key, table->size);
  while (table->entries[index].occupied) {
    if (table->entries[index].key == key) {
      table->entries[index].count++;
      return;
    }
    index = (index + 1) % table->size;
  }

  table->entries[index].occupied = 1;
  table->entries[index].key = key;
  table->entries[index].count = 1;
  table->count++;
}

int lookup(HashTable *table, int key) {
  unsigned int index = hash(key, table->size);
  while (table->entries[index].occupied) {
    if (table->entries[index].key == key) {
      return table->entries[index].count;
    }
    index = (index + 1) % table->size;
  }
  return 0;
}

void free_table(HashTable *table) {
  free(table->entries);
  table->entries = NULL;
  table->size = 0;
  table->count = 0;
}

void merge_exact(HashTable *table, int key, int count) {

  if ((float)table->count / table->size >= LOAD_FACTOR) {
    resize_table(table);
  }

  unsigned int index = hash(key, table->size);
  while (table->entries[index].occupied) {
    if (table->entries[index].key == key) {
      table->entries[index].count += count;
      return;
    }
    index = (index + 1) % table->size;
  }

  table->entries[index].occupied = 1;
  table->entries[index].key = key;
  table->entries[index].count = count;
  table->count++;
}

Pair *convert_table(HashTable *table, int *count) {
  // Turns a hashtable into an array of Pair objects
  // Outputs count through pointer to calling code

  Pair *pairs = malloc(sizeof(Pair) * table->count);
  int index = 0;

  for (int i = 0; i < table->size; i++) {
    if (table->entries[i].occupied) {
      pairs[index].key = table->entries[i].key;
      pairs[index].count = table->entries[i].count;
      index++;
    }
  }

  *count = index;
  return pairs;
}

void clone_table(HashTable *dest, const HashTable *src) {
    init_table(dest);
    for (int i = 0; i < src->size; i++) {
        if (src->entries[i].occupied) {
            merge_exact(dest, src->entries[i].key, src->entries[i].count);
        }
    }
}
