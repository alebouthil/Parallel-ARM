#include "dynamic_hash_table.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#define BUF_SIZE 65536

long get_file_size(const char *filename) {
  // Uses the stat() system call to return the size of filename in bytes
  struct stat st;

  if (stat(filename, &st) == 0) {
    return st.st_size;
  } else {
    perror("Failed to get filesize");
    return -1;
  }
}

void split_file(const char *filename, long split_points[], int p) {
  // Fill array split(char **)NULL_points with byte values
  // Creates split points based on the number of processes

  // Get an estimate of each chunk size
  long fsize = get_file_size(filename);
  int chunk_size = fsize / p;

  // Create split points
  FILE *fp = fopen(filename, "r");
  for (int i = 0; i < p; i++) {
    fseek(fp, chunk_size, SEEK_CUR);

    // Skip to start of next line
    char buffer[1024];
    fgets(buffer, sizeof(buffer), fp);

    // Store byte value of line start
    split_points[i] = ftell(fp);
  }
  split_points[p - 1] = fsize;

  fclose(fp);
  printf("filesize of %li determined \n", fsize);
}

int process_chunk(const char *filename, long split_points[], int rank,
                  HashTable *frequent_table, float global_support) {
  // Process filename from split_points[rank] to split_points[rank + 1]
  // Returns the number of lines processed
  // Fills hashtable with frequent and unique integers as well as their support

  int lines = 0;
  HashTable local_table;
  init_table(&local_table);

  // Cluster file system allows concurrent reads
  FILE *fp = fopen(filename, "r");
  if (rank != 0) {
    fseek(fp, split_points[rank - 1], SEEK_SET);
  }

  printf("Proc %i reading from %li to %li \n", rank, ftell(fp),
         split_points[rank]);
  fflush(stdout);

  // Process file chunk to extract unique integers
  // Fills local_table with all unique integers and their frequency
  char buffer[2048];
  while (ftell(fp) < split_points[rank] &&
         fgets(buffer, sizeof(buffer), fp) != NULL) {
    lines++;

    // Tokenize line for hashing
    char *token;
    token = strtok(buffer, " \t\r\n");

    while (token != NULL) {
      char *endptr;
      long num =
          strtol(token, &endptr, 10); // Convert string read into a number
      if (*endptr == '\0') {
        insert(&local_table, (int)num); // Insert number into hashtable
      }
      token = strtok(NULL, " \t\r\n");
    }
    memset(buffer, 0, sizeof(buffer));
  }

  // Local support threshold is based on the percentage of the file that this
  // proc is responsible for.
  float local_support = global_support * ((split_points[rank] - ftell(fp)) /
                                          get_file_size(filename));
  printf("Proc %i has a local support of %f", rank, local_support);

  // Find support for each unique int, adding it to the output table if high
  // enough
  for (int i = 0; i <= local_table.count; i++) {
    float support = local_table.entries[i].count / lines;
    if (support >= local_support) {
      insert(frequent_table, local_table.entries[i].key);
    }
  }
  return lines; // Return # of transactions processed in local chunk
}
