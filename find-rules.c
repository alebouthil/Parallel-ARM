#include "Apriori.h"
#include "dynamic_hash_table.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "io-processing.c"

int main(int argc, char **argv) {
  int size, rank;
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (argc < 3) {
    if (rank == 0)
      fprintf(stderr, "Usage: %s <inputfile> <support_Threshold> \n", argv[0]);
    MPI_Finalize();
    return 1;
  }

  float global_support = strtof(argv[2], NULL);

  // Split file into roughly equal sized portions for processing
  long *split_points = NULL;
  split_points = (long *)malloc((size - 1) * sizeof(long));
  if (rank == 0) {
    printf("#################### \n");
    printf("Begin input processing \n");
    split_file(argv[1], split_points, size);
    printf("Split points generated \n");
    for (int i = 0; i <= size - 1; i++) {
      printf("split %i is %li \n", i, split_points[i]);
    }
    printf("#################### \n");
  }

  // Send split points to all processors
  MPI_Bcast(split_points, size, MPI_LONG, 0, MPI_COMM_WORLD);

  // Each processor finds frequent ints, support and transaction counts
  HashTable local_table;
  init_table(&local_table, INT_TYPE);
  int local_transaction_count;
  local_transaction_count =
      process_chunk(argv[1], split_points, rank, &local_table, global_support);
  printf("Proc %i has finished processing %i lines in its file chunk \n", rank,
         local_transaction_count);
  free(split_points);

  // Get total number of transactions
  int total_transactions;
  MPI_Reduce(&local_transaction_count, &total_transactions, 1, MPI_INT, MPI_SUM,
             0, MPI_COMM_WORLD);

  if (rank == 0) {
    printf("#################### \n");
    printf("Total of %i transactions found \n", total_transactions);
    //    for (int i = 1; i < size; i++) {
    //      // Get number of unique ints found by each worker
    //      int num_pairs;
    //      MPI_Recv(&num_pairs, 1, MPI_INT, i, 0, MPI_COMM_WORLD,
    //      MPI_STATUS_IGNORE); Pair *tpairs = malloc(num_pairs * sizeof(Pair));
    //      MPI_Recv(tpairs, num_pairs * sizeof(Pair), MPI_BYTE, i, 1,
    //      MPI_COMM_WORLD,
    //               MPI_STATUS_IGNORE);
    //
    //      // Add ints + counts to hashtable on proc 0
    //      for (int j = 0; j < num_pairs; j++) {
    //        merge_exact(&global_table, tpairs[j].key, tpairs[j].value.f);
    //      }
    //      printf("Master has recieved %i pairs from proc %i \n", num_pairs,
    //      i); free(tpairs);
    //    }
    printf("#################### \n");
    printf("Input processing complete \n");
    printf("#################### ");
    printf("Beginning SON by performing Apriori on each processor");
  }

  TriangularMatrix local_tri = *build_tri_matrix(&local_table);

  // Set file buffer for reading
  FILE *fp = fopen(argv[1], "r");
  if (rank != 0) {
    fseek(fp, split_points[rank - 1], SEEK_SET);
  } else {
    fseek(fp, 0, SEEK_SET);
  }

  // Read each line, entering candidate pairs into triangle matrix
  while (ftell(fp) < split_points[rank]) {
    int outcount;
    int *frequent_items = extract_frequent(
        &local_table, get_line(fp, rank, split_points), &outcount);
    if (outcount < 2) {
      continue;
    }
    check_pairs(&local_tri, frequent_items, outcount);
  }

  // Create list of supported pairs from triangle
  int valid_pairs;
  float local_support =
      global_support * ((float)local_transaction_count / total_transactions);
  ItemSet *frequent_pairs;
  frequent_pairs = prune_triangle(&local_tri, local_support,
                                  local_transaction_count, &valid_pairs);
  if (rank == 0) {
    for (int i = 0; i < 10; i++) {
      printf("found pair %i,%i", frequent_pairs[i].elements[i],
             frequent_pairs[i].elements[i]);
    }
  }

  MPI_Finalize();
}
