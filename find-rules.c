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
  init_table(&local_table);
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

  // Merge frequent ints and their supports into a single HashTable on master
  if (rank == 0) {

    // Keep local frequency table available when creating master table
    HashTable global_table;
    clone_table(&global_table, &local_table);

    printf("#################### \n");
    printf("Total of %i transactions found \n", total_transactions);
    for (int i = 1; i < size; i++) {
      // Get number of unique ints found by each worker
      int num_pairs;
      MPI_Recv(&num_pairs, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      Pair *tpairs = malloc(num_pairs * sizeof(Pair));
      MPI_Recv(tpairs, num_pairs * sizeof(Pair), MPI_BYTE, i, 1, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);

      // Add ints + counts to hashtable on proc 0
      for (int j = 0; j < num_pairs; j++) {
        merge_exact(&global_table, tpairs[j].key, tpairs[j].count);
      }
      printf("Master has recieved %i pairs from proc %i \n", num_pairs, i);
      free(tpairs);
    }
    printf("#################### \n");
    printf("Master done merging \n");
    printf("Hashtable contains %i frequent integers from the file \n",
           global_table.count);
    print_sample(&global_table, 20, rank);
    printf("#################### \n");
    printf("Input processing complete");

  } else {
    // Workers convert their tables into KVpairs, send to master
    int local_items;
    Pair *pairs = convert_table(&local_table, &local_items);
    MPI_Send(&local_table.count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(pairs, local_table.count * sizeof(Pair), MPI_BYTE, 0, 1, MPI_COMM_WORLD);
    free(pairs); // Pairs only used for communication, values available in local
                 // table still
  }


  /* At this point, each proc contains a hashtable of all size 1 frequent
   * itemsets in its chunk of the file. Master contains a hashtable of all size
   * 1 frequent itemsets globally. Next we need to have each proc find all
   * frequent itemsets of all size that exist locally in its chunk.
   *
   *
   *
   */
  MPI_Finalize();
}
// How to do initial input processing?
//
// Main reads every line of input file
//  Keeps running count of # of transactions
// Round robin distribution > Send each line to a different proc
//  Each proc recieves line > add all items not yet seen to local list of
//  items, keep running count
//   Use 2 arrays, one for items one for counts. A[i] has count B[i]?
// At end all procs send local list + count for master to merge
//          vs
// Master assigns procs chunks of input file > in order to do so evenly, would
// need to go line by line anyways to get linecount Each proc works through
// it's chunk, creates local list of unique ints and keeps running count #0(n
// + (n^2)/p)#
//  Use same 2 array strategy to keep counts
//          vs
// Master loops through entire file once, creates list of unique ints and
// keeps counts   #0(n^2)#
