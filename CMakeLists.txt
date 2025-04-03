cmake_minimum_required(VERSION 3.10)
project(ParallelARM C)

# If you need MPI:
find_package(MPI REQUIRED)

# Add includes—assuming headers are in the same folder
include_directories(${MPI_INCLUDE_PATH})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# List the source files
set(SOURCE_FILES
    Apriori.c
    dynamic_hash_table.c
    find-rules.c
    io-processing.c
    itemset_hash_table.c
    map_items.c
    map_lookup.c
    rules-aux.c
)

# Create an executable named find-rules
add_executable(find-rules ${SOURCE_FILES})

# Link MPI (if required)
target_link_libraries(find-rules ${MPI_LIBRARIES})

