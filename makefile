MPICC       ?= mpicc
CFLAGS      ?= -O2 -Wall -Wextra
TARGET      = find-rules

SOURCES     = Apriori.c \
              dynamic_hash_table.c \
              find-rules.c \
              io-processing.c \
              itemset_hash_table.c \
              map_items.c \
              map_lookup.c \
              rules-aux.c
OBJECTS     = $(SOURCES:.c=.o)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(MPICC) $(CFLAGS) -o $@ $(OBJECTS)

%.o: %.c
	$(MPICC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

