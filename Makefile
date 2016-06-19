CC         = gcc
CFLAGS     = -c -lm -Wall -o3
GCFLAGS    =  -lm -Wall -o3
SOURCES    = ./CommonFiles/matrix.c        \
	     ./CommonFiles/linked_list.c   \
	     ./CommonFiles/graph.c         \
	     ./CommonFiles/ilup.c          \
	     ./CommonFiles/rcm.c           \
	     ./GMRES/gmres.c               \
	     ./CommonFiles/Vector/Vector.c \
	     program.c
OBJECTS    = $(SOURCES:.c=.o)
EXECUTABLE = program

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(GCFLAGS) -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

debug: CFLAGS += -DDEBUG -g
debug: GCFLAGS += -DDEBUG -g
debug: all

$(DEBUG_EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(GCFLAGS) -o $@



