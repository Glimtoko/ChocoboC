#CC = clang-10
#CFLAGS = -g -std=c17 -ffp-exception-behavior=strict

CC = gcc
CFLAGS = -g -DDEBUG

OBJECTS = main.o hydro.o mesh.o initial_conditions.o output.o

.PHONY: clean

chocoboc: $(OBJECTS)
	$(CC) *.o -o chocoboc -lm

clean:
	-rm *.o
	-rm chocoboc

main.o: hydro.o mesh.o initial_conditions.o output.o
