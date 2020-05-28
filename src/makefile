#CC = clang-10
#CFLAGS = -g -std=c17 -ffp-exception-behavior=strict

CC = clang-10
CFLAGS = -O3

OBJECTS = main.o hydro.o mesh.o initial_conditions.o output.o ini.o input.o

.PHONY: clean

chocoboc: $(OBJECTS)
	$(CC) *.o -o chocoboc -lm

clean:
	-rm *.o
	-rm chocoboc

main.o: hydro.o mesh.o initial_conditions.o output.o input.o

input.o: ini.o
