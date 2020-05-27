CC = clang
CFLAGS = -g -std=c17

OBJECTS = main.o hydro.o mesh.o initial_conditions.o

.PHONY: clean

chocoboc: $(OBJECTS)
	$(CC) *.o -o chocoboc -lm

clean:
	-rm *.o
	-rm chocoboc

main.o: hydro.o mesh.o initial_conditions.o
