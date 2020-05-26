CC = clang
CFLAGS = -g

OBJECTS = main.o hydro.o mesh.o

.PHONY: clean

chocoboc: $(OBJECTS)
	$(CC) *.o -o chocoboc -lm

clean:
	-rm *.o
	-rm chocoboc

main.o: hydro.o mesh.o
