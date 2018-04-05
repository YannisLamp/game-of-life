CC=mpicc 
LIBS = -lm
OBJECTS = game_of_life *.o

all : game_of_life

game_of_life : mpi_main.o mpi_functions.o
	$(CC) -o game_of_life mpi_main.o mpi_functions.o $(LIBS)

mpi_main.o : mpi_main.c
	$(CC) -c mpi_main.c

mpi_functions.o : mpi_functions.c
	$(CC) -c mpi_functions.c

.PHONY: clean
clean :
	rm $(OBJECTS)
