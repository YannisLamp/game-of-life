CC=mpicc 
LIBS = -lm
OBJECTS = game_of_life *.o
FLAGS = -fopenmp

all : game_of_life

game_of_life : ompi_main.o ompi_functions.o
	$(CC) $(FLAGS) -o game_of_life ompi_main.o ompi_functions.o $(LIBS)

ompi_main.o : ompi_main.c
	$(CC) -c ompi_main.c

ompi_functions.o : ompi_functions.c
	$(CC) -c ompi_functions.c

.PHONY: clean
clean :
	rm $(OBJECTS)