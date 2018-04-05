compile:mpicc -fopenmp -o game_of_life ompi_main.c ompi_functions.c -lm
run:mpiexec -np <N> ./game_of_life
