compile:mpicc -o game_of_life mpi_main.c mpi_functions.c -lm
run:mpiexec -np <N> ./game_of_life
