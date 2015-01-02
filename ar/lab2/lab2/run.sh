mpicc -std=c99 -Wall main.c -o prog -lm && mpiexec -n $1 ./prog $2
