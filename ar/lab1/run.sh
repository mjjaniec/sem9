mpicc -std=c99 main.c -o lab1 && mpiexec -n $1 ./lab1 "$2"
