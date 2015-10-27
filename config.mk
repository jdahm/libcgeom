MPI_INC=/usr/include/openmpi-x86_64
MPI_LINK=/usr/lib64/openmpi/lib
CPPC=g++
MPICPPC=mpic++
CFLAGS=-std=c++14 -Wall -Wextra -g -O2 -I$(MPI_INC)
LFLAGS=-L$(MPI_LINK)
