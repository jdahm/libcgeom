#include "Parallel.hh"

int rank()
{
  int rank;
  MPI_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
  return rank;
}

int nProc()
{
  int p;
  MPI_CALL(MPI_Comm_size(MPI_COMM_WORLD, &p));
  return p;
}
