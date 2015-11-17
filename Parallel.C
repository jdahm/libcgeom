#include "Parallel.H"

void initParallel(int *pargc, char **pargv[]) {
        MPI_CALL(MPI_Init(pargc, pargv));
}


void finalizeParallel() { MPI_CALL(MPI_Finalize()); }

unsigned int rank() {
        int rank;
        MPI_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
        return rank;
}


unsigned int numProcs() {
        int p;
        MPI_CALL(MPI_Comm_size(MPI_COMM_WORLD, &p));
        return p;
}
