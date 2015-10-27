#ifndef PARALLEL_HH
#define PARALLEL_HH

#ifdef USE_PARALLEL
#include <mpi.h>

// Macro for calling MPI while debugging
#ifndef NDEBUG
#define MPI_CALL(X) do{ X; } while (0)
#else
#define MPI_CALL(X) do{ int ierr = X; if (ierr != MPI_SUCCESS) throw std::runtime_error("MPI Error"); } while (0)
#endif
#endif // USE_PARALLEL

int rank();

int nProc();


#endif
