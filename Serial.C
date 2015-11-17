#include "Parallel.H"

void initParallel(int *pargc, char **pargv[]) { UNUSED(pargc); UNUSED(pargv); }

void finalizeParallel() { }

unsigned int rank() { return 0; }

unsigned int numProcs() { return 1; }
