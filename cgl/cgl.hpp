#ifndef CGL_HPP
#define CGL_HPP

// Purpose: A couple useful global definitions. Defining a real here removes a
// lot of single typename object templates, since most objects will need the
// concept of a real type.

#include <limits>
#include "par/traits.hpp"

namespace cgl
{

typedef double real;
static constexpr real real_eps = std::numeric_limits<real>::epsilon();

static constexpr MPI_Datatype mpi_real = par::mpi_datatype<real>::type;

enum class ProcTopology {Line, NestedGrid, RCB, RIB};

}

#endif
