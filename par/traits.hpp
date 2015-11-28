#ifndef PAR_TRAITS_HPP
#define PAR_TRAITS_HPP

#include <mpi.h>

namespace par
{

// mpi_datatype
template<typename T>
struct mpi_datatype { static constexpr MPI_Datatype type = MPI_BYTE; };
template<>
struct mpi_datatype<int> { static constexpr MPI_Datatype type = MPI_INT; };
template<>
struct mpi_datatype<float> { static constexpr MPI_Datatype type = MPI_FLOAT; };
template<>
struct mpi_datatype<double> { static constexpr MPI_Datatype type = MPI_DOUBLE; };
template<>
struct mpi_datatype<unsigned int> { static constexpr MPI_Datatype type = MPI_UNSIGNED; };
template<>
struct mpi_datatype<char *> { static constexpr MPI_Datatype type = MPI_CHAR; };

// mpi_op
struct max { };
struct min { };
struct sum { };
struct prod { };

template<typename O>
struct mpi_op { static constexpr MPI_Op type = MPI_NO_OP; };
template<>
struct mpi_op<max> { static constexpr MPI_Op type = MPI_MAX; };
template<>
struct mpi_op<min> { static constexpr MPI_Op type = MPI_MIN; };
template<>
struct mpi_op<sum> { static constexpr MPI_Op type = MPI_SUM; };
template<>
struct mpi_op<prod> { static constexpr MPI_Op type = MPI_PROD; };

// threading modes
enum class threading_modes { single=MPI_THREAD_SINGLE,
                funneled=MPI_THREAD_FUNNELED, 
                serialized=MPI_THREAD_SERIALIZED,
                multiple=MPI_THREAD_MULTIPLE };

} // namespace par

#endif
