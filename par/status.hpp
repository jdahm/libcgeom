#ifndef PAR_STATUS_HPP
#define PAR_STATUS_HPP

#include <mpi.h>

#include "par/traits.hpp"

namespace par
{

class status : private MPI_Status {
public:
        status();
        int source() const;
        int tag() const;
        int error() const;
        bool is_cancelled() const;
        bool is_canceled() const;

        template<typename T>
        int get_count() const {
                int result;
                MPI_Get_count(
                        reinterpret_cast<const MPI_Status *>(this),
                        mpi_datatype<T>::type,
                        &result);
                return result;
        }
};

} // namespace par

#endif
