#include "par/status.hpp"

namespace par
{

status::status()
{
        MPI_Status::MPI_SOURCE=MPI_ANY_SOURCE;
        MPI_Status::MPI_TAG=MPI_ANY_TAG;
        MPI_Status::MPI_ERROR=MPI_SUCCESS;
}

int status::source() const { return MPI_Status::MPI_SOURCE; }

int status::tag() const { return MPI_Status::MPI_TAG; }

int status::error() const { return MPI_Status::MPI_ERROR; }

bool status::is_cancelled() const {
        int result;
        MPI_Test_cancelled(reinterpret_cast<const MPI_Status *>(this), &result);
        return result;
}

bool status::is_canceled() const { return is_cancelled(); }

} // namespace par
