#include <mpi.h>

#include "par/request.hpp"

namespace par
{

request::request() : req(MPI_REQUEST_NULL) { }

request::request(const T &other) : req(other.req) { }

void request::cancel()
{
        MPI_Cancel(&req);
}

std::pair<bool, status> request::test()
{
        int result;
        status s;
        MPI_Test(&req, &result, reinterpret_cast<MPI_Status *>(&s));
        return std::make_pair(static_cast<bool>(result), s);
}

status request::wait()
{
        status s;
        MPI_Wait(&req, reinterpret_cast<MPI_Status *>(&s));
        return s;
}

std::pair<bool, status> request::get_status()
{
        int result;
        status s;
        MPI_Request_get_status(req, &result, reinterpret_cast<MPI_Status *>(&s));
        return std::make_pair(static_cast<bool>(result), s);
}

request::~request()
{
        if (req != MPI_REQUEST_NULL) MPI_Request_free(&req);
}

} // namespace par
