#ifndef PAR_REQUEST_HPP
#define PAR_REQUEST_HPP

#include <utility>
#include <mpi.h>

#include "par/status.hpp"

namespace par
{

class request {
protected:
        MPI_Request req;

public:
        request() : req(MPI_REQUEST_NULL) { }

        request(MPI_Request req) : req(req) { }

        request(const request& other) : req(other.req) { }

        void cancel() {
                MPI_Cancel(&req);
        }

        std::pair<bool, status> test() {
                int result; 
                status s;
                MPI_Test(&req, &result, reinterpret_cast<MPI_Status *>(&s));
                return std::make_pair(static_cast<bool>(result), s);
        }

        status wait() {
                status s;
                MPI_Wait(&req, reinterpret_cast<MPI_Status *>(&s));
                return s;
        }

        std::pair<bool, status> get_status() {
                int result; 
                status s;
                MPI_Request_get_status(req, &result, reinterpret_cast<MPI_Status *>(&s));
                return std::make_pair(static_cast<bool>(result), s);
        }

        ~request() {
                if (req!=MPI_REQUEST_NULL)
                        MPI_Request_free(&req);
        }
};

} // namespace par


#endif
