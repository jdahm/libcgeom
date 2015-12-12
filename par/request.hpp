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
        request();

        request(MPI_Request req);

        request(const request& other);

        void cancel();

        std::pair<bool, status> test();

        status wait();

        std::pair<bool, status> get_status();

        ~request();
};

} // namespace par


#endif
