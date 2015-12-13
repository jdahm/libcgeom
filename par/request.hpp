#ifndef PAR_REQUEST_HPP
#define PAR_REQUEST_HPP

#include <vector>
#include <utility>
#include <mpi.h>

#include "par/status.hpp"

namespace par
{

class request_pool;

class request {
protected:
        MPI_Request req;

public:
        friend class request_pool;

        request();

        request(MPI_Request req);

        request(const request& other);

        void cancel();

        std::pair<bool, status> test();

        status wait();

        std::pair<bool, status> get_status();

        ~request();
};

class request_pool {
protected:
        std::vector<MPI_Request> reqs;
        std::vector<status> stats;
public:
        typedef std::vector<MPI_Request>::size_type size_type;

        ~request_pool();

        size_type size() const;

        bool empty() const;

        const status& get_status(size_type) const;

        void cancel(size_type);

        void cancelall();

        void push(const request&);

        std::pair<bool, size_type> waitany();

        std::pair<bool, size_type> testany();

        void waitall();

        bool testall();
};

} // namespace par


#endif
