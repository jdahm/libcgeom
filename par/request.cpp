#include <mpi.h>

#include "par/request.hpp"

namespace par
{

request::request() : req(MPI_REQUEST_NULL) { }

request::request(const request &other) : req(other.req) { }

request::request(MPI_Request req) : req(req) { }

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


request_pool::~request_pool() {
	for (std::vector<MPI_Request>::iterator i(reqs.begin()), i_end(reqs.end()); i!=i_end; ++i)
                if ((*i)!=MPI_REQUEST_NULL)
                        MPI_Request_free(&(*i));
}


request_pool::size_type request_pool::size() const {
	return reqs.size();
}

bool request_pool::empty() const {
	return reqs.empty();
}

const status & request_pool::get_status(request_pool::size_type i) const {
	return stats[i];
}

void request_pool::cancel(request_pool::size_type i) {
	MPI_Cancel(&reqs[i]);
}

void request_pool::cancelall() {
	for (size_type i=0; i<reqs.size(); ++i)
                cancel(i);
}

void request_pool::push(const request& other) {
	reqs.push_back(other.req);
	stats.push_back(status());
}

std::pair<bool, request_pool::size_type> request_pool::waitany() {
	int index;
	status s;
	MPI_Waitany(size(), &reqs[0], &index, reinterpret_cast<MPI_Status *>(&s));
	if (index!=MPI_UNDEFINED) {
                stats[index]=s;
                return std::make_pair(true, static_cast<size_type>(index));
	}
	return std::make_pair(false, size());
}

std::pair<bool, request_pool::size_type> request_pool::testany() {
	int index, flag;
	status s;
	MPI_Testany(size(), &reqs[0], &index, &flag, reinterpret_cast<MPI_Status *>(&s));
	if (flag and index!=MPI_UNDEFINED) {
                stats[index]=s;
                return std::make_pair(true, static_cast<size_type>(index));
	}
	return std::make_pair(static_cast<bool>(flag), size());
}

void request_pool::waitall() {
	MPI_Waitall(size(), reqs.data(), reinterpret_cast<MPI_Status *>(stats.data()));
}

bool request_pool::testall() {
	int flag;
	MPI_Testall(size(), reqs.data(), &flag, reinterpret_cast<MPI_Status *>(stats.data()));
	return flag;
}
} // namespace par
