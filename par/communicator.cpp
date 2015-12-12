#include "par/status.hpp"
#include "par/communicator.hpp"

namespace par
{

communicator::communicator(MPI_Comm comm) : comm(comm) { }

communicator::communicator(const communicator &other) {
        MPI_Comm_dup(other.comm, &comm);
}

communicator::communicator(communicator &&other) {
        comm=other.comm;
        other.comm=MPI_COMM_NULL;
}

communicator::communicator(const communicator& other, int color, int key) {
        MPI_Comm_split(other.comm, color, key, &comm);
}

communicator::~communicator() {
        if (comm != MPI_COMM_NULL) {
                int result1, result2;
                MPI_Comm_compare(comm, MPI_COMM_WORLD, &result1);
                MPI_Comm_compare(comm, MPI_COMM_SELF, &result2);
                if (result1 != MPI_IDENT && result2 != MPI_IDENT)
                        MPI_Comm_free(&comm);
        }
}

bool communicator::operator==(const communicator &other) const {
        int result;
        MPI_Comm_compare(comm, other.comm, &result);
        return result==MPI_IDENT;
}

bool communicator::operator!=(const communicator &other) const {
        int result;
        MPI_Comm_compare(comm, other.comm, &result);
        return result!=MPI_IDENT;
}

communicator::equality_type communicator::compare(const communicator &other) const {
        int result;
        MPI_Comm_compare(comm, other.comm, &result);
        return static_cast<equality_type>(result);
}

int communicator::size() const {
        int result;
        MPI_Comm_size(comm, &result);
        return result;
}

int communicator::rank() const {
        int result;
        MPI_Comm_rank(comm, &result);
        return result;
}

void communicator::abort(const std::string &msg, int err) const {
        std::cerr << "rank = " << rank() << ": " << msg << std::endl;
        MPI_Abort(comm, err);
}

status communicator::probe(int source, int tag) const {
        status s;
        MPI_Probe(source, tag, comm, reinterpret_cast<MPI_Status *>(&s));
        return s;
}

void communicator::barrier() const {
        MPI_Barrier(comm);
}

} // namespace par
