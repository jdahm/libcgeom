#ifndef PAR_COMMUNICATOR_HPP
#define PAR_COMMUNICATOR_HPP

#include <mpi.h>

#include "par/status.hpp"
#include "par/request.hpp"
#include "par/traits.hpp"

namespace par
{

// Forward declare env
namespace detail { class env; }

class communicator {
protected:
        MPI_Comm comm;
        communicator(MPI_Comm comm=MPI_COMM_NULL);
public:
        enum class equality_type { ident=MPI_IDENT,
                        congruent=MPI_CONGRUENT,
                        similar=MPI_SIMILAR,
                        unequal=MPI_UNEQUAL };

        friend class detail::env;

        communicator(const communicator&);
        communicator(communicator&&);
        communicator(const communicator&, int, int=0);
        ~communicator();

        communicator& operator=(const communicator&) = delete;
        bool operator==(const communicator&) const;
        bool operator!=(const communicator&) const;

        equality_type compare(const communicator&) const;

        int size() const;
        int rank() const;

        void abort(const std::string &msg, int err) const;

        status probe(int source, int tag=0) const;

        void barrier() const;

        MPI_Comm raw() const { return comm; }

        template<typename T>
        void bcast(int root, T* buf, int count = 1) const {
                MPI_Bcast(buf,
                          count,
                          mpi_datatype<T>::type,
                          root,
                          comm);
        }

        template<typename O, typename T>
        void scan(const T* sendbuf, T* recvbuf, O op, int count = 1) const {
                MPI_Scan(sendbuf,
                         recvbuf,
                         count,
                         mpi_datatype<T>::type,
                         mpi_op<O>::type,
                         comm);
                static_cast<void>(op);
        }

        template<typename O, typename T>
        void reduce(const T* sendbuf, T* recvbuf, O op, int count = 1) const {
                MPI_Reduce(&sendbuf,
                           &recvbuf,
                           count,
                           mpi_datatype<T>::type,
                           mpi_op<O>::type,
                           comm);
                static_cast<void>(op);
        }

        template<typename O, typename T>
        void allreduce(const T* sendbuf, T* recvbuf, O op, int count = 1) const {
                MPI_Allreduce(sendbuf,
                              recvbuf,
                              count,
                              mpi_datatype<T>::type,
                              mpi_op<O>::type,
                              comm);
                static_cast<void>(op);
        }

        template<typename T>
        void allgather(const T* sendbuf, int sendcount, T* recvbuf, int recvcount) const {
                MPI_Allgather(sendbuf,
                              sendcount,
                              mpi_datatype<T>::type,
                              recvbuf,
                              recvcount,
                              mpi_datatype<T>::type,
                              comm);
        }

        template<typename T>
        void allgatherv(const T* sendbuf, int sendcount, T* recvbuf,
                        const int* recvcounts, const int* displs) const {
                MPI_Allgatherv(sendbuf,
                               sendcount,
                               mpi_datatype<T>::type,
                               recvbuf,
                               recvcounts,
                               displs,
                               mpi_datatype<T>::type,
                               comm);
        }


        template<typename T>
        void alltoall(const T* sendbuf, unsigned int sendcount,
                      T* recvbuf, unsigned int recvcount) const {
                MPI_Alltoall(sendbuf,
                             sendcount,
                             mpi_datatype<T>::type,
                             recvbuf,
                             recvcount,
                             mpi_datatype<T>::type,
                             comm);
        }

        template<typename T>
        void alltoallv(const T* sendbuf, const int* sendcounts, const int* sdispls,
                       T* recvbuf, const int* recvcounts, const int* rdispls) const {
                MPI_Alltoallv(sendbuf,
                              sendcounts,
                              sdispls,
                              mpi_datatype<T>::type,
                              recvbuf,
                              recvcounts,
                              rdispls,
                              mpi_datatype<T>::type,
                              comm);
        }

        template<typename T>
        request isend(const T* data, int dest, int count = 1, int tag = 0) const {
                MPI_Request req;
                MPI_Isend(data, count, mpi_datatype<T>::type, dest, tag, comm, &req);
                return request(req);
        }

        template<typename T>
        request irecv(T* data, int source, int count = 1, int tag = 0) const {
                MPI_Request req;
                MPI_Irecv(data, count, mpi_datatype<T>::type, source, tag, comm, &req);
                return request(req);
        }

};

} // namespace par

#endif
