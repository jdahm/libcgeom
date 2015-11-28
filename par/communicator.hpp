#ifndef PAR_COMMUNICATOR_HPP
#define PAR_COMMUNICATOR_HPP

#include <mpi.h>

#include "par/status.hpp"
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

        communicator(const communicator &other);
        communicator(communicator &&other);
        ~communicator();

        void operator=(const communicator &)=delete;
        bool operator==(const communicator &other) const;
        bool operator!=(const communicator &other) const;

        equality_type compare(const communicator &other) const;

        int size() const;
        int rank() const;

        void abort(const std::string &msg, int err) const;

        status probe(int source, int tag=0) const;

        void barrier() const;

        MPI_Comm raw() const { return comm; }

        template<typename T>
        void bcast(int root, T* data, int count = 1) const {
                MPI_Bcast(data,
                          count,
                          mpi_datatype<T>::type,
                          root,
                          comm);
        }

        template<typename O, typename T>
        void scan(const T* senddata, T* recvdata, O op, int count = 1) const {
                MPI_Scan(senddata,
                         recvdata,
                         count,
                         mpi_datatype<T>::type,
                         mpi_op<O>::type,
                         comm);
                static_cast<void>(op);
        }

        template<typename O, typename T>
        void reduce(const T* senddata, T* recvdata, O op, int count = 1) const {
                MPI_Reduce(&senddata,
                           &recvdata,
                           count,
                           mpi_datatype<T>::type,
                           mpi_op<O>::type,
                           comm);
                static_cast<void>(op);
        }

        template<typename O, typename T>
        void allreduce(const T* senddata, T* recvdata, O op, int count = 1) const {
                MPI_Allreduce(senddata,
                              recvdata,
                              count,
                              mpi_datatype<T>::type,
                              mpi_op<O>::type,
                              comm);
                static_cast<void>(op);
        }

        template<typename T>
        void allgather(const T* senddata, T* recvdata, int count = 1) const {
                MPI_Allgather(senddata,
                              count,
                              mpi_datatype<T>::type,
                              recvdata,
                              count,
                              mpi_datatype<T>::type,
                              comm);
        }

        template<typename T>
        void alltoall(const T* senddata, unsigned int sendcount,
                      T* recvdata, unsigned int recvcount) const {
                MPI_Alltoall(senddata,
                             sendcount,
                             mpi_datatype<T>::type,
                             recvdata,
                             recvcount,
                             mpi_datatype<T>::type,
                             comm);
        }

        template<typename T, typename V>
        void alltoallv(const T* senddata, const V* sendcounts, const V* sdispls,
                       T* recvdata, const V* recvcounts, const V* rdispls) const {
                MPI_Alltoallv(senddata,
                              sendcounts,
                              sdispls,
                              mpi_datatype<T>::type,
                              recvdata,
                              recvcounts,
                              rdispls,
                              mpi_datatype<T>::type,
                              comm);
        }
};

} // namespace par

#endif
