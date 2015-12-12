#include <iostream>
#include <cmath>
#include <vector>
// #include <array>
#include <algorithm>
#include <numeric>
#include <iomanip>

#include "cgl/load_balance.hpp"
#include "par/environment.hpp"
#include "cgl/geom2d.hpp"

namespace cgl
{

typedef Point2d point_type;
typedef std::list<point_type> list_type;

// Forward declare
void lb_zoltan(PSTopology, LBMethod, unsigned int, list_type&);

struct LBData {
        std::vector<unsigned int> size;
        std::vector<real> cut;
        std::vector<int> bin, new_size;
        std::vector<real> send, recv;

        // MPI vectors for alltoall communications to cut down on reallocations
        // These have size num_proc
        std::vector<int> sendcounts, recvcounts, sdispls, rdispls;

        LBData(const par::communicator& comm, unsigned int s) :
                size(par::comm_world().size()), cut(0),
                bin(par::comm_world().size()),
                new_size(par::comm_world().size()),
                send(0), recv(0),
                sendcounts(par::comm_world().size()),
                recvcounts(par::comm_world().size()),
                sdispls(par::comm_world().size()),
                rdispls(par::comm_world().size()) {
                calculate_size(comm, s);
        }

        void change_comm(const par::communicator& comm, unsigned int s) {
                const unsigned int num_proc = comm.size();
                size.resize(num_proc);
                bin.resize(num_proc);
                new_size.resize(num_proc);
                sendcounts.resize(num_proc);
                recvcounts.resize(num_proc);
                sdispls.resize(num_proc);
                rdispls.resize(num_proc);

                calculate_size(comm, s);
        }

private:
        void calculate_size(const par::communicator& comm, unsigned int s) {
                const unsigned int my_rank = comm.rank();
                size[my_rank] = s;
                comm.allgather(size.data() + my_rank, 1, size.data(), 1);
        }
};

static void lb_create_cut(const par::communicator& comm, unsigned int num_cut,
                          unsigned int dimen, const list_type& pl, LBData& lbd)
{
#ifndef NDEBUG
        if (num_cut > pl.size()-1)
                throw std::runtime_error("Cannot create more cuts than points");
#endif
        const unsigned int my_rank = comm.rank();
        // const unsigned int num_proc = comm.size();

        std::vector<real> cut(num_cut);

        const unsigned int num_point = pl.size();
        unsigned int ipoint = 0;
        list_type::const_iterator pit = pl.begin();
        for (unsigned int icut=0; icut<num_cut; icut++) {
                const unsigned int target = std::round((num_point - ipoint) / (num_cut + 1 - icut));
                // std::cout << "icut=" << icut << ", target = " << target << std::endl;
                const unsigned int visited = ipoint;
                while ((ipoint - visited) < target) {
                        ++pit;
                        ipoint++;
                }
                const point_type& p = *pit;
                cut[icut] = p[dimen];
        }

        // allgather number of cuts on each proc
        lbd.recvcounts[my_rank] = cut.size();
        comm.allgather(lbd.recvcounts.data() + my_rank, 1, lbd.recvcounts.data(), 1);

        const unsigned int glob_num_cut =
                std::accumulate(lbd.recvcounts.begin(), lbd.recvcounts.end(), 0);

        // Resize cut (destroying existing cuts)
        lbd.cut.resize(glob_num_cut);

        lbd.rdispls[0] = 0;
        for (unsigned int i=1; i<lbd.rdispls.size(); i++)
                lbd.rdispls[i] = lbd.rdispls[i-1] + lbd.recvcounts[i-1];

        // allgatherv into lbd.pointcut
        comm.allgatherv(cut.data(), cut.size(), lbd.cut.data(),
                        lbd.recvcounts.data(), lbd.rdispls.data());

        // Sort the cuts
        std::sort(lbd.cut.begin(), lbd.cut.end());
}

static void lb_create_bin(const par::communicator& comm, unsigned int dimen,
                          const list_type& pl, LBData& lbd)
{
        const unsigned int dim = point_type::dim;
        const unsigned int my_rank = comm.rank();
        const unsigned int num_proc = comm.size();

        const std::vector<real>::size_type num_cut = lbd.cut.size();

        // Count the number of points in each global bin
        std::vector<int> bin(num_cut + 1);
        unsigned int i = 0;
        for (const point_type& p : pl) {
                while (p[dimen] >= lbd.cut[i] && i < num_cut)
                        i++;
                if (i >= num_cut) break;
                bin[i]++;
        }
        const unsigned int no_last =
                std::accumulate(bin.begin(), --bin.end(), 0);

        const unsigned int num_point = lbd.size[my_rank];

        bin.back() = num_point - no_last;
        // std::cout << "mybin: ";
        // for (const auto& p : bin) std::cout << p << " ";
        // std::cout << std::endl;

        lbd.bin.resize(num_cut + 1);
        comm.allreduce(bin.data(), lbd.bin.data(), par::sum(), lbd.bin.size());

        // std::cout << "point: ";
        // for (const auto& p : pl) std::cout << p << " ";
        // std::cout << std::endl;
        // std::cout << "bin: ";
        // for (const auto& p : lbd.pointbin) std::cout << p << " ";
        // std::cout << std::endl;
        // std::cout << "cut: ";
        // for (const auto& c : lbd.pointcut) std::cout << c << " ";
        // std::cout << std::endl;

        const list_type::size_type glob_num_point =
                std::accumulate(lbd.size.begin(), lbd.size.end(), 0);

        // Bin up pointbins into procbins
        // const real target_avg = static_cast<real>(glob_num_point) / num_proc;
        // std::cout << "target_avg = " << target_avg << std::endl;
        unsigned int point_so_far = 0;
        unsigned int ibin = 0;
        for (unsigned int iproc=0; iproc<num_proc; iproc++) {
                // Loop until adding the next bin would bring us further from
                // the target number of points
                const unsigned int target =
                        std::round((glob_num_point - point_so_far) / (num_proc - iproc));
                // std::cout << "procbin target = " << target << std::endl;
                unsigned int np = 0, my_np = 0;
                while (std::abs(static_cast<real>(np) - target) >=
                       std::abs(static_cast<real>(np + lbd.bin[ibin]) - target)) {
                        np += lbd.bin[ibin];
                        my_np += bin[ibin];
                        ibin++;
                        if (ibin >= num_cut + 1) break;
                }
                // Record sendcount
                lbd.new_size[iproc] = np;
                lbd.sendcounts[iproc] = my_np * dim;
                point_so_far += np;
        }
}

static real lb_calculate_imbalance(LBData& lbd)
{
        const list_type::size_type num_point =
                std::accumulate(lbd.size.begin(), lbd.size.end(), 0);

        const auto minmax_elem = std::minmax_element(lbd.new_size.begin(), lbd.new_size.end());
        const real imbalance = (*minmax_elem.second - *minmax_elem.first) / static_cast<real>(num_point);

        return imbalance;
}

static void lb_migrate_points(const par::communicator& comm, list_type& pl, LBData& lbd)
{
        const unsigned int num_proc = comm.size();
        const unsigned int dim = point_type::dim;

        // Sendcounts should be already filled
        // Perform alltoall for recvcounts from sendcounts
        lbd.recvcounts.resize(num_proc);
        comm.alltoall(lbd.sendcounts.data(), 1, lbd.recvcounts.data(), 1);

        // Resize the send buffer
        lbd.send.resize(pl.size() * dim);

        // Pack up the send buffer
        {
                unsigned int i = 0;
                for (const point_type& p : pl)
                        for (unsigned int d=0; d<dim; d++, i++)
                                lbd.send[i] = p[d];
        }

        // Offsets in memory (for alltoallv)
        lbd.sdispls[0] = 0;
        for (unsigned int i=1; i<num_proc; i++)
                lbd.sdispls[i] = lbd.sdispls[i-1] + lbd.sendcounts[i-1];

        lbd.rdispls[0] = 0;
        for (unsigned int i=1; i<num_proc; i++)
                lbd.rdispls[i] = lbd.rdispls[i-1] + lbd.recvcounts[i-1];

        // Create receive buffer
        const int total_recv = std::accumulate(lbd.recvcounts.begin(), lbd.recvcounts.end(), 0);
        lbd.recv.resize(total_recv);

        // Communicate points (alltoallv)
        comm.alltoallv(lbd.send.data(), lbd.sendcounts.data(), lbd.sdispls.data(),
                       lbd.recv.data(), lbd.recvcounts.data(), lbd.rdispls.data());

        // Clear out current points and reload from the receive buffer
        pl.clear();
        for (unsigned int i=0; i<lbd.recv.size(); i+=dim)
                pl.emplace_back(point_type(lbd.recv[i], lbd.recv[i+1]));
}


static void lb_s2a21d(unsigned int dimen, list_type& pl)
{
        LBData lbd = LBData(par::comm_world(), pl.size());

        const par::communicator& comm = par::comm_world();
        const unsigned int my_rank = comm.rank();
        const unsigned int num_proc = comm.size();

        const unsigned int num_point = lbd.size[my_rank];

        const list_type::size_type glob_num_point =
                std::accumulate(lbd.size.begin(), lbd.size.end(), 0);

        if (glob_num_point < 2 * num_proc) comm.abort("Try more points", 1);

        // First: sort list
        sort_point2d_list(dimen, pl);

        std::vector<unsigned int>::const_iterator max_elem =
                std::max_element(lbd.size.begin(), lbd.size.end());

        unsigned int my_num_cut = std::floor(
                static_cast<real>(num_point) / (*max_elem) *
                std::min(static_cast<int>(num_proc),
                         std::max(static_cast<int>(num_point)-1, 0)));

        // Iterate until enough cuts are placed so the points can be (almost)
        // equally placed onto the processors - or give up
        real imbalance;
        for (unsigned int i=0; i<3; i++) {
                lb_create_cut(comm, my_num_cut, dimen, pl, lbd);
                lb_create_bin(comm, dimen, pl, lbd);
                imbalance = lb_calculate_imbalance(lbd);
                my_num_cut = std::min(2 * my_num_cut, num_point - 1);
                int at_limit = static_cast<int>(my_num_cut == num_point - 1);
                int all_at_limit;
                comm.allreduce(&at_limit, &all_at_limit, par::min(), 1);
                if (imbalance < 0.05 || all_at_limit) break;
        }
        if (my_rank == 0) std::cout << "Imbalance = " << imbalance << std::endl;

        // Migrate the points
        lb_migrate_points(comm, pl, lbd);
}

static void rcb1d_recurse(const par::communicator& comm, unsigned int dimen,
                          list_type& pl, LBData& lbd)
{
        const unsigned int my_rank = comm.rank();
        const unsigned int num_proc = comm.size();

        // First: sort list
        sort_point2d_list(dimen, pl);

        if (num_proc == 1) return;
        lbd.change_comm(comm, pl.size());

        const list_type::size_type orig_num_point = lbd.size[my_rank];
        const unsigned int dim = point_type::dim;

        const list_type::size_type glob_num_point =
                std::accumulate(lbd.size.begin(), lbd.size.end(), 0);

        lbd.cut.resize(num_proc);
        lbd.bin.resize(num_proc);

        const real ratio = static_cast<real>((num_proc+1) / 2) / num_proc;
        list_type::size_type this_offset = ratio * orig_num_point;
        const unsigned int rank_midpt = ratio * num_proc;
        {
                real cut = 0;
                int use = 0;
                if (orig_num_point > 0) {
                        const point_type& p = *std::next(pl.begin(), this_offset);
                        cut = p[dimen];
                        use = 1;
                }
                comm.allgather(&cut, 1, lbd.cut.data(), 1);
                comm.allgather(&use, 1, lbd.bin.data(), 1);
        }

        std::cout << "this_offset = " << this_offset << std::endl;
        std::cout << "ratio = " << ratio << std::endl;
        std::cout << "cuts: ";
        for (auto& p : lbd.cut) std::cout << p << " ";
        std::cout << std::endl;
        std::cout << "use: ";
        for (auto& p : lbd.bin) std::cout << p << " ";
        std::cout << std::endl;

        real midpt = 0.;
        std::cout << "my_rank = " << my_rank << std::endl;
        std::cout << "num_proc = " << num_proc << std::endl;
        for (unsigned int i=0; i<num_proc; i++)
                if (lbd.bin[i] != 0) // include this point
                        midpt += lbd.cut[i] *
                                static_cast<real>(lbd.size[i]) /
                                static_cast<real>(glob_num_point);

        list_type other;
        list_type::const_iterator pit = std::find_if(pl.begin(), pl.end(),
                                                     [midpt,dimen](const point_type& p)
                                                     { return p[dimen] >= midpt; });
        list_type::size_type offset;
        if (my_rank < rank_midpt) {
                other.splice(other.begin(), pl, pit, pl.end());
                offset = pl.size();
        }
        else {
                other.splice(other.begin(), pl, pl.begin(), pit);
                offset = other.size();
        }
        std::cout << "pl: ";
        for (auto& p : pl) std::cout << p << " ";
        std::cout << std::endl;
        std::cout << "other: ";
        for (auto& p : other) std::cout << p << " ";
        std::cout << std::endl;
        std::cout << "offset = " << offset << std::endl;
        lbd.send.resize(other.size() * dim);
        {
                list_type::iterator it=other.begin();
                list_type::size_type i = 0;
                while (it != other.end()) {
                        const point_type& p = *it;
                        for (unsigned int d=0; d<dim; d++, i++)
                                lbd.send[i] = p[d];
                        other.erase(it++);
                }
        }
        std::cout << "send: ";
        for (auto& p : lbd.send) std::cout << p << " ";
        std::cout << std::endl;

        std::cout << "midpt = " << midpt << std::endl;
        std::cout << "rank_midpt = " << rank_midpt << std::endl;
        std::cout << "this_offset = " << this_offset << std::endl;

        // Pack up points
        lbd.sendcounts.resize(num_proc);

        // Zero out
        std::fill(lbd.sendcounts.begin(), lbd.sendcounts.end(), 0);

        if (my_rank < rank_midpt) {
                // Left side
                for (unsigned int i=offset; i<orig_num_point; i++) {
                        const unsigned int inc = (i-offset) % (num_proc - rank_midpt);
                        lbd.sendcounts[rank_midpt + inc] += dim;
                }
        }
        else {
                // Right side
                for (unsigned int i=0; i<offset; i++) {
                        lbd.sendcounts[i % rank_midpt] += dim;
                }
        }

        std::cout << "sendcounts: ";
        for (auto& p : lbd.sendcounts) std::cout << p << " ";
        std::cout << std::endl;

        // alltoall recvcounts
        comm.alltoall(lbd.sendcounts.data(), 1, lbd.recvcounts.data(), 1);

        std::cout << "recvcounts: ";
        for (auto& p : lbd.recvcounts) std::cout << p << " ";
        std::cout << std::endl;

        // Offsets in memory (for alltoallv)
        lbd.sdispls[0] = 0;
        for (unsigned int i=1; i<num_proc; i++)
                lbd.sdispls[i] = lbd.sdispls[i-1] + lbd.sendcounts[i-1];

        lbd.rdispls[0] = 0;
        for (unsigned int i=1; i<num_proc; i++)
                lbd.rdispls[i] = lbd.rdispls[i-1] + lbd.recvcounts[i-1];

        // Create receive buffer
        const int total_recv = std::accumulate(lbd.recvcounts.begin(), lbd.recvcounts.end(), 0);
        lbd.recv.resize(total_recv);

        // Communicate points (alltoallv)
        comm.alltoallv(lbd.send.data(), lbd.sendcounts.data(), lbd.sdispls.data(),
                       lbd.recv.data(), lbd.recvcounts.data(), lbd.rdispls.data());

        // Clear out current points and reload from the receive buffer
        for (unsigned int i=0; i<lbd.recv.size(); i+=dim)
                pl.emplace_back(point_type(lbd.recv[i], lbd.recv[i+1]));

        std::cout << "points: ";
        for (auto& p : pl) std::cout << p << " ";
        std::cout << std::endl;
        std::cout << "--------------" << std::endl;

        par::communicator halfcomm(comm, my_rank < rank_midpt);
        rcb1d_recurse(halfcomm, dimen, pl, lbd);
}

static void lb_rcb1d(unsigned int dimen, list_type& pl)
{
        LBData lbd(par::comm_world(), pl.size());
        rcb1d_recurse(par::comm_world(), dimen, pl, lbd);
}

void balance_set(PSTopology top, LBMethod method, unsigned int dimen, list_type& pl)
{
        const par::communicator& comm = par::comm_world();
        if (comm.size() == 1) return;

        const double start_time = par::wtime();

        if (top == PSTopology::Unary) {
                if (method == LBMethod::S2A2) lb_s2a21d(dimen, pl);
                else if (method == LBMethod::RCB) lb_rcb1d(dimen, pl);
                else throw std::runtime_error("Unknown load balancing method");
        }
        else if (top == PSTopology::Zoltan) {
                lb_zoltan(top, method, dimen, pl);
        }

        comm.barrier();
        const double end_time = par::wtime();

        if (comm.rank() == 0)
                std::cout << "LB time = "
                          << std::setprecision(6)
                          << (end_time - start_time) << std::endl;
}

} // namepsace cgl
