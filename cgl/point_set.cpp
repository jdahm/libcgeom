#include <zoltan.h>
#include <cmath>
#include <vector>
#include <array>
#include <numeric>
#include <algorithm>
#include <iomanip>

#include "cgl/point_set.hpp"
#include "par/environment.hpp"

namespace cgl
{

// Zoltan query functions
namespace detail
{

static int ps_num_obj(void *data, int *ierr) {
        PointSet *ps = static_cast<PointSet *>(data);
        (*ierr) = ZOLTAN_OK;
        return ps->size();
}

static void ps_obj_list(void *data, int num_gid_entries,
                        int num_lid_entries,
                        ZOLTAN_ID_PTR global_ids,
                        ZOLTAN_ID_PTR local_ids,
                        int wgt_dim,
                        float *obj_wgts,
                        int *ierr) {
        PointSet *ps = static_cast<PointSet *>(data);
        PointSet::size_type global_offset = ps->global_offset();
        for (PointSet::size_type i=0; i<ps->size(); i++) {
                global_ids[i] = global_offset + i;
                // skip local_ids
        }
        // skip setting obj_wgts
        (*ierr) = ZOLTAN_OK;
        static_cast<void>(num_gid_entries);
        static_cast<void>(num_lid_entries);
        static_cast<void>(local_ids);
        static_cast<void>(wgt_dim);
        static_cast<void>(obj_wgts);
}

static int ps_num_geom(void *data, int *ierr) {
        PointSet *ps = static_cast<PointSet *>(data);
        (*ierr) = 0;
        return ps->dim();
}

static void ps_geom_multi(void *data, int num_gid_entries,
                          int num_lid_entries,
                          int num_obj,
                          ZOLTAN_ID_PTR global_ids,
                          ZOLTAN_ID_PTR local_ids,
                          int num_dim,
                          double *geom_vec,
                          int *ierr) {
        PointSet *ps = static_cast<PointSet *>(data);
        PointSet::size_type global_offset = ps->global_offset();
        for (int i=0; i<num_obj; i++) {
                // PointSet::const_iterator it = ps->begin() + i;
                PointSet::const_iterator it = std::next(ps->begin(), i);
                for (PointSet::size_type d=0; d<PointSet::dim(); d++)
                        geom_vec[i*ps->dim() + d] = (*it)[d];
        }
        (*ierr) = 0;
        static_cast<void>(num_gid_entries);
        static_cast<void>(num_lid_entries);
        static_cast<void>(global_ids);
        static_cast<void>(local_ids);
        static_cast<void>(num_dim);
}

// Migration query functions
static void ps_mid_migrate_pp(void *data,
                              int num_gid_entries,
                              int num_lid_entries,
                              int num_import,
                              ZOLTAN_ID_PTR import_global_ids,
                              ZOLTAN_ID_PTR import_local_ids,
                              int *import_procs,
                              int *import_to_part,
                              int num_export,
                              ZOLTAN_ID_PTR export_global_ids,
                              ZOLTAN_ID_PTR export_local_ids,
                              int *export_procs,
                              int *export_to_part,
                              int *ierr) {
        PointSet *ps = static_cast<PointSet *>(data);
        PointSet::size_type global_offset = ps->global_offset();
        PointSet::size_type i = 0, j = 0, nexport = num_export;
        PointSet::iterator it = ps->begin();
        while (it != ps->end()) {
                if ((j < nexport) &&
                    // (it - ps->begin() + j ==
                    (i == export_global_ids[j] - global_offset)) {
                        it = ps->erase(it);
                        j++;
                }
                else {
                        ++it;
                }
                i++;
        }
        (*ierr) = ZOLTAN_OK;
        static_cast<void>(num_gid_entries);
        static_cast<void>(num_lid_entries);
        static_cast<void>(num_import);
        static_cast<void>(import_global_ids);
        static_cast<void>(import_local_ids);
        static_cast<void>(import_procs);
        static_cast<void>(import_to_part);
        static_cast<void>(export_local_ids);
        static_cast<void>(export_procs);
        static_cast<void>(export_to_part);
}

static void ps_obj_size_multi(void *data, int num_gid_entries,
                              int num_lid_entries,
                              int num_ids,
                              ZOLTAN_ID_PTR global_id,
                              ZOLTAN_ID_PTR local_id,
                              int *sizes,
                              int *ierr) {
        PointSet *ps = static_cast<PointSet *>(data);
        for (int i=0; i<num_ids; i++) sizes[i] = ps->dim() * sizeof(real);
        (*ierr) = ZOLTAN_OK;
        static_cast<void>(num_gid_entries);
        static_cast<void>(num_lid_entries);
        static_cast<void>(global_id);
        static_cast<void>(local_id);
}

static void ps_pack_obj_multi(void *data, int num_gid_entries,
                                    int num_lid_entries,
                                    int num_ids,
                                    ZOLTAN_ID_PTR global_ids,
                                    ZOLTAN_ID_PTR local_ids,
                                    int *dest,
                                    int *sizes,
                                    int *idx,
                                    char *buf,
                                    int *ierr) {
        PointSet *ps = static_cast<PointSet *>(data);
        PointSet::size_type global_offset = ps->global_offset();
        for (int i=0; i<num_ids; i++) {
                PointSet::const_iterator it =
                        std::next(ps->begin(), global_ids[i] - global_offset);
                // PointSet::const_iterator it = ps->begin() + global_ids[i] - global_offset;
                for (PointSet::size_type d=0; d<PointSet::dim(); d++)
                        memcpy(buf + idx[i] + d*sizeof(real), &(*it)[d], sizeof(real));
        }
        (*ierr) = ZOLTAN_OK;
        static_cast<void>(num_gid_entries);
        static_cast<void>(num_lid_entries);
        static_cast<void>(local_ids);
        static_cast<void>(dest);
        static_cast<void>(sizes);
}

static void ps_unpack_obj_multi(void *data, int num_gid_entries,
                                      int num_ids,
                                      ZOLTAN_ID_PTR global_ids,
                                      int *sizes,
                                      int *idx,
                                      char *buf,
                                      int *ierr) {
        PointSet *ps = static_cast<PointSet *>(data);
        for (int i=0; i<num_ids; i++) {
                real x, y;
                memcpy(&x, buf + idx[i] + 0, sizeof(real));
                memcpy(&y, buf + idx[i] + sizeof(real), sizeof(real));
                ps->add(Point2d(x, y));
        }
        (*ierr) = ZOLTAN_OK;
        static_cast<void>(num_gid_entries);
        static_cast<void>(global_ids);
        static_cast<void>(sizes);
}

} // namepsace detail

void PointSet::z_init()
{
        const par::communicator &comm_world = par::comm_world();

        // Allocate the Zoltan stack
        zz = Zoltan_Create(comm_world.raw());

        // Set some default sane parameters
        Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
        // Zoltan_Set_Param(zz, "KEEP_CUTS", "1");
        // Zoltan_Set_Param(zz, "LB_APPROACH", "REPARTITION");
        // Zoltan_Set_Param(zz, "MIGRATE_ONLY_PROC_CHANGES", "1");
        Zoltan_Set_Param(zz, "AUTO_MIGRATE", "TRUE");
        // Set higher for more debugging output
        Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");

        // Set partition query methods
        Zoltan_Set_Num_Obj_Fn(zz, detail::ps_num_obj,
                              static_cast<void *>(this));
        Zoltan_Set_Obj_List_Fn(zz, detail::ps_obj_list,
                               static_cast<void *>(this));
        Zoltan_Set_Num_Geom_Fn(zz, detail::ps_num_geom,
                               static_cast<void *>(this));
        Zoltan_Set_Geom_Multi_Fn(zz, detail::ps_geom_multi,
                                 static_cast<void *>(this));

        // Migration query methods
        Zoltan_Set_Mid_Migrate_PP_Fn(zz, detail::ps_mid_migrate_pp,
                                     static_cast<void *>(this));
        Zoltan_Set_Obj_Size_Multi_Fn(zz, detail::ps_obj_size_multi,
                                     static_cast<void *>(this));
        Zoltan_Set_Pack_Obj_Multi_Fn(zz, detail::ps_pack_obj_multi,
                                     static_cast<void *>(this));
        Zoltan_Set_Unpack_Obj_Multi_Fn(zz, detail::ps_unpack_obj_multi,
                                       static_cast<void *>(this));
}

void PointSet::z_balance(ProcTopology top) {
        int
                ierr,
                changes,
                num_gid_entries, num_lid_entries,
                num_import,
                *import_procs,
                *import_to_part,
                num_export,
                *export_procs,
                *export_to_part;
        unsigned int
                *import_global_ids,
                *import_local_ids,
                *export_global_ids,
                *export_local_ids;

        // Set method
        if (top == ProcTopology::RCB)
                Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
        else if (top == ProcTopology::RIB)
                Zoltan_Set_Param(zz, "LB_METHOD", "RIB");
        else
                throw std::runtime_error("No such Zoltan method");

        ierr = Zoltan_LB_Partition(zz,
                                   &changes,
                                   &num_gid_entries,
                                   &num_lid_entries,
                                   &num_import,
                                   &import_global_ids,
                                   &import_local_ids,
                                   &import_procs,
                                   &import_to_part,
                                   &num_export,
                                   &export_global_ids,
                                   &export_local_ids,
                                   &export_procs,
                                   &export_to_part);

        const par::communicator &comm_world = par::comm_world();
        if (ierr != ZOLTAN_OK) comm_world.abort("Zoltan error", 1);
}

void PointSet::z_destroy() { Zoltan_Destroy(&zz); }

PointSet::PointSet() : point(), global_valid(false),
                       zz(0), sorted_dimension(-1) { }

PointSet::PointSet(const container_type& p) :
        point(p), global_valid(false),
        zz(0), sorted_dimension(-1) { }

PointSet::PointSet(container_type&& p) :
        point(p), global_valid(false),
        zz(0), sorted_dimension(-1) { }

PointSet::~PointSet() { if (zz != 0) z_destroy(); }

PointSet::PointSet(const PointSet& other) :
        point(other.point), global_valid(other.global_valid),
        global_offset_(other.global_offset_) { }

PointSet& PointSet::operator=(const PointSet& other) {
        point = other.point;
        global_valid = other.global_valid;
        global_offset_ = other.global_offset_;
        // Do not copy Zoltan data
        return *this;
}

PointSet PointSet::split(PointSet::iterator it) {
        // Create a new container
        container_type pointR;

        // Take the right side of the list
        pointR.splice(pointR.begin(), point, it, point.end());

        // Return this->begin() in it
        it = begin();

        // Return a PointSet
        return PointSet(pointR);
}

PointSet PointSet::halve() {
        return split(std::next(begin(), (size() + 1) / 2));
}

void PointSet::sort(unsigned int dimension) {
        // Sort
        if (dimension == 0)
                point.sort([](const point_type& a, const point_type& b)
                           { if (std::abs(a[0] - b[0]) < real_eps) return a[1] < b[1];
                                   else return a[0] < b[0]; });
        else
                point.sort([](const point_type& a, const point_type& b)
                           { if (std::abs(a[1] - b[1]) < real_eps) return a[0] < b[0];
                                   else return a[1] < b[1]; });

        // Remove duplicates (has to be already sorted)
        point.unique([](const point_type& a, const point_type& b)
                     { return (std::abs(a[0] - b[0]) < real_eps) &&
                                     (std::abs(a[1] - b[1]) < real_eps); });

        sorted_dimension = dimension;
}

void PointSet::distribute(ProcTopology top)
{
        if (top == ProcTopology::Line) {
                partition_1d(0);
        }
        else if (top == ProcTopology::NestedGrid) {
                throw std::runtime_error("No such proc topology");
        }
        else {
                // Zoltan methods
                if (zz == 0) z_init();
                z_balance(top);
        }
        // Invalidate global_offset
        global_valid = false;
}

PointSet::size_type PointSet::size() const { return point.size(); }

PointSet::iterator PointSet::begin() { return point.begin(); }
PointSet::iterator PointSet::end() { return point.end(); }
PointSet::const_iterator PointSet::begin() const { return point.begin(); }
PointSet::const_iterator PointSet::end() const { return point.end(); }

const PointSet::point_type& PointSet::front() const { return point.front(); }
const PointSet::point_type& PointSet::back() const { return point.back(); }


void PointSet::add(const point_type& p) { point.push_back(p); sorted_dimension = -1; }

PointSet::iterator PointSet::erase(iterator it) { return point.erase(it); }

PointSet::size_type PointSet::global_offset()
{
        // Recompute global_offset if invalid and return it
        if (!global_valid) recompute_global_offset();
        return global_offset_;
}

void PointSet::recompute_global_offset()
{
        // Perform a scan(+)
        const par::communicator &comm_world = par::comm_world();
        const size_type s = size();
        comm_world.scan(&s, &global_offset_, par::sum());
        // Scan was inclusive, so must subtract off this value
        global_offset_ -= size();
        // Set valid flag
        global_valid = true;
}

template <typename T1, typename T2, typename T3>
std::vector< std::list<T1> > split_multi(std::list<T1>& original_list,
                                         const std::vector<T2>& cut,
                                         T3 pred)
{
        // Concept from @Yakk answer here:
        // http://stackoverflow.com/questions/23794679/splitting-an-stl-list-based-on-a-condition
        std::vector< std::list<T1> > result;

        std::list<T1> current;
        typedef typename std::list<T1>::const_iterator const_iterator;
        typename std::vector<T2>::size_type i=0;
        for (const_iterator it=original_list.begin(); it!=original_list.end();/* nothing */) {
                ++it; // about to become invalid/in wrong list
                if (i < cut.size())
                        if (pred(*it, cut[i])) {
                                result.emplace_back(std::move(current));
                                i++;
                        }
                current.splice(current.end(), original_list, original_list.begin());
        }
        result.emplace_back(std::move(current));

        for (unsigned int i=result.size()-1; i<cut.size(); i++)
                result.emplace_back(std::list<T1>());

        return result;
}

void PointSet::partition_1d(unsigned int dimen)
{
        typedef std::vector<unsigned int> id_container;
        const par::communicator &comm = par::comm_world();
        const unsigned int num_proc = comm.size();
        const unsigned int my_rank = comm.rank();

        // Return early if this doesn't make sense
        if (num_proc == 1) return;

        // Sort if not already sorted
        if (sorted_dimension != static_cast<int>(dimen)) sort(dimen);

        // for (auto& p: point) std::cout << p << std::endl;

        std::vector<int> glob_np(num_proc);
        std::vector<int> sendcounts(num_proc), recvcounts(num_proc);
        std::vector<int> sdispls(num_proc), rdispls(num_proc);

        real imbalance;
        unsigned int counter = 0;

        glob_np[my_rank] = size();
        const unsigned int my_num_point = glob_np[my_rank];
        comm.allgather(glob_np.data() + my_rank, 1, glob_np.data(), 1);

        const size_type glob_num_point =
                std::accumulate(glob_np.begin(), glob_np.end(), 0);        

        const auto minmax_np = std::minmax_element(glob_np.begin(), glob_np.end());
        unsigned int my_num_cut = std::floor(static_cast<real>(my_num_point) / *minmax_np.second * std::min(num_proc, my_num_point));

        do {
                const size_type target_npb = (my_num_point > 0) ? std::floor(static_cast<real>(my_num_point) / (my_num_cut + 1)) : 0;
                // std::cout << my_num_cut << " " << target_npb << std::endl;
                // points are sorted, so next determine cuts
                std::vector<real> cut(my_num_cut, 0);
                {
                        const_iterator pit = begin();
                        size_type np = 0;
                        for (unsigned int i=0; i<my_num_cut; i++) {
                                while (np < target_npb) {
                                        ++pit;
                                        np++;
                                }
                                const point_type& p1 = *pit;
                                cut[i] = p1[dimen];
                                np = 0;
                        }
                }

                // std::cout << "cuts:";
                // for (auto& c : cut) std::cout << c << " ";
                // std::cout << std::endl;

                // allgather number of cuts on each proc
                recvcounts[my_rank] = cut.size();
                comm.allgather(recvcounts.data() + my_rank, 1, recvcounts.data(), 1);

                const unsigned int glob_num_cut =
                        std::accumulate(recvcounts.begin(), recvcounts.end(), 0);

                if (glob_num_cut + 1 < num_proc)
                        comm.abort("Global number of bins < number of procs", 1);

                // std::cout << "global number of cuts = " << glob_num_cut << std::endl;
                std::vector<real> glob_cut(glob_num_cut);

                rdispls[0] = 0;
                for (unsigned int i=1; i<rdispls.size(); i++)
                        rdispls[i] = rdispls[i-1] + recvcounts[i-1];

                // allgatherv into glob_cut
                comm.allgatherv(cut.data(), cut.size(), glob_cut.data(),
                                recvcounts.data(), rdispls.data());

                // std::cout << "global cuts: ";
                // for (auto& c : glob_cut) std::cout << c << " ";
                // std::cout << std::endl;

                // Count the number of points in each global bin
                std::vector<unsigned int> my_pointbin(glob_num_cut + 1);
                {
                        unsigned int i = 0;
                        for (const point_type& p : point) {
                                if (i < glob_num_cut)
                                        while (p[dimen] > glob_cut[i] && i < glob_num_cut) i++;
                                my_pointbin[i]++;
                        }
                }
                // std::cout << "my_pointbin: ";
                // for (auto& c : my_pointbin) std::cout << c << " ";
                // std::cout << std::endl;
                std::vector<unsigned int> pointbin(my_pointbin);
                comm.allreduce(my_pointbin.data(), pointbin.data(), par::sum(), pointbin.size());

                // std::cout << "pointbin: ";
                // for (auto& c : pointbin) std::cout << c << " ";
                // std::cout << std::endl;

                const real target_avg = static_cast<real>(glob_num_point) / num_proc;
                // std::cout << "target_avg = " << target_avg << std::endl;
                for (id_container::size_type my_np=0, np=0, j=0, i=0; i<num_proc-1; i++) {
                        // Loop until adding the next bin would bring us further from
                        // the target number of points
                        while (std::abs(static_cast<real>(np) - target_avg) >=
                               std::abs(static_cast<real>(np + pointbin[j]) - target_avg)) {
                                np += pointbin[j];
                                my_np += my_pointbin[j];
                                j++;
                        }
                        // Record sendcount
                        sendcounts[i] = my_np * dim();
                        glob_np[i] = np;
                        // Reset counter
                        np = 0;
                        my_np = 0;
                }

                {
                        const int sc_so_far = std::accumulate(sendcounts.begin(), --sendcounts.end(), 0);
                        sendcounts.back() = my_num_point * dim() - sc_so_far;

                        const int gnp_so_far = std::accumulate(glob_np.begin(), --glob_np.end(), 0);
                        glob_np.back() = glob_num_point - gnp_so_far;
                }

                // std::cout << "glob_np: ";
                // for (auto& c : glob_np) std::cout << c << " ";
                // std::cout << std::endl;

                const auto minmax_after = std::minmax_element(glob_np.begin(), glob_np.end());

                imbalance = (*minmax_after.second - *minmax_after.first) / static_cast<real>(glob_num_point);
                counter++;
                my_num_cut = std::min(2 * my_num_cut + 1, my_num_point - 1);
                // std::cout << "imbalance = " << imbalance << std::endl;

        } while (imbalance > 0.05 && counter < 3);

        // Perform alltoall for recvcounts from sendcounts
        comm.alltoall(sendcounts.data(), 1, recvcounts.data(), 1);

        // std::cout << "glob_np: ";
        // for (auto& c : glob_np) std::cout << c << " ";
        // std::cout << std::endl;

        // std::cout << "sendcounts: ";
        // for (auto& c : sendcounts) std::cout << c << " ";
        // std::cout << std::endl;

        // std::cout << "recvcounts: ";
        // for (auto& c : recvcounts) std::cout << c << " ";
        // std::cout << std::endl;

        // Offsets in memory (for alltoallv)
        sdispls[0] = 0;
        for (id_container::size_type i=1; i<sdispls.size(); i++)
                sdispls[i] = sdispls[i-1] + sendcounts[i-1];

        rdispls[0] = 0;
        for (id_container::size_type i=1; i<rdispls.size(); i++)
                rdispls[i] = rdispls[i-1] + recvcounts[i-1];

        std::vector<real> send_points(dim() * size(), 0);
        {
                unsigned int i = 0;
                for (const point_type& p : point)
                        for (unsigned int d=0; d<dim(); d++, i++)
                                send_points[i] = p[d];
        }

        // Create a receive buffer
        const int total_recv = std::accumulate(recvcounts.begin(), recvcounts.end(), 0);
        std::vector<real> recv_points(total_recv, 0);

        // Communicate points (alltoallv)
        comm.alltoallv(send_points.data(), sendcounts.data(), sdispls.data(),
                       recv_points.data(), recvcounts.data(), rdispls.data());

        // Clear out current points and reload from the receive buffer
        point.clear();
        for (std::vector<real>::size_type i=0; i<recv_points.size(); i+=dim())
                point.emplace_back(point_type(recv_points[i], recv_points[i+1]));

        // for (auto& p: point) std::cout << p << std::endl;

        // // // Sort after receiving points
        // // sort(dimension);
}

void write_csv(const PointSet &ps, const std::string &prefix)
{
        const par::communicator &comm_world = par::comm_world();

        std::stringstream ss;
        ss << comm_world.rank();

        std::string rank;
        ss >> rank;
        std::ofstream fst(prefix + "_" + rank + ".csv", std::ios::out);

        for (PointSet::const_iterator it=ps.begin(); it!=ps.end(); ++it) {
                const Point2d &a = *it;
                fst << a[0] << "," << a[1] << std::endl;
        }

        fst.close();
}

void read_csv(PointSet& ps, const std::string& prefix)
{
        const par::communicator &comm_world = par::comm_world();

        std::stringstream ss;
        ss << comm_world.rank();

        std::string rank;
        ss >> rank;
        std::ifstream fst(prefix + "_" + rank + ".csv", std::ios::in);

        while (fst.good()) {
                real x, y;
                fst >> x;
                fst.ignore();
                fst >> y;
                if (fst.eof()) break;
                ps.add(PointSet::point_type(x, y));
        }

        fst.close();
}

} // namespace cgl
