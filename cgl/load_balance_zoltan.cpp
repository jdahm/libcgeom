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

namespace cgl {

typedef Point2d point_type;
typedef std::list<point_type> list_type;

class ZData {
public:
        list_type list;
        const par::communicator& comm;

        ZData(const par::communicator& comm, list_type&& l) :
                list(l), comm(comm), glob_offset_(0), offset_valid(false) { }

        list_type::size_type global_offset() {
                if (!offset_valid) recompute_offset();
                return glob_offset_;
        }

private:
        void recompute_offset() {
                // Perform a scan(+)
                const list_type::size_type s = list.size();
                comm.scan(&s, &glob_offset_, par::sum());
                // Scan was inclusive, so must subtract off this value
                glob_offset_ -= s;
                // Set valid flag
                offset_valid = true;                
        }

        list_type::size_type glob_offset_;
        bool offset_valid;
};

// Zoltan query functions
static int pl_num_obj(void *data, int *ierr) {
        ZData* zd = static_cast<ZData *>(data);
        (*ierr) = ZOLTAN_OK;
        static_cast<void>(data);
        return zd->list.size();
}

static void pl_obj_list(void *data, int num_gid_entries,
                        int num_lid_entries,
                        ZOLTAN_ID_PTR global_ids,
                        ZOLTAN_ID_PTR local_ids,
                        int wgt_dim,
                        float *obj_wgts,
                        int *ierr) {
        ZData* zd = static_cast<ZData *>(data);
        list_type::size_type offset = zd->global_offset();
        const list_type::size_type s = zd->list.size();
        for (list_type::size_type i=0; i<s; i++) {
                global_ids[i] = offset + i;
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

static int pl_num_geom(void *data, int *ierr) {
        (*ierr) = ZOLTAN_OK;
        static_cast<void>(data);
        return point_type::dim;
}

static void pl_geom_multi(void *data, int num_gid_entries,
                          int num_lid_entries,
                          int num_obj,
                          ZOLTAN_ID_PTR global_ids,
                          ZOLTAN_ID_PTR local_ids,
                          int num_dim,
                          double *geom_vec,
                          int *ierr) {
        ZData *zd = static_cast<ZData *>(data);
        // ZData::size_type global_offset = pl->global_offset();
        const unsigned int dim = point_type::dim;
        for (int i=0; i<num_obj; i++) {
                // ZData::const_iterator it = pl->begin() + i;
                list_type::const_iterator it = std::next(zd->list.begin(), i);
                for (list_type::size_type d=0; d<dim; d++)
                        geom_vec[i*dim + d] = (*it)[d];
        }
        (*ierr) = 0;
        static_cast<void>(num_gid_entries);
        static_cast<void>(num_lid_entries);
        static_cast<void>(global_ids);
        static_cast<void>(local_ids);
        static_cast<void>(num_dim);
}

// Migration query functions
static void pl_mid_migrate_pp(void *data,
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
        ZData *zd = static_cast<ZData *>(data);
        list_type::size_type offset = zd->global_offset();
        list_type::size_type i = 0, j = 0, nexport = num_export;
        list_type::iterator it = zd->list.begin();
        while (it != zd->list.end()) {
                if ((j < nexport) &&
                    // (it - pl->begin() + j ==
                    (i == export_global_ids[j] - offset)) {
                        it = zd->list.erase(it);
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

static void pl_obj_size_multi(void *data, int num_gid_entries,
                              int num_lid_entries,
                              int num_ids,
                              ZOLTAN_ID_PTR global_id,
                              ZOLTAN_ID_PTR local_id,
                              int *sizes,
                              int *ierr) {
        for (int i=0; i<num_ids; i++) sizes[i] = point_type::dim * sizeof(real);
        (*ierr) = ZOLTAN_OK;
        static_cast<void>(data);
        static_cast<void>(num_gid_entries);
        static_cast<void>(num_lid_entries);
        static_cast<void>(global_id);
        static_cast<void>(local_id);
}

static void pl_pack_obj_multi(void *data, int num_gid_entries,
                              int num_lid_entries,
                              int num_ids,
                              ZOLTAN_ID_PTR global_ids,
                              ZOLTAN_ID_PTR local_ids,
                              int *dest,
                              int *sizes,
                              int *idx,
                              char *buf,
                              int *ierr) {
        ZData *zd = static_cast<ZData *>(data);
        list_type::size_type offset = zd->global_offset();
        const unsigned int dim = point_type::dim;
        for (int i=0; i<num_ids; i++) {
                list_type::const_iterator it =
                        std::next(zd->list.begin(), global_ids[i] - offset);
                // ZData::const_iterator it = pl->begin() + global_ids[i] - global_offset;
                for (list_type::size_type d=0; d<dim; d++)
                        memcpy(buf + idx[i] + d*sizeof(real), &(*it)[d], sizeof(real));
        }
        (*ierr) = ZOLTAN_OK;
        static_cast<void>(num_gid_entries);
        static_cast<void>(num_lid_entries);
        static_cast<void>(local_ids);
        static_cast<void>(dest);
        static_cast<void>(sizes);
}

static void pl_unpack_obj_multi(void *data, int num_gid_entries,
                                int num_ids,
                                ZOLTAN_ID_PTR global_ids,
                                int *sizes,
                                int *idx,
                                char *buf,
                                int *ierr) {
        ZData *zd = static_cast<ZData *>(data);
        for (int i=0; i<num_ids; i++) {
                real x, y;
                memcpy(&x, buf + idx[i] + 0, sizeof(real));
                memcpy(&y, buf + idx[i] + sizeof(real), sizeof(real));
                zd->list.emplace_back(Point2d(x, y));
        }
        (*ierr) = ZOLTAN_OK;
        static_cast<void>(num_gid_entries);
        static_cast<void>(global_ids);
        static_cast<void>(sizes);
}

void lb_zoltan(PSTopology top, LBMethod method, unsigned int dimen, list_type& pl)
{
        const par::communicator& comm = par::comm_world();

        float ver;
        Zoltan_Initialize(0, 0, &ver);

        struct Zoltan_Struct *zz;

        // Create ZData (moves pl into zd)
        ZData zd(par::comm_world(), std::move(pl));

        // Allocate the Zoltan data
        zz = Zoltan_Create(zd.comm.raw());

        // Set some default sane parameters
        if (method == LBMethod::RCB)
                Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
        else
                throw std::runtime_error("Unknown load balancing method");

        // Zoltan_Set_Param(zz, "KEEP_CUTS", "1");
        // Zoltan_Set_Param(zz, "LB_APPROACH", "REPARTITION");
        // Zoltan_Set_Param(zz, "MIGRATE_ONLY_PROC_CHANGES", "1");
        Zoltan_Set_Param(zz, "AUTO_MIGRATE", "TRUE");
        // Set higher for more debugging output
        Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");

        // Set partition query methods
        Zoltan_Set_Num_Obj_Fn(zz, pl_num_obj,
                              static_cast<void *>(&zd));
        Zoltan_Set_Obj_List_Fn(zz, pl_obj_list,
                               static_cast<void *>(&zd));
        Zoltan_Set_Num_Geom_Fn(zz, pl_num_geom,
                               static_cast<void *>(&zd));
        Zoltan_Set_Geom_Multi_Fn(zz, pl_geom_multi,
                                 static_cast<void *>(&zd));

        // Migration query methods
        Zoltan_Set_Mid_Migrate_PP_Fn(zz, pl_mid_migrate_pp,
                                     static_cast<void *>(&zd));
        Zoltan_Set_Obj_Size_Multi_Fn(zz, pl_obj_size_multi,
                                     static_cast<void *>(&zd));
        Zoltan_Set_Pack_Obj_Multi_Fn(zz, pl_pack_obj_multi,
                                     static_cast<void *>(&zd));
        Zoltan_Set_Unpack_Obj_Multi_Fn(zz, pl_unpack_obj_multi,
                                       static_cast<void *>(&zd));

        int
                zerr,
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

        zerr = Zoltan_LB_Partition(zz,
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

        if (zerr != ZOLTAN_OK) comm.abort("Zoltan error", 1);

        // Move the data back again out of the struct
        pl = std::move(zd.list);

        Zoltan_Destroy(&zz);
}

} // namespace cgl
