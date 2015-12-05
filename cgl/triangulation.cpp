#include <iomanip>
#include <string>
#include "cgl/triangulation.hpp"
#include "par/environment.hpp"

namespace cgl
{

void Triangulation::swap(const edge_type& e)
{
        edge_type a = e->Oprev();
        edge_type b = e->Sym()->Oprev();
        splice(e, a);
        splice(e->Sym(), b);
        splice(e, a->Lnext());
        splice(e->Sym(), b->Lnext());
        e->set_id(a->Dest(), b->Dest());
}

Delaunay::Delaunay(PointSet& ps) : bounding_box(ps.front()) {
        if (ps.size() == 1) throw std::runtime_error("Need more than one point.");

        // Sort the points again in the x-direction
        ps.sort(0);

        // Create merge stack (recursiveness happens here)
        create_merge_stack();

        // Loop to add points to the bounding box
        for (const point_type& p : ps) bounding_box.add_point(p);

        // Delaunay the processor portion
        edge_type el, er;
        init_dc(ps, el, er, 0);

        // // Compute global ids (used in the processor merges below)
        // fill_global_ids();

        while (!merge_stack.empty()) {
                // Get the next element
                MergeInfo& mi = merge_stack.top();
                // If active, perform the merge
                proc_merge(mi.neighbor,
                           mi.neighbor_dir,
                           mi.neighbor_dir == Direction::Left ? el : er);
                // Remove the merge info element
                merge_stack.pop();
        }
}


void Delaunay::merge(edge_type& ldo, const edge_type& ldi,
                     const edge_type& rdi, edge_type& rdo)
{
        // Create a first cross edge basel from rdi.Org to ldi.Org
        edge_type basel = connect_edges(rdi->Sym(), ldi);

        if (ldi->Org() == ldo->Org()) ldo = basel->Sym();
        if (rdi->Org() == rdo->Org()) rdo = basel;

        // This is the merge loop
        while (true) {
                // Locate the first L point (lcand.Dest) to be encountered by
                // the rising bubble, and delete L edges out of basel.Dest that
                // fail the circle test.
                edge_type lcand = basel->Sym()->Onext();
                if (right_of(lcand->Dest(), basel))
                        while (in_circle(get_point(basel->Dest()),
                                         get_point(basel->Org()),
                                         get_point(lcand->Dest()),
                                         get_point(lcand->Onext()->Dest()))) {
                                edge_type t = lcand->Onext();
                                remove_edge(lcand);
                                lcand = t;
                        }

                // Symmetrically, locate the first R point to be hit, and delete
                // R edges
                edge_type rcand = basel->Oprev();
                if (right_of(rcand->Dest(), basel))
                        while (in_circle(get_point(basel->Dest()),
                                         get_point(basel->Org()),
                                         get_point(rcand->Dest()),
                                         get_point(rcand->Oprev()->Dest()))) {
                                edge_type t = rcand->Oprev();
                                remove_edge(rcand);
                                rcand = t;
                        }

                // If both lcand and rcand are invalid, then basel is the upper
                // common tangent
                if (!right_of(lcand->Dest(), basel) &&
                    !right_of(rcand->Dest(), basel)) break;

                // The next cross edge is to be connected to either lcand.Dest
                // or rcand.Dest If both are valid, then choose the appropriate
                // one using the in_circle test.
                if (!right_of(lcand->Dest(), basel) ||
                    (right_of(rcand->Dest(), basel) &&
                     in_circle(
                             get_point(lcand->Dest()),
                             get_point(lcand->Org()),
                             get_point(rcand->Org()),
                             get_point(rcand->Dest())))) {
                        // Add cross edge basel from rcand.Dest to basel.Dest
                        basel = connect_edges(rcand, basel->Sym());
                }
                else {
                        // Add cross edge basel from basel.Org to lcand.Dest
                        basel = connect_edges(basel->Sym(), lcand->Sym());
                }
        } // Merge loop
}

void Delaunay::init_dc(PointSet& ps, edge_type& el, edge_type& er, int i)
{
        if (ps.size() == 2) {
                edge_type ea = add_edge(ps.front(), ps.back());
                el = ea;
                er = ea->Sym();
        }
        else if (ps.size() == 3) {
                PointSet::const_iterator it = ps.begin();
                const point_type& a = *it;
                const point_type& b = *++it;
                const point_type& c = *++it;
                edge_type ea = add_edge(a, b);
                edge_type eb = extend_edge(ea, c);
                // Now close the triangle
                if (ccw(a, b, c)) {
                        connect_edges(eb, ea);
                        el = ea;
                        er = eb->Sym();
                }
                else if (ccw(a, c, b)) {
                        edge_type ec = connect_edges(eb, ea);
                        el = ec->Sym();
                        er = ec;
                }
                else {
                        el = ea;
                        er = eb->Sym();
                }
        }
        else {
                // Halve the point set
                PointSet psr = ps.halve();
                // Existing ps has the left half
                // PointSet& psl = ps;

                // Compute Delaunay recusively until one of the conditions above is met
                edge_type ldo, ldi, rdi, rdo;
                init_dc(ps, ldo, ldi, i + 1);
                init_dc(psr, rdi, rdo, i + 1);

                // Traverse the convex hulls of left and right point sets in preparation
                // to compute the base tangent
                while (true) {
                        if (left_of(rdi->Org(), ldi)) ldi = ldi->Lnext();
                        else if (right_of(ldi->Org(), rdi)) rdi = rdi->Rprev();
                        else break;
                }

                // Merge the two triangulations
                merge(ldo, ldi, rdi, rdo);

                // Return [ldo, rdo]
                el = ldo;
                er = rdo;
        }
}

void Delaunay::create_merge_stack()
{
        const par::communicator &comm_world = par::comm_world();

        const int my_rank = comm_world.rank();
        const int num_proc = comm_world.size();

        if (num_proc == 1) return;

        if (my_rank % 2) {
                const int neighbor = my_rank - 1;
                if (neighbor >= 0)
                        merge_stack.emplace(MergeInfo({neighbor, Direction::Left}));
        }
        else {
                const int neighbor = my_rank + 1;
                if (neighbor < num_proc)
                        merge_stack.emplace(MergeInfo({neighbor, Direction::Right}));
        }

        if (num_proc == 2) return;

        if (my_rank % 2) {
                const int neighbor = my_rank + 1;
                if (neighbor < num_proc)
                        merge_stack.emplace(MergeInfo({neighbor, Direction::Right}));
        }
        else {
                const int neighbor = my_rank - 1;
                if (neighbor >= 0)
                        merge_stack.emplace(MergeInfo({neighbor, Direction::Left}));
        }
}

void Delaunay::proc_merge(unsigned int neighbor, Direction dir, edge_type& e)
{
        // // 1. Exchange facing part of bounding box
        // std::array<real, 2> bounds = bounding_box.project(static_cast<int>(dir) % 2);

        // std::array<real, 2> neighbor_bounds;

        // comm.isend(bounds.data(), neighbor, 2);
        // { par::request rr = comm.irecv(neighbor_bounds.data(), neighbor, 2); rr.wait(); }

        // // Have bounds
        // std::cout << "neighbor = " << neighbor << ": " << neighbor_bounds[0] << " " << neighbor_bounds[1] << std::endl;

        // Determine a bounding point outside the AABB2d
        point_type p = bounding_box.bounding_point(dir);

        // Create a subdivision of the part of the convex hull facing the neighbor
        Subdivision s = extract_hull(dir == Direction::Left ? e : e->Sym(), p);

        // // // send the subdivision to neighbor
        // // send_subdivision(s, comm, neighbor);

        // // // Receive subdivision from neighbor
        // // Subdivision s = recv_subdivision(comm, neighbor);

        // // // Merge the subdivisions
        // // Delaunay& D = *this;
        // // Delaunay& O = dynamic_cast<Delaunay>(s);

        // // D += O;
}


} // namespace cgl
