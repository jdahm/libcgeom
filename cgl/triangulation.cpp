#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
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


void Triangulation::write_txt(const std::string& filePrefix) {
        const par::communicator& comm_world = par::comm_world();

        std::stringstream ss;
        ss << comm_world.rank();
        std::string rank;
        ss >> rank;
        std::ofstream fst(filePrefix + "_" + rank + ".txt", std::ios::out);

        fst << num_points() << " " << num_edges() << std::endl;

        for (Point2d& p : point)
                fst << std::scientific << std::setprecision(16) << p[0] << " "
                    << std::scientific << std::setprecision(16) << p[1]
                    << std::endl;

        edge_iterator it = edge.begin_loop();
        while (edge.valid(it)) {
                fst << it->e->Org() << " " << it->e->Dest() << std::endl;
                it = edge.iterate(it);
        }

        fst.close();
}

Delaunay::Delaunay(PointSet& ps) : bounding_box(ps.front()) {
        // Sort the points again in the x-direction
        ps.sort(0);

        // Create merge stack (recursiveness happens here)
        create_merge_stack(ps, par::comm_world());

        // Loop to add points to the bounding box
        std::cout << "point: " << std::endl;
        std::cout << std::endl;
        for (auto& p : ps) std::cout << p << " ";
        for (const point_type& p : ps) bounding_box.add_point(p);
        std::cout << bounding_box << std::endl;

        // std::stack<MergeInfo> my_stack = merge_stack;
        // std::cout << par::comm_world().rank() << std::endl;
        // while(!my_stack.empty()) //body
        // {
        //         auto &p = my_stack.top();
        //         if (p.active)
        //                 std::cout << p.comm.rank() << " " << p.neighbor << std::endl;
        //         my_stack.pop();
        // }

        // Delaunay the processor portion
        edge_type el, er;
        init_dc(ps, el, er, 0);
        std::cout << "el: " << get_point(el->Org()) << "," << get_point(el->Dest()) << std::endl;
        std::cout << "er: " << get_point(er->Org()) << "," << get_point(er->Dest()) << std::endl;

        // // Compute global ids (used in the processor merges below)
        // fill_global_ids();

        while (!merge_stack.empty()) {
                // Get the next element
                MergeInfo& mi = merge_stack.top();
                // If active, perform the merge
                if (mi.active) proc_merge(mi.comm, mi.neighbor, mi.neighbor_dir);
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
                        while (in_circle(
                                       get_point(basel->Dest()),
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
                        while (in_circle(
                                       get_point(basel->Dest()),
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
        else if (ps.size() == 1) {
                throw std::runtime_error("Cannot handle 1 point. This happening means there's a bug somewhere else.");
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

void Delaunay::create_merge_stack(const PointSet& ps, const par::communicator& comm)
{
        unsigned int my_rank = comm.rank();
        unsigned int num_proc = comm.size();

        // Stop the recursion if num_proc == 1
        if (num_proc == 1) return;

        unsigned int half = num_proc / 2;
        // Left = half - 1, Right = half
        // std::cout << my_rank << "/" << num_proc-1 << " " << half << " " << static_cast<int>(my_rank<half) << std::endl;
        par::communicator new_comm(comm, my_rank < half);

        merge_stack.emplace(MergeInfo({comm, 0, Direction::Left, (my_rank == half - 1) || (my_rank == half)}));
        MergeInfo& m = merge_stack.top();

        if (my_rank == half - 1) {
                m.neighbor = half;
                m.neighbor_dir = Direction::Right;
                m.active = true;
        }
        else if (my_rank == half) {
                m.neighbor = half - 1;
                m.neighbor_dir = Direction::Left;
                m.active = true;
        }

        create_merge_stack(ps, new_comm);
}

void Delaunay::proc_merge(const par::communicator& comm, unsigned int neighbor, Direction dir)
{
        // 1. Exchange facing part of bounding box
        std::array<real, 2> bounds = bounding_box.project(static_cast<int>(dir) % 2);

        std::array<real, 2> neighbor_bounds;

        comm.isend(bounds.data(), neighbor, 2);
        { par::request rr = comm.irecv(neighbor_bounds.data(), neighbor, 2); rr.wait(); }

        // Have bounds
        std::cout << "neighbor = " << neighbor << ": " << neighbor_bounds[0] << " " << neighbor_bounds[1] << std::endl;
        
        // // Determine a bounding point outside the AABB2d
        // point_type p = bounding_box.bounding_point(dir);

        // // Create a subdivision of the part of the convex hull facing the neighbor
        // Subdivision s = facing_hull(p);

        // // // Send the subdivision to neighbor
        // // send_subdivision(s, comm, neighbor);

        // // // Receive subdivision from neighbor
        // // Subdivision s = recv_subdivision(comm, neighbor);

        // // // Merge the subdivisions
        // // Delaunay& D = *this;
        // // Delaunay& O = dynamic_cast<Delaunay>(s);

        // // D += O;
}


} // namespace cgl
