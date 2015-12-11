#include <string>

#include "cgl/triangulation.hpp"
#include "par/environment.hpp"

namespace cgl
{

Hull::Hull(std::vector<real> raw) : base(point_type(raw[0], raw[1])), hlist()
{
        container_type::size_type i = 2;
        while (i < raw.size()) {
#ifndef NDEBUG
                if (i+4 > raw.size()) throw std::runtime_error("Raw hull size error");
#endif
                hlist.emplace_back(data_type(point_type(raw[i+0], raw[i+1]),
                                             point_type(raw[i+2], raw[i+3])));
                i += 4;
        }
}

void Hull::swap_new(iterator it, std::vector<real> raw)
{
#ifndef NDEBUG
        if (raw.size() != 4) throw std::runtime_error("Raw data has incorrect size");
#endif
        data_type& d = *it;
        iterator newit = hlist.insert(++it, d);
        d.second = point_type(raw[0], raw[1]);
        data_type& newd = *newit;
        newd.second = point_type(raw[2], raw[3]);
}

Hull::iterator Hull::begin() { return hlist.begin(); }
Hull::iterator Hull::end() { return hlist.end(); }
Hull::reverse_iterator Hull::rbegin() { return hlist.rbegin(); }
Hull::reverse_iterator Hull::rend() { return hlist.rend(); }

const Hull::point_type& Hull::org() const { return base; }


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

Delaunay::Delaunay(PointSet& ps) {
        if (ps.size() == 1) throw std::runtime_error("Need more than one point.");

        // Sort the points again in the x-direction
        ps.sort(0);

        std::cout << ps.size() << std::endl;

        // Create merge stack (recursiveness happens here)
        create_merge_stack();

        // Delaunay the processor portion
        edge_type el, er;
        init_dc(ps, el, er, 0);

        // // Compute global ids (used in the processor merges below)
        // fill_global_ids();

        while (!merge_stack.empty()) {
                // Get the next element
                for (const Neighbor& n : merge_stack.top()) {
                        if (n.dir == Direction::Left)
                                merge_proc(n.neighbor, n.dir, er, el);
                        else if (n.dir == Direction::Right)
                                merge_proc(n.neighbor, n.dir, el, er);
                        else
                                throw std::runtime_error("Cannot merge this direction");
                }
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

                std::cout << "Before traverse" << std::endl;
                std::cout << get_point(ldi->Org()) << "," << get_point(ldi->Dest()) << std::endl;
                std::cout << get_point(rdi->Org()) << "," << get_point(rdi->Dest()) << std::endl;
                // Traverse the convex hulls of left and right point sets in preparation
                // to compute the base tangent
                while (true) {
                        if (left_of(rdi->Org(), ldi)) ldi = ldi->Lnext();
                        else if (right_of(ldi->Org(), rdi)) rdi = rdi->Rprev();
                        else break;
                }
                std::cout << "After traverse" << std::endl;
                std::cout << get_point(ldi->Org()) << "," << get_point(ldi->Dest()) << std::endl;
                std::cout << get_point(rdi->Org()) << "," << get_point(rdi->Dest()) << std::endl;

                // Merge the two triangulations
                merge(ldo, ldi, rdi, rdo);

                // Return [ldo, rdo]
                el = ldo;
                er = rdo;
        }
}

void Delaunay::merge_hull(edge_type& ldo, const edge_type& ldi, Hull& h)
{
        Hull::iterator hit = h.begin();
        edge_type basel = extend_edge(ldi->Sym(), h.org())->Sym();
        std::cout << get_point(basel->Org()) << " " << get_point(basel->Dest()) << std::endl;
        std::cout << get_point(basel->Oprev()->Org()) << " " << get_point(basel->Oprev()->Dest()) << std::endl;
        std::cout << get_point(basel->Dnext()->Sym()->Org()) << " " << get_point(basel->Dnext()->Sym()->Dest()) << std::endl;

        // Keep this here for old times' sake
        if (ldi->Org() == ldo->Org()) ldo = basel->Sym();

        while (true) {
                edge_type lcand = basel->Sym()->Onext();
                if (right_of(lcand->Dest(), basel))
                        while (in_circle(get_point(basel->Dest()),
                                         get_point(basel->Org()),
                                         get_point(lcand->Dest()),
                                         get_point(lcand->Onext()->Dest()))) {
                                edge_type t = lcand->Onext();
                                remove_edge(lcand);
                                lcand = t;
                                throw std::runtime_error("Need to transfer here");
                        }

                // Symmetrically, locate the first R point to be hit, and delete
                // R edges
                const Hull::data_type& d = *hit;
                if (right_of(d.first, basel))
                        while (in_circle(get_point(basel->Dest()),
                                         get_point(basel->Org()),
                                         d.first,
                                         d.second)) {
                                throw std::runtime_error("Need to receive here");
                        }

                // If both lcand and rcand are invalid, then basel is the upper
                // common tangent
                if (!right_of(lcand->Dest(), basel) &&
                    !right_of(d.first, basel)) break;

                std::cout << static_cast<int>(!right_of(lcand->Dest(), basel)) << " || "  << static_cast<int>(right_of(d.first, basel)) << " && " << static_cast<int>(in_circle(get_point(lcand->Dest()), get_point(lcand->Org()), get_point(basel->Org()), d.first)) << std::endl;
                std::cout << get_point(lcand->Dest()) << "," << get_point(lcand->Org()) << "," << get_point(basel->Org()) << "," << d.first << std::endl;
                std::cout << "lcand here = " << get_point(lcand->Org()) << "," << get_point(lcand->Dest()) << std::endl;

                // The next cross edge is to be connected to either lcand.Dest
                // or rcand.Dest If both are valid, then choose the appropriate
                // one using the in_circle test.
                if (!right_of(lcand->Dest(), basel) ||
                    (right_of(d.first, basel) &&
                     in_circle(
                             get_point(lcand->Dest()),
                             get_point(lcand->Org()),
                             get_point(basel->Org()),
                             d.first))) {
                        // Add cross edge basel from rcand.Dest to basel.Dest
                        edge_type rcand = extend_edge(basel->Sym(), d.first);
                        std::cout << "rcand = " << get_point(rcand->Org()) << "," << get_point(rcand->Dest()) << std::endl;
                        std::cout << "basel = " << get_point(basel->Org()) << "," << get_point(basel->Dest()) << std::endl;
                        basel = connect_edges(rcand, basel->Sym());
                        std::cout << "new basel = " << get_point(basel->Org()) << "," << get_point(basel->Dest()) << std::endl;
                        ++hit;
                }
                else {
                        // Add cross edge basel from basel.Org to lcand.Dest
                        std::cout << "lcand = " << get_point(lcand->Org()) << "," << get_point(lcand->Dest()) << std::endl;
                        std::cout << "basel = " << get_point(basel->Org()) << "," << get_point(basel->Dest()) << std::endl;
                        basel = connect_edges(basel->Sym(), lcand->Sym());
                        std::cout << "new basel = " << get_point(basel->Org()) << "," << get_point(basel->Dest()) << std::endl;
                }

        } // Merge loop
}

void Delaunay::merge_hull(Hull& h, const edge_type& rdi, edge_type& rdo)
{
        Hull::iterator hit = h.begin();

        edge_type basel = extend_edge(rdi->Sym(), h.org());
        std::cout << get_point(basel->Org()) << " " << get_point(basel->Dest()) << std::endl;
        std::cout << get_point(basel->Oprev()->Org()) << " " << get_point(basel->Oprev()->Dest()) << std::endl;
        std::cout << get_point(basel->Dnext()->Sym()->Org()) << " " << get_point(basel->Dnext()->Sym()->Dest()) << std::endl;

        // Keep this here for old times' sake
        if (rdi->Org() == rdo->Org()) rdo = basel;

        edge_type lcand_former = basel;

        // This is the merge loop
        while (true) {
                // Locate the first L point (lcand.Dest) to be encountered by
                // the rising bubble, and delete L edges out of basel.Dest that
                // fail the circle test.
                const Hull::data_type& d = *hit;
                if (right_of(d.first, basel))
                        while (in_circle(get_point(basel->Dest()),
                                         get_point(basel->Org()),
                                         d.first,
                                         d.second)) {
                                throw std::runtime_error("Need to receive here");
                        }

                // Symmetrically, locate the first R point to be hit, and delete
                // R edges
                edge_type rcand = basel->Oprev();
                if (right_of(rcand->Dest(), basel)) {
                        while (in_circle(get_point(basel->Dest()),
                                         get_point(basel->Org()),
                                         get_point(rcand->Dest()),
                                         get_point(rcand->Oprev()->Dest()))) {
                                edge_type t = rcand->Oprev();
                                remove_edge(rcand);
                                rcand = t;
                                throw std::runtime_error("Need to transfer here");
                        }
                }

                // If both lcand and rcand are invalid, then basel is the upper
                // common tangent
                if (!right_of(d.first, basel) &&
                    !right_of(rcand->Dest(), basel)) break;
                std::cout << static_cast<int>(!right_of(d.first, basel)) << " || "  << static_cast<int>(right_of(rcand->Dest(), basel)) << " && " << static_cast<int>(in_circle(d.first, get_point(basel->Dest()), get_point(rcand->Org()), get_point(rcand->Dest()))) << std::endl;
                std::cout << d.first << "," << get_point(basel->Dest()) << "," << get_point(rcand->Org()) << "," << get_point(rcand->Dest()) << std::endl;

                // The next cross edge is to be connected to either lcand.Dest
                // or rcand.Dest If both are valid, then choose the appropriate
                // one using the in_circle test.
                if (!right_of(d.first, basel) ||
                    (right_of(rcand->Dest(), basel) &&
                     in_circle(
                             d.first,
                             get_point(basel->Dest()),
                             get_point(rcand->Org()),
                             get_point(rcand->Dest())))) {
                        // Add cross edge basel from rcand.Dest to basel.Dest
                        std::cout << "rcand = " << get_point(rcand->Org()) << "," << get_point(rcand->Dest()) << std::endl;
                        std::cout << "basel = " << get_point(basel->Org()) << "," << get_point(basel->Dest()) << std::endl;
                        basel = connect_edges(rcand, basel->Sym());
                        std::cout << "new basel = " << get_point(basel->Org()) << "," << get_point(basel->Dest()) << std::endl;
                }
                else {
                        // Add cross edge basel from basel.Org to lcand.Dest
                        edge_type lcand = extend_edge(lcand_former, d.first);
                        std::cout << "lcand = " << get_point(lcand->Org()) << "," << get_point(lcand->Dest()) << std::endl;
                        std::cout << "basel = " << get_point(basel->Org()) << "," << get_point(basel->Dest()) << std::endl;
                        basel = connect_edges(basel->Sym(), lcand->Sym());
                        std::cout << "new basel = " << get_point(basel->Org()) << "," << get_point(basel->Dest()) << std::endl;
                        ++hit;
                        lcand_former = lcand;
                }
        } // Merge loop
}


void Delaunay::create_merge_stack()
{
        // LIMITATION: This is very specific to a 1d merge of processors
        // in a "unary" layout
        const par::communicator &comm_world = par::comm_world();

        const int my_rank = comm_world.rank();
        const int num_proc = comm_world.size();

        if (num_proc == 1) return;

        merge_stack.push(std::vector<Neighbor>());

        if (my_rank % 2) {
                const int neighbor = my_rank - 1;
                if (neighbor >= 0)
                        merge_stack.top().emplace_back(Neighbor({neighbor, Direction::Left}));
        }
        else {
                const int neighbor = my_rank + 1;
                if (neighbor < num_proc)
                        merge_stack.top().emplace_back(Neighbor({neighbor, Direction::Right}));
        }

        if (num_proc == 2) return;

        merge_stack.push(std::vector<Neighbor>());
        if (my_rank % 2) {
                const int neighbor = my_rank + 1;
                if (neighbor < num_proc)
                        merge_stack.top().emplace_back(Neighbor({neighbor, Direction::Right}));
        }
        else {
                const int neighbor = my_rank - 1;
                if (neighbor >= 0)
                        merge_stack.top().emplace_back(Neighbor({neighbor, Direction::Left}));
        }
}

void Delaunay::merge_proc(unsigned int neighbor, Direction dir, edge_type& eo, edge_type& ei)
{
        // LIMITATION: This can only merge convex regions that are bounded "well" by an AABB2d

        const par::communicator& comm = par::comm_world();

        // // Determine a bounding point outside the AABB2d
        // point_type p = bounding_box.bounding_point(dir);

        // ei seem to be CW when ei takes the rdi and CCW when ei takes the
        // place of ldi. While this is fine, it's easier to extract the hull and
        // do the merge if the face is exactly opposite. Inside the merge the
        // edges are flipped anyway before any merge occurs. Flip it here instead.
        std::cout << "before flip = " << get_point(ei->Org()) << "," << get_point(ei->Dest()) << std::endl;
        ei = ei->Sym();

        // Create a subdivision of the part of the convex hull facing the neighbor
        std::vector<real> myhull;
        if (dir == Direction::Left)  myhull = extract_left_hull(ei);
        if (dir == Direction::Right) myhull = extract_right_hull(ei);

        const std::vector<real>::size_type myhullsize = myhull.size();
        comm.isend(&myhullsize, neighbor);

        std::vector<real>::size_type hullsize;
        {
                par::request rr = comm.irecv(&hullsize, neighbor);
                rr.wait();
        }
        std::cout << "hullsize = " << hullsize << std::endl;
        // Next, send the hull
        std::vector<real> hull(hullsize);
        comm.isend(myhull.data(), neighbor, myhullsize);
        {
                par::request rr = comm.irecv(hull.data(), neighbor, hullsize);
                rr.wait();
        }

        for (auto& p : hull) std::cout << p << " ";
        std::cout << std::endl;

                // while (true) {
                //         if (left_of(rdi->Org(), ldi)) ldi = ldi->Lnext();
                //         else if (right_of(ldi->Org(), rdi)) rdi = rdi->Rprev();
                //         else break;
                // }

        Hull h(hull);
        if (dir == Direction::Right) {
                // std::cout << "going in " << get_point(ei->Org()) << "," << get_point(ei->Dest()) << std::endl;
                // for (Hull::reverse_iterator it=h.rbegin(); it!=h.rend(); ++it) {
                //         while (left_of(it->first, ei)) ei = ei->Lnext();
                // }
                // while (left_of(h.org(), ei)) ei = ei->Lnext();

                std::cout << "going in " << get_point(ei->Org()) << "," << get_point(ei->Dest()) << std::endl;
                merge_hull(eo, ei, h);
        }
        else if (dir == Direction::Left) {
                // std::cout << "going in " << get_point(ei->Org()) << "," << get_point(ei->Dest()) << std::endl;
                // std::cout << "hull " << h.org() << std::endl;
                // if (right_of(h.org(), ei)) throw std::runtime_error("Need to implement");
                // else {
                //         std::cout << static_cast<int>(left_of(h.org(), ei)) << std::endl;
                //         std::cout << static_cast<int>(right_of(h.org(), ei)) << std::endl;
                // for (Hull::reverse_iterator it=h.rbegin(); it!=h.rend(); ++it) {
                //         while (right_of(it->first, ei)) ei = ei->Rprev();
                // }
                // while (right_of(h.org(), ei)) ei = ei->Rprev();
                // ei = ei->Sym();
                        // }

                std::cout << "going in " << get_point(ei->Org()) << "," << get_point(ei->Dest()) << std::endl;
                merge_hull(h, ei, eo);
        }
}


} // namespace cgl
