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

Hull::iterator Hull::swap(iterator it, std::vector<real> raw)
{
#ifndef NDEBUG
        if (raw.size() != 4) throw std::runtime_error("Raw data has incorrect size");
#endif
        // Copy the current data pointed to by it, to before it
        iterator newit = hlist.insert(it, *it);
        // The primary point is the candidate of *it
        newit->first = it->second;
        // The candidate is the first point transferred
        newit->second = point_type(raw[0], raw[1]);
        // The candidate for the next point (under *it) is the second point transferred
        it->second = point_type(raw[2], raw[3]);
        // Return an iterator to the new element
        return newit;
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
        const par::communicator& comm = par::comm_world();

        if (ps.size() == 1) throw std::runtime_error("Need more than one point.");

        const real sstart = par::wtime();

        // Sort the points again in the x-direction
        ps.sort(0);

        comm.barrier();
        const real send = par::wtime();

        if (comm.rank() == 0)
                std::cout << "Sort time = " << send - sstart << std::endl;

        comm.barrier();
        const real sdsstart = par::wtime();

        // Delaunay the processor portion
        edge_type el, er;
        init_dc(ps, el, er, 0);

        comm.barrier();
        const real sdsend = par::wtime();
        if (comm.rank() == 0)
                std::cout << "SelfDT time = " << sdsend - sdsstart << std::endl;

        comm.barrier();
        const real ndsstart = par::wtime();

        // Create merge stack
        create_merge_stack(comm);

        // return;
        while (!merge_stack.empty()) {
                // Get the next element
                const Neighbor& n = merge_stack.top();
                if (n.dir == Direction::Left)
                        merge_proc(comm, n.neighbor, n.dir, er, el);
                else if (n.dir == Direction::Right)
                        merge_proc(comm, n.neighbor, n.dir, el, er);
                else
                        throw std::runtime_error("Cannot merge this direction");
                // Remove the merge info element
                merge_stack.pop();
        }

        comm.barrier();
        const real ndsend = par::wtime();
        if (comm.rank() == 0)
                std::cout << "ParDT time = " << ndsend - ndsstart << std::endl;
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

void Delaunay::merge_hull(const par::communicator& comm, unsigned int neighbor,
                          edge_type& ldo, const edge_type& ldi, Hull& h)
{
        Hull::iterator hit = h.begin();
        edge_type basel = extend_edge(ldi->Sym(), h.org())->Sym();

        // Keep this here for old times' sake
        if (ldi->Org() == ldo->Org()) ldo = basel->Sym();

        std::vector<real> trans(4);
        while (true) {
#ifndef NDEBUG
                if (hit == h.end()) comm.abort("Prematurely hit end of hull", 1);
#endif
                edge_type lcand = basel->Sym()->Onext();
                if (right_of(lcand->Dest(), basel))
                        while (in_circle(get_point(basel->Dest()),
                                         get_point(basel->Org()),
                                         get_point(lcand->Dest()),
                                         get_point(lcand->Onext()->Dest()))) {
                                edge_type t = lcand->Onext();
                                remove_edge(lcand);
                                lcand = t;
                                const point_type& p1 = get_point(lcand->Dprev()->Org());
                                const point_type& p2 = get_point(lcand->Dnext()->Dnext()->Org());
                                trans[0] = p1[0];
                                trans[1] = p1[1];
                                trans[2] = p2[0];
                                trans[3] = p2[1];
                                comm.send(trans.data(), neighbor, 2 * point_type::dim);
                        }

                // Symmetrically, locate the first R point to be hit, and delete
                // R edges
                if (right_of(hit->first, basel))
                        while (in_circle(get_point(basel->Dest()),
                                         get_point(basel->Org()),
                                         hit->first,
                                         hit->second)) {
                                comm.recv(trans.data(), neighbor, 2 * point_type::dim);
                                hit = h.swap(hit, trans);
                        }

                // If both lcand and rcand are invalid, then basel is the upper
                // common tangent
                if (!right_of(lcand->Dest(), basel) &&
                    !right_of(hit->first, basel)) break;

                // The next cross edge is to be connected to either lcand.Dest
                // or rcand.Dest If both are valid, then choose the appropriate
                // one using the in_circle test.
                if (!right_of(lcand->Dest(), basel) ||
                    (right_of(hit->first, basel) &&
                     in_circle(
                             get_point(lcand->Dest()),
                             get_point(lcand->Org()),
                             get_point(basel->Org()),
                             hit->first))) {
                        // Add cross edge basel from rcand.Dest to basel.Dest
                        edge_type rcand = extend_edge(basel->Sym(), hit->first);
                        basel = connect_edges(rcand, basel->Sym());
                        ++hit;
                }
                else {
                        // Add cross edge basel from basel.Org to lcand.Dest
                        basel = connect_edges(basel->Sym(), lcand->Sym());
                }

        } // Merge loop
}

void Delaunay::merge_hull(const par::communicator& comm, unsigned int neighbor,
                          Hull& h, const edge_type& rdi, edge_type& rdo)
{
        Hull::iterator hit = h.begin();

        edge_type basel = extend_edge(rdi->Sym(), h.org());

        // Keep this here for old times' sake
        if (rdi->Org() == rdo->Org()) rdo = basel;

        std::vector<real> trans(4);
        // This is the merge loop
        while (true) {
#ifndef NDEBUG
                if (hit == h.end()) comm.abort("Prematurely hit end of hull", 1);
#endif
                // Locate the first L point (lcand.Dest) to be encountered by
                // the rising bubble, and delete L edges out of basel.Dest that
                // fail the circle test.
                if (right_of(hit->first, basel))
                        while (in_circle(get_point(basel->Dest()),
                                         get_point(basel->Org()),
                                         hit->first,
                                         hit->second)) {
                                comm.recv(trans.data(), neighbor, 4);
                                hit = h.swap(hit, trans);
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
                                const point_type& p1 = get_point(rcand->Dnext()->Org());
                                const point_type& p2 = get_point(rcand->Dprev()->Dprev()->Org());
                                trans[0] = p1[0];
                                trans[1] = p1[1];
                                trans[2] = p2[0];
                                trans[3] = p2[1];
                                comm.send(trans.data(), neighbor, 4);
                        }
                }

                // If both lcand and rcand are invalid, then basel is the upper
                // common tangent
                if (!right_of(hit->first, basel) &&
                    !right_of(rcand->Dest(), basel)) break;

                // The next cross edge is to be connected to either lcand.Dest
                // or rcand.Dest If both are valid, then choose the appropriate
                // one using the in_circle test.
                if (!right_of(hit->first, basel) ||
                    (right_of(rcand->Dest(), basel) &&
                     in_circle(
                             hit->first,
                             get_point(basel->Dest()),
                             get_point(rcand->Org()),
                             get_point(rcand->Dest())))) {
                        // Add cross edge basel from rcand.Dest to basel.Dest
                        basel = connect_edges(rcand, basel->Sym());
                }
                else {
                        // Add cross edge basel from basel.Org to lcand.Dest
                        edge_type lcand = extend_edge(basel, hit->first);
                        basel = connect_edges(basel->Sym(), lcand->Sym());
                        ++hit;
                }
        } // Merge loop
}


void Delaunay::create_merge_stack(const par::communicator& comm)
{
        // LIMITATION: This is very specific to a 1d merge of processors
        // in a "unary" layout
        const int my_rank = comm.rank();
        const int num_proc = comm.size();

        if (num_proc == 1) return;

        if (my_rank % 2) {
                const int neighbor = my_rank - 1;
                if (neighbor >= 0)
                        merge_stack.emplace(Neighbor({neighbor, Direction::Left}));
        }
        else {
                const int neighbor = my_rank + 1;
                if (neighbor < num_proc)
                        merge_stack.emplace(Neighbor({neighbor, Direction::Right}));
        }

        if (num_proc == 2) return;

        if (my_rank % 2) {
                const int neighbor = my_rank + 1;
                if (neighbor < num_proc)
                        merge_stack.emplace(Neighbor({neighbor, Direction::Right}));
        }
        else {
                const int neighbor = my_rank - 1;
                if (neighbor >= 0)
                        merge_stack.emplace(Neighbor({neighbor, Direction::Left}));
        }
}

void Delaunay::merge_proc(const par::communicator& comm, unsigned int neighbor,
                          Direction dir, edge_type& eo, edge_type& ei)
{
        // LIMITATION: This can only merge convex regions

        // ei seem to be CW when ei takes the rdi and CCW when ei takes the
        // place of ldi. While this is fine, it's easier to extract the hull and
        // do the merge if the face is exactly opposite. Inside the merge the
        // edges are flipped anyway before any merge occurs. Flip it here instead.
        ei = ei->Sym();

        // Create a subdivision of the part of the convex hull facing the
        // neighbor. Note that the extract_*_hull functions require that the
        // edge on the left is CCW and CW on the right
        std::vector<real> myhull;
        if (dir == Direction::Left)  myhull = extract_left_hull(ei);
        else if (dir == Direction::Right) myhull = extract_right_hull(ei);
        else comm.abort("Unknown merge direction", 1);

        const std::vector<real>::size_type myhullsize = myhull.size();

        comm.isend(&myhullsize, neighbor);

        std::vector<real>::size_type hullsize;
        par::request rs = comm.irecv(&hullsize, neighbor);
        rs.wait();

        // Next, send the hull
        par::request rhs = comm.isend(myhull.data(), neighbor, myhullsize);

        std::vector<real> hull(hullsize);
        par::request rr = comm.irecv(hull.data(), neighbor, hullsize);
        rr.wait();

        Hull h(hull);
        rhs.wait();
        if (dir == Direction::Right)
                merge_hull(comm, neighbor, eo, ei, h);
        else if (dir == Direction::Left)
                merge_hull(comm, neighbor, h, ei, eo);
        else
                comm.abort("Unknown merge direction", 1);
}


} // namespace cgl
