#include <iomanip>
#include <fstream>
#include "cgl/triangulation.hpp"

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
        std::ofstream fst(filePrefix + ".txt", std::ios::out);

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



Delaunay::Delaunay(PointSet& ps) {
        edge_type el, er;
        ps.sort(0);
        init_dc(ps, el, er);
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
                if (right_of(get_point(lcand->Dest()), basel))
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
                if (right_of(get_point(rcand->Dest()), basel))
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
                if (!right_of(get_point(lcand->Dest()), basel) &&
                    !right_of(get_point(rcand->Dest()), basel)) break;

                // The next cross edge is to be connected to either lcand.Dest
                // or rcand.Dest If both are valid, then choose the appropriate
                // one using the in_circle test.
                if (!right_of(get_point(lcand->Dest()), basel) ||
                    (right_of(get_point(rcand->Dest()), basel) &&
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

void Delaunay::init_dc(PointSet& ps, edge_type& el, edge_type& er)
{
        if (ps.size() == 2) {
                edge_type e = add_edge(ps.front(), ps.back());
                el = e;
                er = e->Sym();
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
                init_dc(ps, ldo, ldi);
                init_dc(psr, rdi, rdo);

                // Traverse the convex hulls of left and right point sets in preparation
                // to compute the base tangent
                while (left_of(get_point(rdi->Org()), ldi)) { ldi = ldi->Lnext(); }
                while (right_of(get_point(ldi->Org()), rdi)) { rdi = rdi->Rprev(); }

                // Merge the two triangulations
                merge(ldo, ldi, rdi, rdo);

                // Return [ldo, rdo]
                el = ldo;
                er = rdo;
        }
}

} // namespace cgl
