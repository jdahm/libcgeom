#include "cgl/quad_edge.hpp"

namespace cgl
{

QuadEdge::QuadEdge() {
        e[0].num = 0, e[1].num = 1, e[2].num = 2, e[3].num = 3;
        e[0].next = &(e[0]); e[1].next = &(e[3]);
        e[2].next = &(e[2]); e[3].next = &(e[1]);
}


/*********************** Basic Topological Operators ****************/

// Factory function for edges
Edge* make_edge() {
        QuadEdge *ql = new QuadEdge();
        return ql->e;
}

Edge* make_edge(Edge::id_type org, Edge::id_type dest) {
        Edge *e = make_edge();
        e->set_id(org, dest);
        return e;
}


// This operator affects the two edge rings around the origins of a and b,
// and, independently, the two edge rings around the left faces of a and b.
// In each case, (i) if the two rings are distinct, Splice will combine
// them into one; (ii) if the two are the same ring, Splice will break it
// into two separate pieces.
// Thus, Splice can be used both to attach the two edges together, and
// to break them apart. See Guibas and Stolfi (1985) p.96 for more details
// and illustrations.
// For the rings a Org and b Org, the cuts will occur immediately
// after a and b (in counterclockwise order); for the rings aLeft and
// bLeft, the cut will occur immediately before a Rot and b Rot.
void splice(Edge* a, Edge* b) {
        Edge* alpha = a->Onext()->Rot();
        Edge* beta  = b->Onext()->Rot();
        Edge* t1 = b->Onext();
        Edge* t2 = a->Onext();
        Edge* t3 = beta->Onext();
        Edge* t4 = alpha->Onext();
        a->next = t1;
        b->next = t2;
        alpha->next = t3;
        beta->next = t4;
}


void delete_edge(Edge* e) {
        splice(e, e->Oprev());
        splice(e->Sym(), e->Sym()->Oprev());
        delete e->Qedge();
}

} // namespace cgl
