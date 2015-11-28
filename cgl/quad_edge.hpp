#ifndef CGL_QUAD_EDGE_HH
#define CGL_QUAD_EDGE_HH

#include "cgl/geom2d.hpp"

namespace cgl
{

class QuadEdge;

class Edge {
public:
        typedef long int id_type;

        /************** Edge Algebra ****************/

        // Return the dual of the current edge, directed from its right to its left.
        inline Edge* Rot() { return (num < 3) ? this + 1 : this - 3; }

        // Return the dual of the current edge, directed from its left to its right.
        inline Edge* invRot() { return (num > 0) ? this - 1 : this + 3; }

        // Return the edge from the destination to the origin of the current edge.
        inline Edge* Sym() { return (num < 2) ? this + 2 : this - 2; }

        // Return the next ccw edge around (from) the origin of the current edge.
        inline Edge* Onext() { return next; }

        // Return the next cw edge around (from) the origin of the current edge.
        inline Edge* Oprev() { return Rot()->Onext()->Rot(); }

        // Return the next ccw edge around (into) the destination of the current edge.
        inline Edge* Dnext() { return Sym()->Onext()->Sym(); }

        // Return the next cw edge around (into) the destination of the current edge.
        inline Edge* Dprev() { return invRot()->Onext()->invRot(); }

        // Return the ccw edge around the left face following the current edge.
        inline Edge* Lnext() { return invRot()->Onext()->Rot(); }

        // Return the ccw edge around the left face before the current edge.
        inline Edge* Lprev() { return Onext()->Sym(); }

        // Return the edge around the right face ccw following the current edge.
        inline Edge* Rnext() { return Rot()->Onext()->invRot(); }

        // Return the edge around the right face ccw before the current edge.
        inline Edge* Rprev() { return Sym()->Onext(); }


        /************** Access to data pointers ****************/
        inline id_type Org() { return id; }

        inline id_type Dest() { return Sym()->id; }

        inline const QuadEdge* Qedge() const {
                return reinterpret_cast<const QuadEdge*>(this - num);
        }

        inline void set_id(id_type o, id_type d) {
                id = o;
                Sym()->id = d;
        }

private:
        friend QuadEdge;
        friend void splice(Edge*, Edge*);
        friend Edge* make_edge(id_type, id_type);

        unsigned int num;
        Edge* next;
        id_type id;
};

class QuadEdge {
public:
        QuadEdge();
private:
        static constexpr unsigned int quad = 4;

        friend Edge* make_edge();
        friend Edge* make_edge(Edge::id_type, Edge::id_type);

        Edge e[quad];
};


/*********************** Basic Topological Operators ****************/
Edge* make_edge();
Edge* make_edge(Edge::id_type, Edge::id_type);

void splice(Edge* a, Edge* b);

void delete_edge(Edge* e);

} // namespace cgl

#endif
