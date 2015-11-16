#ifndef QUADEDGE_HH
#define QUADEDGE_HH

#include "Geom2d.hh"

class QuadEdge;

class Edge {
        friend QuadEdge;
        friend void splice(Edge*, Edge*);

private:
        int num;
        Edge    *next;
        Point2d *data;
public:
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
        // Origin
        inline Point2d* Org() { return data; }

        // Destination
        inline Point2d* Dest() { return Sym()->data; }

        inline const Point2d& Org2d() const { return *data; }

        inline const Point2d& Dest2d() const {
                return (num < 2) ? *((this + 2)->data) : *((this - 2)->data);
        }


        inline void end_points(Point2d* o, Point2d* d) {
                data = o; Sym()->data =
                        d;
        }


        inline QuadEdge* Qedge() {
                return reinterpret_cast<QuadEdge*>(this -
                                                   num);
        }


};

class QuadEdge {
        friend Edge* make_edge();

private:
        Edge e[4];
public:
        QuadEdge();
};


/*********************** Basic Topological Operators ****************/
Edge* make_edge();

void splice(Edge* a, Edge* b);

void delete_edge(Edge* e);

#endif
