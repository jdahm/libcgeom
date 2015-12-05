#ifndef CGL_SUBDIVISION_HPP
#define CGL_SUBDIVISION_HPP

#include <cstddef> // std::size_t
#include <vector>
#include <list>
#include <array>
#include <string>

#include "cgl/quad_edge.hpp"
#include "cgl/geom2d.hpp"

namespace cgl
{
namespace detail
{

class EdgeHashTable {
public:
        typedef std::size_t size_type;
        typedef Edge::id_type id_type;
        typedef Edge* value_type;
        static constexpr size_type nsort = 2;

        typedef std::array<id_type, nsort> key_type;

private:
        struct EdgeEntry {
                id_type j;
                value_type e;
        };
        typedef std::vector< std::list<EdgeEntry> > table_type;

public:
        typedef typename std::list<EdgeEntry>::iterator iterator;

        EdgeHashTable();

        size_type size() const;

        iterator add(const key_type&, const value_type&);

        value_type& get(const key_type&);

        void remove(const key_type&);

        bool exists(const key_type&);

        iterator begin_loop();

        bool valid(iterator);

        iterator iterate(iterator);

        iterator end();

private:
        key_type get_vkey(const key_type&) const;

        bool valid_vkey(const key_type&) const;

        iterator find(const key_type&);

        void ensure_row_exists(size_type);

        table_type table;
        size_type index;
        size_type num;
};

} // namespace detail


class Subdivision {
public:
        typedef Point2d point_type;
        typedef detail::EdgeHashTable edge_hash_type;
        typedef typename edge_hash_type::id_type id_type;
        typedef typename edge_hash_type::value_type edge_type;
        typedef typename edge_hash_type::key_type edge_key_type;
        typedef typename edge_hash_type::iterator edge_iterator;
        typedef std::vector<point_type> point_container_type;
        typedef std::vector<point_type> point_index_container_type;
        typedef typename point_container_type::size_type size_type;
        typedef typename point_container_type::const_iterator point_iterator;

        friend class EdgeHashTable;

        ~Subdivision();

        // Delete some operators (for now)
        Subdivision operator=(const Subdivision& other) = delete;

        size_type num_edges() const;
        size_type num_points() const;

        edge_type add_edge(const point_type&, const point_type&);
        edge_type extend_edge(const edge_type&, const point_type&);
        edge_type connect_edges(const edge_type&, const edge_type&);

        void remove_edge(const edge_type&);

        edge_type last_edge() const;

        edge_type locate(const point_type&) const;

        const point_type& get_point(size_type i) const;

        bool right_of(const point_type&, const edge_type&) const;
        bool right_of(size_type, const edge_type&) const;

        bool left_of(const point_type&, const edge_type&) const;
        bool left_of(size_type, const edge_type&) const;

        Subdivision extract_hull(edge_type, const point_type&) const;

        void write_txt(const std::string&);

protected:
        int add_point(const point_type&);

        void unlink_and_delete_edge(const edge_type&);

        point_container_type point;
        edge_hash_type edge;
        edge_type startedge;
        point_index_container_type global_id;
};

} // namespace cgl

#endif
