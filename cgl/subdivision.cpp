#include <stdexcept> // std::runtime_error

#include "cgl/subdivision.hpp"

namespace cgl
{
namespace detail
{

EdgeHashTable::EdgeHashTable() : table(), index(0), num(0) { }

EdgeHashTable::size_type EdgeHashTable::size() const
{
        return num;
}

EdgeHashTable::iterator EdgeHashTable::add(const EdgeHashTable::key_type& key,
                                                 const EdgeHashTable::value_type& val)
{
        const key_type vkey = get_vkey(key);
        id_type i = vkey[0], j = vkey[1];
#ifndef NDEBUG
        if (!valid_vkey(vkey)) {
                std::cout << "KEY = " << i << "," << j << std::endl;
                throw std::runtime_error("Invalid vkey");
        }
#endif

        ensure_row_exists(i);

        // Assumes the entry does not already exist
        table[i].push_back({j, val});

        // Increment number of entries
        num++;

        return --table[i].end();
}

EdgeHashTable::value_type& EdgeHashTable::get(const EdgeHashTable::key_type& key)
{
        key_type vkey = get_vkey(key);
#ifndef NDEBUG
        if (!valid_vkey(vkey)) throw std::runtime_error("Invalid vkey");
#endif
        iterator it = find(vkey);
#ifndef NDEBUG
        if (!valid(it)) throw std::runtime_error("Does not exist in hash");
#endif
        return it->e;

}

void EdgeHashTable::remove(const EdgeHashTable::key_type& key)
{
        key_type vkey = get_vkey(key);
#ifndef NDEBUG
        if (!valid_vkey(vkey)) throw std::runtime_error("Invalid vkey");
#endif
        iterator it = find(vkey);
#ifndef NDEBUG
        if (!valid(it)) throw std::runtime_error("Does not exist in hash");
#endif
        table[vkey[0]].erase(it);

        // Update number of edges
        num--;
}

bool EdgeHashTable::exists(const EdgeHashTable::key_type& key)
{
        key_type vkey = get_vkey(key);
#ifndef NDEBUG
        if (!valid_vkey(vkey)) throw std::runtime_error("Invalid vkey");
#endif
        return find(vkey) != end();
}

EdgeHashTable::iterator EdgeHashTable::begin_loop()
{
        index = 0;
        return table[index].begin();
}

bool EdgeHashTable::valid(EdgeHashTable::iterator it)
{
        return (index < table.size() && it != end());
}

EdgeHashTable::iterator EdgeHashTable::iterate(EdgeHashTable::iterator it)
{
#ifndef NDEBUG
        if (!valid(it)) throw std::runtime_error("Invalid iterator");
#endif

        ++it;
        if (it == table[index].end()) {
                while (true) {
                        index++;
                        if (index >= table.size()) break;
                        if (table[index].size() != 0) break;
                }
                // if (index < table.size() - 1) {
                //         while (table[index].size())
                // }
                // do {index++; } while (table[index].size() == 0 || index < table.size() - 1);
                // do { if (index < table.size() - 1) index++; } while (table[index].size() == 0);
                it = (index < table.size()) ? table[index].begin() : end();
        }

        return it;
}

EdgeHashTable::iterator EdgeHashTable::end()
{
        return table.back().end();
}


EdgeHashTable::key_type EdgeHashTable::get_vkey(const EdgeHashTable::key_type& key) const
{
        return (key[0] > key[1]) ?
                key_type({key[1], key[0]}) :
                key_type({key[0], key[1]});
}

bool EdgeHashTable::valid_vkey(const key_type& vkey) const
{
        if (vkey[0] < 0 || vkey[1] < 0) return false;
        if (vkey[0] >= vkey[1]) return false;
        return true;
}

EdgeHashTable::iterator EdgeHashTable::find(const EdgeHashTable::key_type& key)
{
        id_type i = key[0], j = key[1];

#ifndef NDEBUG
        if (table.size() <= static_cast<size_type>(i)) return end();
#endif

        for (iterator it=table[i].begin(); it!=table[i].end(); ++it)
                if (it->j == j) return it;

        return end();
}

void EdgeHashTable::ensure_row_exists(EdgeHashTable::size_type i)
{
        if (table.size() <= i)
                while (table.size() <= i)
                        table.push_back(std::list<EdgeEntry>());
}

} // namespace detail


Subdivision::~Subdivision()
{
        // EdgeHashTable does not delete the underlying edge
        edge_iterator it = edge.begin_loop();
        while (edge.valid(it)) {
                unlink_and_delete_edge(it->e);
                it = edge.iterate(it);
        }
}

Subdivision::size_type Subdivision::num_edges() const
{
        return edge.size();
}

Subdivision::size_type Subdivision::num_points() const
{
        return point.size();
}

Subdivision::edge_type Subdivision::add_edge(const Subdivision::point_type& org,
                                             const Subdivision::point_type& dest)
{
        id_type o = add_point(org);
        id_type d = add_point(dest);
        edge_type e = make_edge(o, d);
        const edge_key_type ekey = {o, d};
        edge.add(ekey, e);
        startedge = e;
        return e;
}

Subdivision::edge_type Subdivision::extend_edge(const Subdivision::edge_type& e0,
                                                const Subdivision::point_type& dest)
{
        // Add an edge to the existing subdivision with a new point
        id_type d = add_point(dest);
        edge_type e = make_edge(e0->Dest(), d);
        const edge_key_type ekey = {e0->Dest(), d};
        edge.add(ekey, e);
        splice(e0->Sym(), e);
        startedge = e;
        return e;
}

Subdivision::edge_type Subdivision::connect_edges(const Subdivision::edge_type& a,
                                                  const Subdivision::edge_type& b)
{
        const edge_key_type ekey = {a->Dest(), b->Org()};
        edge_type e = make_edge(a->Dest(), b->Org());
        edge.add(ekey, e);
        splice(e, a->Lnext());
        splice(e->Sym(), b);
        startedge = e;
        return e;
}

void Subdivision::remove_edge(const Subdivision::edge_type& e)
{
        // Remove it from the table
        const edge_key_type ekey = {e->Org(), e->Dest()};
        edge.remove(ekey);

        // Pick a new startedge... Hopefully this one is valid
        if (startedge == e) startedge = e->Lnext();

        // Deallocate e
        unlink_and_delete_edge(e);
}

Subdivision::edge_type Subdivision::last_edge() const
{
        return startedge;
}

Subdivision::edge_type Subdivision::locate(const point_type& p) const
{
        edge_type e = startedge;
        while (true) {
                if (p == get_point(e->Org()) || p == get_point(e->Dest()))
                        return e;
                else if (right_of(p, e))
                        e = e->Sym();
                else if (!right_of(p, e->Onext()))
                        e = e->Onext();
                else if (!right_of(p, e->Dprev()))
                        e = e->Dprev();
                else return e;
        }
}

const Subdivision::point_type& Subdivision::get_point(size_type i) const
{
        return point[i];
}


bool Subdivision::right_of(const point_type& x, const edge_type& e) const
{
        return ccw(x, get_point(e->Dest()), get_point(e->Org()));
}


bool Subdivision::left_of(const point_type& x, const edge_type& e) const
{
        return ccw(x, get_point(e->Org()), get_point(e->Dest()));
}

int Subdivision::add_point(const Subdivision::point_type& p)
{
        point.push_back(p);
        return point.size() - 1;
}

void Subdivision::unlink_and_delete_edge(const Subdivision::edge_type& e)
{
        delete_edge(e);
}

} // namespace cgl
