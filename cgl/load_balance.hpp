#ifndef CGL_LOAD_BALANCE_HPP
#define CGL_LOAD_BALANCE_HPP

#include <vector>

#include "cgl/cgl.hpp"
#include "cgl/geom2d.hpp"

namespace cgl
{

enum class LBMethod {S2A2, RCB};

enum class PSTopology {Unary, Binary, Zoltan};

void balance_set(PSTopology, LBMethod, unsigned int, std::list<Point2d>&);

} // namespace cgl

#endif
