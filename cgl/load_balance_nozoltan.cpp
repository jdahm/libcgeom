#include "cgl/load_balance.hpp"

namespace cgl {

typedef Point2d point_type;
typedef std::list<point_type> list_type;

void lb_zoltan(PSTopology top, LBMethod method, unsigned int dimen, list_type& pl)
{
        std::cerr << "No Zoltan detected." << std::endl;
}

}
