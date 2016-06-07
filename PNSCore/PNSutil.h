#ifndef __PNSUTIL_H__
#define __PNSUTIL_H__
#include <Eigen/Dense>
#include <cmath>
using std::floor;
using std::atan2;

namespace statismo {
    //const double PI	=	3.14159265358979323846; // TODO: comment this out later when it is integrated to statismo
    typedef enum itype_t {
        SEQTEST,
        SMALL,
        GREAT
    } itype;

    inline double modBy2PI( const double& x ) { 
        return ( x - (2*PI)*floor( x / (2*PI) ) ); 
    } ;

    template <typename Scalar> struct MakeAtan2Op{
        EIGEN_EMPTY_STRUCT_CTOR(MakeAtan2Op);
        Scalar operator() (const Scalar& a, const Scalar& b) const { return atan2(a,b); }
    };
} // namespace statismo


#endif //__PNSUTIL_H__
