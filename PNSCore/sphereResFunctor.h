#ifndef __SPHERERESFUNCTOR_H__
#define __SPHERERESFUNCTOR_H__
#include "GenericFunctor.h"
//TODO: replace cmath header with CommonType.h and replace every M_PI with PI
#include <cmath> // To be replace with CommonType.h

class sphereResFunctor : public GenericFunctor<double> {
    public:
        typedef enum mode_t{
            SEQTEST,
            SMALL,
            GREAT
        } mode_t;
        mode_t m_mode;
        sphereResFunctor(int sizeX, int sizeY, const ValueType& y, int mode): GenericFunctor<double>(sizeX, sizeY, y), m_mode( static_cast<mode_t>(mode) ) {}; 
        double operator() (const VectorXd& x, VectorXd& fvec) const {
            assert(x.size()==inputs());
            assert(fvec.size()==values());

            MatrixXd inputMatrix(inputs(), values());
            MatrixXd auxMatrix(inputs(), values());
            double radius(0);

            inputMatrix = x.replicate(1, values());
            auxMatrix   = ( (m_y.array() - inputMatrix.array()).square() ).matrix();

            VectorXd di = (auxMatrix.colwise().sum()).array().sqrt().matrix();
            if(m_mode != GREAT) {
                radius = di.array().mean();
            }
            else {
                radius = M_PI / 2; // To be replaced with PI
            }

            fvec = ( di.array() - radius ).matrix();
            return radius;
        }
        int df(const VectorXd& x, MatrixXd& fjac) const {

            assert(fjac.rows() == values() && fjac.cols() == inputs());
            VectorXd rC(values());
            (*this)(x,rC);
            double TOL = 1e-10;
            double dx = 0.25*TOL;
            VectorXd rd(values());
            VectorXd xd(inputs());

            for(int i=0;i<inputs();++i) { 
                xd = x;
                xd(i) += dx;
                rd.setZero();
                (*this)(xd, rd);
                fjac.col(i) = (rd.array() - rC.array())/dx;
            }
            return 0;
        }
};

#endif //__SPHERERESFUNCTOR_H__
