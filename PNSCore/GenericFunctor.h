#ifndef __GENERICFUNCTOR_H__
#define __GENERICFUNCTOR_H__
#include <iostream>
#include <unsupported/Eigen/NonLinearOptimization>
using Eigen::MatrixXd;
using Eigen::VectorXd;

template <typename _Scalar, int NX=Eigen::Dynamic, int NY=Eigen::Dynamic>
class GenericFunctor {
    public:
        typedef _Scalar Scalar;
        enum{
            InputsAtCompileTime = NX,
            ValuesAtCompileTime = NY
        };

        typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
        typedef Eigen::Matrix<Scalar,InputsAtCompileTime, ValuesAtCompileTime> ValueType;
        //        typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;
        int m_inputs, m_values;

        //VectorXd m_y;
        ValueType m_y;
        GenericFunctor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
        GenericFunctor(int inputs, int values, const ValueType& y) : m_inputs(inputs), m_values(values), m_y(y) {}

        int inputs() const { return m_inputs; }
        int values() const { return m_values; }
        // you should define that in the subclass :
        //  void operator() (const InputType& x, ValueType* v, JacobianType* _j=0) const;
};



#endif /*__GENERICFUNCTOR_H__*/
