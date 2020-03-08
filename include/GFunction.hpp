#ifndef G_FUNCTION_HPP
#define G_FUNCTION_HPP

#include "numerics/Constants.hpp"
#include <complex>

namespace PCFT
{
    class GFunction
    {
    public:
        GFunction(double r, double sigma, double dTau) :
            mDTau(dTau),
            mQuadraticTerm{-sigma*sigma*2*numerics::PI*numerics::PI, 0.0},
            mLinearTerm{0.0, (2*r - sigma*sigma)*numerics::PI},
            mConstantTerm{-r, 0.0}
        { }

        std::complex<double> operator()(double omega) const
        {
            std::complex<double> omegaSquared = omega*omega;
            std::complex<double> logTerm = (omegaSquared*mQuadraticTerm + omega*mLinearTerm + mConstantTerm)*mDTau;
            
            return std::exp(logTerm);
        }

    private:
        double mDTau;
        std::complex<double> mQuadraticTerm;
        std::complex<double> mLinearTerm;
        std::complex<double> mConstantTerm;
    };
}
#endif // !G_FUNCTION_HPP