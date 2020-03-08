#ifndef DOMAIN_PARAMETERS_HPP
#define DOMAIN_PARAMETERS_HPP

namespace PCFT
{
    static double getP(double xMax, double xMin)
    {
        return xMax - xMin;
    }

    static double getDx(double xMax, double xMin, int N)
    {
        return getP(xMax, xMin)/N;
    }

    static double getDTau(double T, int M)
    {
        return T/M;
    }

    struct DomainParameters
    {
        int M;
        int N;

        double xMin;
        double xMax;
        double T;

        double epsilon1;
        double epsilon2;
    };
}


#endif // !DOMAIN_PARAMETERS_HPP