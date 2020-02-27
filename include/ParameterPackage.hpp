#ifndef PARAMETER_PACKAGE_HPP
#define PARAMETER_PACKAGE_HPP

namespace PCFT
{
    class ParameterPackage
    {
    public:
        int M;
        int N;

        double xMax;
        double xMin;
        double T;

        double epsilon1;
        double epsilon2;

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
        };
}


#endif // !PARAMETER_PACKAGE_HPP