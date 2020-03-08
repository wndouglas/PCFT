// PCFT.cpp : Defines the entry point for the application.
//

#include "tools/Timer.hpp"
#include "numerics/FTFactory.hpp"
#include "numerics/Preprocessor.hpp"
#include "DomainParameters.hpp"
#include <complex>
#include <cmath>
#include <thread>
#include <iostream>

using namespace PCFT;

void do_transform();

int main()
{
    using namespace tools;
    Timer t;
    t.start();
    do_transform();
    t.stop();
    Timer::milliseconds dur = t.duration();

	return 0;
}

void do_transform()
{
    using namespace numerics;
    using namespace tools;

    int N = 100;
    int M = 1;

    // Parameters
    double K = 100.0;
    double r = 0.05;
    double sigma = 0.2;
    double T = 1.0;

    double xMin = std::log(0.1);
    double xMax = std::log(200);
    double xMinExtended = xMin - 0.5*(xMax - xMin);
    double xMaxExtended = xMax + 0.5*(xMax - xMin);
    double PExtended = xMaxExtended - xMinExtended;
    int NExtended = 2*N;

    // Test usage of preprocessor
    DomainParameters pPackage 
    {
        M,
        NExtended,

        xMinExtended,
        xMaxExtended,
        T,

        1e-8,
        1e-8
    };

    GFunction greensFunctionTransform(r, sigma, DomainParameters::getDTau(pPackage.T, pPackage.M));
    Preprocessor preprocessor(FTFactory::instance(NExtended), greensFunctionTransform, pPackage);
    IFourierTransformer::ComplexVec GOut = preprocessor.execute();

    // For now we try a European call
    double dx = DomainParameters::getDx(xMax, xMin, N);
    std::vector<double> SGrid(N);
    std::vector<double> SGridExtended(NExtended);
    std::vector<double> vGrid(N);
    std::vector<double> vGridExtended(NExtended);

    double x = xMin;
    for(int i = 0; i < N; i++)
    {
        double S = std::exp(x);
        SGrid[i] = S;
        vGrid[i] = std::max(0.0, S - K);
        x += dx;
    }

    x = xMinExtended;
    for(int i = 0; i < NExtended; i++)
    {
        double S = std::exp(x);
        vGridExtended[i] = std::max(0.0, S - K);
        SGridExtended[i] = S;
        x += dx;
    }
    
    std::unique_ptr<IFourierTransformer> transformer = FTFactory::instance(NExtended);
    IFourierTransformer::ComplexVec vTempGrid(NExtended);
    for(int i = 0; i < NExtended; i++)
    {
        vTempGrid[i].real(vGridExtended[i]);
    }

    for(int m = 0; m < M; m++)
    {
        transformer->shiftedFft(vTempGrid, vTempGrid);
        for(int i = 0; i < NExtended; i++)
        {
            vTempGrid[i] = vTempGrid[i] * GOut[i];
        }
        transformer->shiftedIfft(vTempGrid, vTempGrid);

        for(int i = 0; i < NExtended; i++)
        {
            if(i < N/2)
            {
                vTempGrid[i] = vTempGrid[N/2 + 1];
            }
            else if(i > 3*N/2)
            {
                vTempGrid[i] = SGridExtended[i] - K*exp(-r*T);
            }
        }
    }
    std::cout << "Done";
}

