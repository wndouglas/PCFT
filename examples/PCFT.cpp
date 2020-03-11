// PCFT.cpp : Defines the entry point for the application.
//

#include "tools/Timer.hpp"
#include "numerics/FTFactory.hpp"
#include "numerics/Preprocessor.hpp"
#include "DomainParameters.hpp"
#include "PCFTExecutor.hpp"
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

    // // Parameters
    // double K = 100.0;
    double r = 0.05;
    double sigma = 0.2;
    double T = 1.0;

    double xMin = std::log(0.1);
    double xMax = std::log(200);

    // Test usage of preprocessor
    DomainParameters pPackage 
    {
        M,
        N,

        xMin,
        xMax,
        T,

        1e-8,
        1e-8
    };
    
    PCFTExecutor executor(
        FTFactory::instance(N*2),
        std::make_unique<Preprocessor>(GFunction(r, sigma, DomainParameters::getDTau(pPackage.T, pPackage.M))),
        pPackage);

    OutputPage outputPackage;
    executor.execute(outputPackage);
}

