// PCFT.cpp : Defines the entry point for the application.
//

#include "tools/Timer.hpp"
#include "numerics/FTFactory.hpp"
#include "numerics/Preprocessor.hpp"
#include "DomainParameters.hpp"
#include <complex>
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

    const int N = 1'000;
    double r = 0.05;
    double sigma = 0.2;

    // Test usage of preprocessor
    DomainParameters pPackage 
    {
        1,
        N,

        -1.0,
        1.0,
        1.0,

        0.0,
        0.0
    };

    GFunction greensFunctionTransform(r, sigma, getDTau(pPackage.T, pPackage.M));
    Preprocessor preprocessor(FTFactory::instance(FTFactory::TransformType::FFT, N), pPackage);

    Timer t;
    t.start();
    IFourierTransformer::RealVec l2Err = preprocessor.execute(greensFunctionTransform);
    t.stop();
    auto timeTaken = t.duration();
}

