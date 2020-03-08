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

    int N = 100;
    std::unique_ptr<IFourierTransformer> fftwTransformer = FTFactory::instance(FTFactory::TransformType::Naive, N);

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

    double r = 0.05;
    double sigma = 0.2;
    GFunction greensFunctionTransform(r, sigma, DomainParameters::getDTau(pPackage.T, pPackage.M));
    Preprocessor preprocessor(FTFactory::instance(N), greensFunctionTransform, pPackage);

    std::vector<double> preprocessorOutputVec = preprocessor.execute();
}

