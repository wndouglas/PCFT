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

    int N = 100'000;
    std::unique_ptr<IFourierTransformer> fftwTransformer = FTFactory::instance(FTFactory::TransformType::FFT, N);

    IFourierTransformer::ComplexVec fftwOutVec(N);
    IFourierTransformer::RealVec fftwRealInVec(N);

    const double PI = 2*std::asin(1.0);
    int count = 0;
    for (double& element : fftwRealInVec)
    {
        double x = count * PI / (N - 1.0);
        element = 1/(x+1);
        count++;
    }

    Timer t;
    t.start();
    fftwTransformer->fft(fftwRealInVec, fftwOutVec);
    fftwTransformer->ifft(fftwOutVec, fftwRealInVec);
    t.stop();
    Timer::milliseconds dur0 = t.duration();

    // Test usage of preprocessor
    DomainParameters pPackage 
    {
        1,
        N,

        0.0,
        0.0,
        1.0,

        0.0,
        0.0
    };

    double r = 0.05;
    double sigma = 0.2;
    GFunction greensFunctionTransform(r, sigma, DomainParameters::getDTau(pPackage.T, pPackage.M));
    Preprocessor preprocessor(FTFactory::instance(N), greensFunctionTransform, pPackage);

    std::vector<double> preprocessorOutputVec = preprocessor.execute();

    std::vector<double> resultsVec(N);
    double l2Error = 0.0;
    for(size_t i = 0; i < N; i++)
    {
        resultsVec[i] = pow(fftwRealInVec[i] - preprocessorOutputVec[i], 2);
        l2Error += resultsVec[i];
    }
    l2Error = sqrt(l2Error);
    std::cout << "Results of preprocessor - l2 error: " << l2Error << std::endl;
}

