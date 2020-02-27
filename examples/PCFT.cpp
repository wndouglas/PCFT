// PCFT.cpp : Defines the entry point for the application.
//

#include "tools/Timer.hpp"
#include "numerics/FTFactory.hpp"
#include "numerics/Preprocessor.hpp"
#include "ParameterPackage.hpp"
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
    ParameterPackage pPackage 
    {
        0,
        N,

        0.0,
        0.0,
        0.0,

        0.0,
        0.0
    };

    Preprocessor preprocessor(FTFactory::instance(FTFactory::TransformType::FFT, N), pPackage);

    std::vector<double> preprocessorOutputVec(N);
    preprocessor.execute(fftwRealInVec, preprocessorOutputVec);

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

