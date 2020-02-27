// PCFT.cpp : Defines the entry point for the application.
//

#include "tools/Timer.hpp"
#include "numerics/FTFactory.hpp"
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
    std::unique_ptr<IFourierTransformer> fftwTransformer = FTFactory::instance(FTFactory::TransformType::FFT, N);
    std::unique_ptr<IFourierTransformer> naiveTransformer = FTFactory::instance(FTFactory::TransformType::Naive, N);

    IFourierTransformer::ComplexVec fftwOutVec(N), naiveOutVec(N);
    IFourierTransformer::RealVec fftwRealInVec(N), naiveRealInVec(N);

    const double PI = 2*std::asin(1.0);
    int count = 0;
    const double x = PI;
    for (double& element : fftwRealInVec)
    {
        //const double x = count * PI / (N - 1.0);
        //element = 1/(x+1);
        element = x;
        count++;
    }

    count = 0;
    for (double& element : naiveRealInVec)
    {
        //const double PI = 2 * std::asin(1.0);
        //const double x = count * PI / (N - 1.0);
        //element = 1/(x+1);
        element = x;

        count++;
    }

    Timer t;
    t.start();
    fftwTransformer->fft(fftwRealInVec, fftwOutVec);
    fftwTransformer->ifft(fftwOutVec, fftwRealInVec);
    t.stop();
    Timer::milliseconds dur0 = t.duration();

    t.start();
    naiveTransformer->fft(naiveRealInVec, naiveOutVec);
    naiveTransformer->ifft(naiveOutVec, naiveRealInVec);
    t.stop();
    Timer::milliseconds dur1 = t.duration();

    
    std::vector<double> resultsReal(N);
    std::vector<double> resultsComplex(N);
    double resultsComplexL2Error = 0.0;
    double resultsRealL2Error = 0.0;
    for (size_t i = 0; i < N; i++)
    {
        resultsReal[i] = pow(fftwRealInVec[i] - naiveRealInVec[i], 2);
        resultsRealL2Error += resultsReal[i];

        resultsComplex[i] = pow(fftwOutVec[i].real() - naiveOutVec[i].real(), 2) + pow(fftwOutVec[i].imag() - naiveOutVec[i].imag(), 2);
        resultsComplexL2Error += resultsComplex[i];
    }
    resultsRealL2Error = sqrt(resultsRealL2Error/N);
    resultsComplexL2Error = sqrt(resultsComplexL2Error/N);

    std::cout << "Complex l2 error is: " << resultsComplexL2Error << std::endl;
    std::cout << "Real l2 error is: " << resultsRealL2Error << std::endl;

    // ************************ Test changing the transformer on the fly ****************************************
    N = 1000;

    fftwOutVec.resize(N);
    naiveOutVec.resize(N);
    fftwRealInVec.resize(N);
    naiveRealInVec.resize(N);

    count = 0;
    for (double& element : fftwRealInVec)
    {
        element = x;
        count++;
    }

    count = 0;
    for (double& element : naiveRealInVec)
    {
        element = x;
        count++;
    }

    t.start();
    fftwTransformer->fft(fftwRealInVec, fftwOutVec);
    fftwTransformer->ifft(fftwOutVec, fftwRealInVec);
    t.stop();
    Timer::milliseconds dur2 = t.duration();

    t.start();
    naiveTransformer->fft(naiveRealInVec, naiveOutVec);
    naiveTransformer->ifft(naiveOutVec, naiveRealInVec);
    t.stop();
    Timer::milliseconds dur3 = t.duration();
}

