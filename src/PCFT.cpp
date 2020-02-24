// PCFT.cpp : Defines the entry point for the application.
//

#include "tools/Timer.hpp"
#include "numerics/IFourierTransformer.hpp"
#include "numerics/FFTWTransformer.hpp"
#include "numerics/NaiveTransformer.hpp"
#include <complex>
#include <thread>

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

    const int N = 2;
    std::unique_ptr<IFourierTransformer> fftwTransformer(new FFTWTransformer);
    std::unique_ptr <IFourierTransformer> naiveTransformer(new NaiveTransformer);

    IFourierTransformer::ComplexVec fftwInVec(N), fftwOutVec(N), naiveInVec(N), naiveOutVec(N);

    int count = 0;
    for (std::complex<double> &element : fftwInVec)
    {
        const double PI = 2*std::asin(1.0);
        const double x = count * PI / (N - 1.0);

        element = {0.0,  1/(x+1) };
         
        count++;
    }

    count = 0;
    for (std::complex<double>& element : naiveInVec)
    {
        const double PI = 2 * std::asin(1.0);
        const double x = count * PI / (N - 1.0);

        element = {0.0,  1/(x+1) };

        count++;
    }

    Timer t;
    t.start();
    fftwTransformer->fft(fftwInVec, fftwOutVec);
    fftwTransformer->ifft(fftwOutVec, fftwInVec);
    t.stop();
    Timer::milliseconds dur1 = t.duration();

    t.start();
    naiveTransformer->fft(naiveInVec, naiveOutVec);
    naiveTransformer->ifft(naiveOutVec, naiveInVec);
    t.stop();
    Timer::milliseconds dur2 = t.duration();

    std::vector<double> results(N);
    for (size_t i = 0; i < N; i++)
    {
        double realPart = fftwInVec[i].real() - naiveInVec[i].real();
        double imagPart = fftwInVec[i].imag() - naiveInVec[i].imag();
        results[i] = std::sqrt(realPart * realPart + imagPart * imagPart);
    }
}

