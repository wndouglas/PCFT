// PCFT.cpp : Defines the entry point for the application.
//

#include "PCFT.hpp"
#include "tools/Timer.hpp"
#include "numerics/IFourierTransformer.hpp"
#include "numerics/FFTWTransformer.hpp"
#include "numerics/NaiveTransformer.hpp"
#include <complex>

void do_transform();

int main()
{
    PCFT::tools::Timer<> t;
    do_transform();
    auto out = t.tick();

	return 0;
}

void do_transform()
{
    using namespace PCFT::numerics;

    const int N = 10000;
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

    PCFT::tools::Timer<> t;
    fftwTransformer->fft(fftwInVec, fftwOutVec);
    fftwTransformer->ifft(fftwOutVec, fftwInVec);
    auto fftwOut = t.tick();

    PCFT::tools::Timer<> tNew;
    naiveTransformer->fft(naiveInVec, naiveOutVec);
    naiveTransformer->ifft(naiveOutVec, naiveInVec);
    auto naiveOut = tNew.tick();

    std::vector<double> results(N);
    for (size_t i = 0; i < N; i++)
    {
        double realPart = fftwInVec[i].real() - naiveInVec[i].real();
        double imagPart = fftwInVec[i].imag() - naiveInVec[i].imag();
        results[i] = std::sqrt(realPart * realPart + imagPart * imagPart);
    }
}

