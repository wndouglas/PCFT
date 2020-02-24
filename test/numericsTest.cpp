#include "catch.hpp"
#include "../src/numerics/fourierTransform/FFTWTransformer.hpp"
#include "../src/numerics/fourierTransform/NaiveTransformer.hpp"
#include "../src/numerics/fourierTransform/FourierTransformerFactory.hpp"
#include "../src/tools/timer/TimerFactory.hpp"
#include <memory>

using namespace PCFT::numerics;
using namespace PCFT::tools;
typedef std::unique_ptr<IFourierTransformer> pTransformer; 

void arrange_transforms(pTransformer& fftwTransformer,
pTransformer& naiveTransformer,
IFourierTransformer::ComplexVec& fftwInVec,IFourierTransformer::ComplexVec& fftwOutVec,
IFourierTransformer::ComplexVec& naiveInVec, IFourierTransformer::ComplexVec& naiveOutVec,
size_t N)
{
    //std::unique_ptr<ITimer> pTimer = TimerFactory::instance();

    typedef std::complex<double> CDouble;
    fftwTransformer.reset(new FFTWTransformer());
    naiveTransformer.reset(new NaiveTransformer());

    fftwInVec.resize(N);
    fftwOutVec.resize(N);
    naiveInVec.resize(N);
    naiveOutVec.resize(N);
    int count = 0;
    for (CDouble &element : fftwInVec)
    {
         const double PI = 2*std::asin(1.0);
         const double x = count * PI / (N - 1.0);

        element = {0.0,  1/(x+1) };
         
        count++;
    }

    count = 0;
    for (CDouble& element : naiveInVec)
    {
        const double PI = 2 * std::asin(1.0);
        const double x = count * PI / (N - 1.0);

        element = {0.0,  1/(x+1) };

        count++;
    }
}

TEST_CASE( "Test FFT Fourier transform", "[transform]" )
{
    // Arrange
    pTransformer fftwTransformer, naiveTransformer;
    IFourierTransformer::ComplexVec fftwInVec, fftwOutVec;
    IFourierTransformer::ComplexVec naiveInVec, naiveOutVec;
    size_t N = 100;
    arrange_transforms(fftwTransformer, naiveTransformer, fftwInVec, fftwOutVec, 
        naiveInVec, naiveOutVec, N);

    // Act
    fftwTransformer->fft(fftwInVec, fftwOutVec);
    fftwTransformer->ifft(fftwOutVec, fftwInVec);

    naiveTransformer->fft(naiveInVec, naiveOutVec);
    naiveTransformer->ifft(naiveOutVec, naiveInVec);

    // Assert
    bool toleranceMet = true;
    for (size_t i = 0; i < N; i++)
    {
        double realPart = fftwInVec[i].real() - naiveInVec[i].real();
        double imagPart = fftwInVec[i].imag() - naiveInVec[i].imag();
        double result = std::sqrt(realPart * realPart + imagPart * imagPart);
        if(fabs(result) >= 1e-8)
        {
            toleranceMet = false;
            break;
        }
    }
    REQUIRE(toleranceMet);
}