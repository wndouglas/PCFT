#include "numerics/FTFactory.hpp"
#include "catch.hpp"

using namespace PCFT::numerics;
typedef std::unique_ptr<IFourierTransformer> pTransformer; 

void arrange_transforms(pTransformer& fftwTransformer,
pTransformer& naiveTransformer,
IFourierTransformer::ComplexVec& fftwInVec,IFourierTransformer::ComplexVec& fftwOutVec,
IFourierTransformer::ComplexVec& naiveInVec, IFourierTransformer::ComplexVec& naiveOutVec,
size_t N)
{
    typedef std::complex<double> CDouble;
    fftwTransformer = FTFactory::instance(FTFactory::TransformType::FFT);
    naiveTransformer = FTFactory::instance(FTFactory::TransformType::FFT);

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
    size_t N = 10000;
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