#include "numerics/FTFactory.hpp"
#include "catch.hpp"

using namespace PCFT::numerics;
typedef std::unique_ptr<IFourierTransformer> pTransformer; 

void arrange_transforms(pTransformer& fftwTransformer,
pTransformer& naiveTransformer,
IFourierTransformer::RealVec& fftwInVec,IFourierTransformer::ComplexVec& fftwOutVec,
IFourierTransformer::RealVec& naiveInVec, IFourierTransformer::ComplexVec& naiveOutVec,
size_t N)
{
    fftwTransformer = FTFactory::instance(FTFactory::TransformType::FFT, N);
    naiveTransformer = FTFactory::instance(FTFactory::TransformType::FFT, N);

    fftwInVec.resize(N);
    fftwOutVec.resize(N);
    naiveInVec.resize(N);
    naiveOutVec.resize(N);
    int count = 0;
    for (double& element : fftwInVec)
    {
        const double PI = 2*std::asin(1.0);
        const double x = count * PI / (N - 1);

        element = std::sin(x);
         
        count++;
    }

    count = 0;
    for (double& element : naiveInVec)
    {
        const double PI = 2 * std::asin(1.0);
        const double x = count * PI / (N - 1.0);

        element = std::sin(x);

        count++;
    }
}

TEST_CASE( "Fourier transform followed by its reverse should be identity operation", "[transform]" )
{
    // Arrange
    pTransformer fftwTransformer, naiveTransformer;
    IFourierTransformer::RealVec fftwInVec, naiveInVec;
    IFourierTransformer::ComplexVec fftwOutVec, naiveOutVec;
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
        double result = fftwInVec[i] - naiveInVec[i];
        if(fabs(result) >= 1e-8)
        {
            toleranceMet = false;
            break;
        }
    }
    REQUIRE(toleranceMet);
}