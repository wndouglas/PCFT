#include "NaiveTransformer.hpp"

using namespace PCFT::numerics;

typedef IFourierTransformer::ComplexVec CVec;
typedef IFourierTransformer::RealVec RVec;

namespace
{
    const double PI = atan(1.0) * 4.0;
    CVec complexConjugate(const CVec& cVec)
    {
        const size_t N = cVec.size();
        CVec newVec(N);
        for(size_t i = 0; i < N; i++)
        {
            newVec[i].imag(-cVec[i].imag());
            newVec[i].real(cVec[i].real());
        }
        return newVec;
    }

    void complex_fft(const CVec& inputVector, CVec& outputVector)
    {
        const size_t N = inputVector.size();

        for (size_t k = 0; k < N; k++)
        {
            std::complex<double>& outputElement = outputVector[k];
            outputElement = { 0.0, 0.0 };
            for (size_t j = 0; j < N; j++)
            {
                const std::complex<double>& inputElement = inputVector[j];
                double x1 = inputElement.real();
                double y1 = inputElement.imag();

                double x2 = std::cos(-2 * PI * j * k / N);
                double y2 = std::sin(-2 * PI * j * k / N);

                outputElement.real(outputElement.real() + x1 * x2 - y1 * y2);
                outputElement.imag(outputElement.imag() + x1 * y2 + x2 * y1);
            }
            outputElement.real(outputElement.real() / std::sqrt(N));
            outputElement.imag(outputElement.imag() / std::sqrt(N));
        }
    }
}

void NaiveTransformer::fft(const CVec& inputVector, CVec& outputVector) const
{
    const size_t N = inputVector.size();
    if (outputVector.size() != N)
    {
        throw std::runtime_error("Invalid input vectors");
    }
    complex_fft(inputVector, outputVector);
}

void NaiveTransformer::ifft(const CVec& inputVector, CVec& outputVector) const
{
    const size_t N = inputVector.size();
    if (outputVector.size() != N)
    {
        throw std::runtime_error("Invalid input vectors");
    }

    CVec conjugateInputVec = complexConjugate(inputVector);

    double small_val = 1e-9;
    bool close = true;

    CVec tempOutputVector(N);
    complex_fft(conjugateInputVec, tempOutputVector);
    tempOutputVector = complexConjugate(tempOutputVector);

    for(size_t i = 0; i < N; i++)
    {
        outputVector[i] = tempOutputVector[i];
    } 
}