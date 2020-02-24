#include "FFTWTransformer.hpp"
#include "fftw3.h"
#include <exception>

using namespace PCFT::numerics;

typedef FFTWTransformer::ComplexVec CVec;

void FFTWTransformer::fft(const CVec& inputVector, CVec& outputVector) const
{
    const size_t N = inputVector.size();
    if (outputVector.size() != N)
    {
        throw std::runtime_error("Invalid input vectors");
    }

    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);

    for (size_t i = 0; i < N; i++)
    {
        fftw_complex& rIn = in[i];
        rIn[0] = inputVector[i].real();
        rIn[1] = inputVector[i].imag();
    }

    const char* filename = "wisdom_output";
    const int WISDOM_IMPORTED  = fftw_import_wisdom_from_filename(filename);
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_MEASURE);
    const int WISDOM_OUTPUTTED = fftw_export_wisdom_to_filename(filename);

    fftw_execute(p); /* repeat as needed */;

    fftw_destroy_plan(p);

    for (size_t i = 0; i < N; i++)
    {
        const fftw_complex& rOut = out[i];
        std::complex<double> outputValue = { rOut[0] / sqrt(N), rOut[1] / sqrt(N) };
        outputVector[i] = outputValue;
    }

    fftw_free(in);
    fftw_free(out);
}

void FFTWTransformer::ifft(const CVec& inputVector, CVec& outputVector) const
{
    const size_t N = inputVector.size();
    if (outputVector.size() != N)
    {
        throw std::runtime_error("Invalid input vectors");
    }

    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);

    for (size_t i = 0; i < N; i++)
    {
        fftw_complex& rIn = in[i];
        rIn[0] = inputVector[i].real();
        rIn[1] = inputVector[i].imag();
    }

    const char* filename = "wisdom_output";
    const int WISDOM_IMPORTED = fftw_import_wisdom_from_filename(filename);
    p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_MEASURE);
    const int WISDOM_OUTPUTTED = fftw_export_wisdom_to_filename(filename);

    fftw_execute(p); /* repeat as needed */;

    fftw_destroy_plan(p);

    for (size_t i = 0; i < N; i++)
    {
        const fftw_complex& rOut = out[i];
        std::complex<double> outputValue = { rOut[0] / sqrt(N), rOut[1] / sqrt(N) };
        outputVector[i] = outputValue;
    }

    fftw_free(in);
    fftw_free(out);
}
