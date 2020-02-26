#include "FFTWTransformer.hpp"
#include "fftw3.h"
#include <exception>

#ifdef _OPENMP
    #include <omp.h>
#endif

using namespace PCFT::numerics;

typedef FFTWTransformer::ComplexVec CVec;
typedef FFTWTransformer::RealVec RVec;

class FFTWTransformer::FFTWTransformerImpl
{
public:
    FFTWTransformerImpl(const int N) : mNumElements(N)
    {
        #ifdef _OPENMP
        // Initialisize multithreaded FFT automatically if OpenMP is
        // available
        FFT::init_multithread(omp_get_max_threads());
        #endif
        
        mIn = fftw_alloc_real(N);
        mOut = fftw_alloc_complex(N);

        const char* filename = "/wisdom_output_r";
        const int WISDOM_IMPORTED  = fftw_import_wisdom_from_filename(filename);     
        mPlanFw = fftw_plan_dft_r2c_1d(N, mIn, mOut, FFTW_MEASURE);
        mPlanBw = fftw_plan_dft_c2r_1d(N, mOut, mIn, FFTW_MEASURE);
        const int WISDOM_OUTPUTTED = fftw_export_wisdom_to_filename(filename);
    }

    ~FFTWTransformerImpl()
    {
        fftw_free(mIn);
        fftw_free(mOut);

        fftw_destroy_plan(mPlanFw);
        fftw_destroy_plan(mPlanBw);
    }

    void fft(const RVec& inputVec, CVec& outputVec)
    {
        for(int i = 0; i < mNumElements; i++)
        {
            mIn[i] = inputVec[i];
        }

        fftw_execute(mPlanFw);
            
        double sqrtN = sqrt(mNumElements);
        for(int i = 0; i < mNumElements; ++i)
        {
            outputVec[i].real(mOut[i][0]/sqrtN);
            outputVec[i].imag(mOut[i][1]/sqrtN);
        }
    }

    void ifft(const CVec& inputVec, RVec& outputVec)
    {
        for(int i = 0; i < mNumElements; ++i)
        {
            mOut[i][0] = inputVec[i].real();
            mOut[i][1] = inputVec[i].imag();
        }

        fftw_execute(mPlanBw);
            
        double sqrtN = sqrt(mNumElements);
        for(int i = 0; i < mNumElements; ++i)
        {
            outputVec[i] = mIn[i]/sqrtN;
        }
    }

    fftw_plan mPlanFw;
    fftw_plan mPlanBw;
    double* mIn;
    fftw_complex* mOut;
    const int mNumElements;
};

FFTWTransformer::FFTWTransformer(const int N) : mImpl(std::make_unique<FFTWTransformerImpl>(N)) { }

FFTWTransformer::~FFTWTransformer() = default;

int FFTWTransformer::getNumPoints() const
{
    return mImpl->mNumElements;
}

void FFTWTransformer::fft(const RVec& inputVector, CVec& outputVector) const
{
    mImpl->fft(inputVector, outputVector);
    // const size_t N = inputVector.size();
    // if (outputVector.size() != N)
    // {
    //     throw std::runtime_error("Invalid input vectors");
    // }

    // double *in;
    // fftw_complex *out;
    // fftw_plan p;

    // in = fftw_alloc_real(N);
    // out = fftw_alloc_complex(N);

    // for (size_t i = 0; i < N; i++)
    // {
    //     double& rIn = in[i];
    //     rIn = inputVector[i];
    // }

    // const char* filename = "/wisdom_output_r";
    // const int WISDOM_IMPORTED  = fftw_import_wisdom_from_filename(filename);
    // p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
    // const int WISDOM_OUTPUTTED = fftw_export_wisdom_to_filename(filename);

    // fftw_execute(p); /* repeat as needed */;

    // fftw_destroy_plan(p);

    // for (size_t i = 0; i < N; i++)
    // {
    //     const fftw_complex& rOut = out[i];
    //     std::complex<double> outputValue = { rOut[0] / sqrt(N), rOut[1] / sqrt(N) };
    //     outputVector[i] = outputValue;
    // }

    // fftw_free(in);
    // fftw_free(out);
}

void FFTWTransformer::ifft(const CVec& inputVector, RVec& outputVector) const
{
    mImpl->ifft(inputVector, outputVector);
    // const size_t N = inputVector.size();
    // if (outputVector.size() != N)
    // {
    //     throw std::runtime_error("Invalid input vectors");
    // }

    // fftw_complex *in;
    // double *out;
    // fftw_plan p;

    // in = fftw_alloc_complex(N);
    // out = fftw_alloc_real(N);

    // for (size_t i = 0; i < N; i++)
    // {
    //     fftw_complex& rIn = in[i];
    //     rIn[0] = inputVector[i].real();
    //     rIn[1] = inputVector[i].imag();
    // }

    // const char* filename = "wisdom_output_r";
    // const int WISDOM_IMPORTED = fftw_import_wisdom_from_filename(filename);
    // p = fftw_plan_dft_c2r_1d(N, in, out, FFTW_ESTIMATE);
    // const int WISDOM_OUTPUTTED = fftw_export_wisdom_to_filename(filename);

    // fftw_execute(p); /* repeat as needed */;

    // fftw_destroy_plan(p);

    // for (size_t i = 0; i < N; i++)
    // {
    //     const double& rOut = out[i];
    //     double outputValue = rOut / sqrt(N);
    //     outputVector[i] = outputValue;
    // }

    // fftw_free(in);
    // fftw_free(out);
}

// void FFTWTransformer::fft(const CVec& inputVector, CVec& outputVector) const
// {
//     const size_t N = inputVector.size();
//     if (outputVector.size() != N)
//     {
//         throw std::runtime_error("Invalid input vectors");
//     }

//     fftw_complex *in, *out;
//     fftw_plan p;

//     in = fftw_alloc_complex(N);
//     out = fftw_alloc_complex(N);

//     for (size_t i = 0; i < N; i++)
//     {
//         fftw_complex& rIn = in[i];
//         rIn[0] = inputVector[i].real();
//         rIn[1] = inputVector[i].imag();
//     }

//     const char* filename = "/wisdom_output";
//     const int WISDOM_IMPORTED  = fftw_import_wisdom_from_filename(filename);
//     p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
//     const int WISDOM_OUTPUTTED = fftw_export_wisdom_to_filename(filename);

//     fftw_execute(p); /* repeat as needed */;

//     fftw_destroy_plan(p);

//     for (size_t i = 0; i < N; i++)
//     {
//         const fftw_complex& rOut = out[i];
//         std::complex<double> outputValue = { rOut[0] / sqrt(N), rOut[1] / sqrt(N) };
//         outputVector[i] = outputValue;
//     }

//     fftw_free(in);
//     fftw_free(out);
// }


// void FFTWTransformer::ifft(const CVec& inputVector, CVec& outputVector) const
// {
//     const size_t N = inputVector.size();
//     if (outputVector.size() != N)
//     {
//         throw std::runtime_error("Invalid input vectors");
//     }

//     fftw_complex *in, *out;
//     fftw_plan p;

//     in = fftw_alloc_complex(N);
//     out = fftw_alloc_complex(N);

//     for (size_t i = 0; i < N; i++)
//     {
//         fftw_complex& rIn = in[i];
//         rIn[0] = inputVector[i].real();
//         rIn[1] = inputVector[i].imag();
//     }

//     const char* filename = "wisdom_output";
//     const int WISDOM_IMPORTED = fftw_import_wisdom_from_filename(filename);
//     p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
//     const int WISDOM_OUTPUTTED = fftw_export_wisdom_to_filename(filename);

//     fftw_execute(p); /* repeat as needed */;

//     fftw_destroy_plan(p);

//     for (size_t i = 0; i < N; i++)
//     {
//         const fftw_complex& rOut = out[i];
//         std::complex<double> outputValue = { rOut[0] / sqrt(N), rOut[1] / sqrt(N) };
//         outputVector[i] = outputValue;
//     }

//     fftw_free(in);
//     fftw_free(out);
// }