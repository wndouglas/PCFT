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
    FFTWTransformerImpl(size_t N) : mNumElements(N), mNumElementsChanged(false)
    {
        #ifdef _OPENMP
        // Initialisize multithreaded FFT automatically if OpenMP is
        // available
        FFT::init_multithread(omp_get_max_threads());
        #endif
        
        mIn = fftw_alloc_real(N);
        mOut = fftw_alloc_complex(N/2 + 1);
        
        const char* filename = "wisdom_output_r";
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
        if(mNumElementsChanged)
        {
            fftw_free(mIn);
            fftw_free(mOut);

            mIn = fftw_alloc_real(mNumElements);
            mOut = fftw_alloc_complex(mNumElements/2 + 1);

            mPlanFw = fftw_plan_dft_r2c_1d(mNumElements, mIn, mOut, FFTW_ESTIMATE);
            mPlanBw = fftw_plan_dft_c2r_1d(mNumElements, mOut, mIn, FFTW_ESTIMATE);
            mNumElementsChanged = false;
        }

        for(int i = 0; i < mNumElements; i++)
        {
            mIn[i] = inputVec[i];
        }

        fftw_execute(mPlanFw);
            
        double sqrtN = sqrt(mNumElements);
        for(int i = 0; i < mNumElements/2 + 1; ++i)
        {
            outputVec[i].real(mOut[i][0]/sqrtN);
            outputVec[i].imag(mOut[i][1]/sqrtN);
        }
        for(int i = mNumElements/2 + 1; i < mNumElements; i++)
        {
            int j = mNumElements - i;
            outputVec[i].real(outputVec[j].real());
            outputVec[i].imag(-outputVec[j].imag());
        }
    }

    void ifft(const CVec& inputVec, RVec& outputVec)
    {
        if(mNumElementsChanged)
        {
            fftw_free(mIn);
            fftw_free(mOut);

            mIn = fftw_alloc_real(mNumElements);
            mOut = fftw_alloc_complex(mNumElements/2 + 1);

            mPlanFw = fftw_plan_dft_r2c_1d(mNumElements, mIn, mOut, FFTW_ESTIMATE);
            mPlanBw = fftw_plan_dft_c2r_1d(mNumElements, mOut, mIn, FFTW_ESTIMATE);
            mNumElementsChanged = false;
        }

        for(int i = 0; i < mNumElements/2 + 1; ++i)
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
    size_t mNumElements;
    bool mNumElementsChanged;
};

FFTWTransformer::FFTWTransformer(const size_t N) : mImpl(std::make_unique<FFTWTransformerImpl>(N)) { }

FFTWTransformer::~FFTWTransformer() = default;

size_t FFTWTransformer::getNumPoints() const
{
    return mImpl->mNumElements;
}

void FFTWTransformer::setNumPoints(size_t N)
{
    mImpl->mNumElements = N;
    mImpl->mNumElementsChanged = true;
}

void FFTWTransformer::fft(const RVec& inputVector, CVec& outputVector) const
{
    if(inputVector.size() != mImpl->mNumElements)
    {
        mImpl->mNumElements = inputVector.size();
        mImpl->mNumElementsChanged = true;
    }

    mImpl->fft(inputVector, outputVector);
}

void FFTWTransformer::ifft(const CVec& inputVector, RVec& outputVector) const
{
    if(inputVector.size() != mImpl->mNumElements)
    {
        mImpl->mNumElements = inputVector.size();
        mImpl->mNumElementsChanged = true;
    }

    mImpl->ifft(inputVector, outputVector);
}