#include "numerics/FTFactory.hpp"
#include "NaiveTransformer.hpp"

using namespace PCFT::numerics;
using namespace std;

#ifdef USING_FFTW3
#include "FFTWTransformer.hpp"

unique_ptr<IFourierTransformer> FTFactory::instance()
{
    return instance(TransformType::FFT);
}

unique_ptr<IFourierTransformer> FTFactory::instance(FTFactory::TransformType transformType)
{
    unique_ptr<IFourierTransformer> transformer;
    switch(transformType)
    {
        case TransformType::FFT:
            transformer = make_unique<FFTWTransformer>();
            break;
        default:
            transformer = make_unique<NaiveTransformer>();
            break;
    }
    return transformer;
}
#else
unique_ptr<IFourierTransformer> FTFactory::instance()
{
    return instance(TransformType::Naive);
}

unique_ptr<IFourierTransformer> FTFactory::instance(FTFactory::TransformType transformType)
{
    unique_ptr<IFourierTransformer> transformer;
    switch(transformType)
    {
        case TransformType::FFT:
            throw runtime_error("FFT type transform has no implementation present.");
            break;
        default:
            transformer = make_unique<NaiveTransformer>();
            break;
    }
    return transformer;
}
#endif