#include "IFourierTransformer.hpp"

namespace PCFT
{
    namespace numerics
    {
        class FTFactory
        {
            public:
                enum class TransformType
                {
                    FFT,
                    Naive
                };

            static std::unique_ptr<IFourierTransformer> instance();
            static std::unique_ptr<IFourierTransformer> instance(TransformType transformType);
        };
    }
}