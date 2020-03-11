#ifndef PCFT_EXECUTOR_HPP
#define PCFT_EXECUTOR_HPP

#include "OutputPage.hpp"
#include "DomainParameters.hpp"
#include "numerics/IFourierTransformer.hpp"
#include "numerics/Preprocessor.hpp"

namespace PCFT
{
    class PCFTExecutor
    {
    public:
        PCFTExecutor(
            std::unique_ptr<numerics::IFourierTransformer> transformer,
            std::unique_ptr<numerics::Preprocessor> preprocessor,
            DomainParameters pPackage) :
                mTransformer(std::move(transformer)),
                mPreprocessor(std::move(preprocessor)),
                mParams(pPackage) {}

        void execute(OutputPage& outputPage);

    private:
        std::unique_ptr<numerics::IFourierTransformer> mTransformer;
        std::unique_ptr<numerics::Preprocessor> mPreprocessor;
        DomainParameters mParams;
    };
}


#endif // !PCFT_EXECUTOR_HPP