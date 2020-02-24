#ifndef FFTW_TRANSFORMER_HPP
#define FFTW_TRANSFORMER_HPP

#include "numerics/IFourierTransformer.hpp"

namespace PCFT
{
	namespace numerics
	{
		class FFTWTransformer : public IFourierTransformer
		{
		public:
			void fft(const ComplexVec& inputVector, ComplexVec& outputVector) const override;
			void ifft(const ComplexVec& inputVector, ComplexVec& outputVector) const override;
		};
	}
}

#endif // !FFTW_TRANSFORMER_HPP