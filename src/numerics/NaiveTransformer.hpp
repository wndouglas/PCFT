#ifndef NAIVE_TRANSFORMER_HPP
#define NAIVE_TRANSFORMER_HPP

#include "numerics/IFourierTransformer.hpp"

namespace PCFT
{
	namespace numerics
	{
		class NaiveTransformer : public IFourierTransformer
		{
		public:
			void fft(const ComplexVec& inputVector, ComplexVec& outputVector) const override;
			void ifft(const ComplexVec& inputVector, ComplexVec& outputVector) const override;
		};
	}
}

#endif // !NAIVE_TRANSFORMER_HPP