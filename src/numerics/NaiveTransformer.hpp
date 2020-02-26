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
			NaiveTransformer(const int N) : mNumElements(N) { }

			void fft(const RealVec& inputVector, ComplexVec& outputVector) const override;
			void ifft(const ComplexVec& inputVector, RealVec& outputVector) const override;
		private:
			const int mNumElements;
		};
	}
}

#endif // !NAIVE_TRANSFORMER_HPP