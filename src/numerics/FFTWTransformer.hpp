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
			FFTWTransformer(const int N);
			~FFTWTransformer();

			void fft(const RealVec& inputVector, ComplexVec& outputVector) const override;
			void ifft(const ComplexVec& inputVector, RealVec& outputVector) const override;

			int getNumPoints() const;

		private:
			class FFTWTransformerImpl;
			std::unique_ptr<FFTWTransformerImpl> mImpl;
		};
	}
}

#endif // !FFTW_TRANSFORMER_HPP