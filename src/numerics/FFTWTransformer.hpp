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
			FFTWTransformer(const size_t N);
			~FFTWTransformer();

			void fft(const ComplexVec& inputVector, ComplexVec& outputVector) const override;
			void ifft(const ComplexVec& inputVector, ComplexVec& outputVector) const override;

			size_t getNumPoints() const;
			void setNumPoints(size_t N);

		private:
			class FFTWTransformerImpl;
			std::unique_ptr<FFTWTransformerImpl> mImpl;
		};
	}
}

#endif // !FFTW_TRANSFORMER_HPP