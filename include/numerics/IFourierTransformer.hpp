#ifndef IFOURIER_TRANSFORMER_HPP
#define IFOURIER_TRANSFORMER_HPP

#include <vector>
#include <complex>

namespace PCFT
{
namespace numerics
{
	class IFourierTransformer
	{
	public:
		typedef std::vector<std::complex<double> > ComplexVec;

		virtual void fft(const ComplexVec& inputVector, ComplexVec& outputVector) const = 0;
		virtual void ifft(const ComplexVec& inputVector, ComplexVec& outputVector) const = 0;

		virtual ~IFourierTransformer() { };
	};
}
}

#endif // !IFOURIER_TRANSFORMER_HPP
