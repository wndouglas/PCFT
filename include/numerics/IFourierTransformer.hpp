#ifndef IFOURIER_TRANSFORMER_HPP
#define IFOURIER_TRANSFORMER_HPP

#include <vector>
#include <complex>

namespace PCFT
{
namespace numerics
{
	// The following class provides the abstract interface for defining a Fourier transform method
	// which takes a complex input and returns a complex output in forward, 'fft' mode, and takes a complex
	// input and returns a complex output in 'ifft' mode.
	class IFourierTransformer
	{
	public:
		typedef std::vector<std::complex<double> > ComplexVec;
		typedef std::vector<double> RealVec;

		virtual void fft(const ComplexVec& inputVector, ComplexVec& outputVector) const = 0;
		virtual void ifft(const ComplexVec& inputVector, ComplexVec& outputVector) const = 0;

		virtual void shiftedFft(const ComplexVec& inputVector, ComplexVec& outputVector) const = 0;
		virtual void shiftedIfft(const ComplexVec& inputVector, ComplexVec& outputVector) const = 0;


		virtual ~IFourierTransformer() { };
	};
}
}

#endif // !IFOURIER_TRANSFORMER_HPP
