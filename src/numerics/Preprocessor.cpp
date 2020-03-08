#include "numerics/Preprocessor.hpp"

using namespace PCFT;
using namespace PCFT::numerics;

typedef IFourierTransformer::RealVec RVec;
typedef IFourierTransformer::ComplexVec CVec;

Preprocessor::Preprocessor(std::unique_ptr<IFourierTransformer> transformer,
				GFunction greensFunctionTransform,
				DomainParameters pPackage) :
                    mTransformer(std::move(transformer)),
					mGreensFunctionTransform(greensFunctionTransform),
					mN(pPackage.N), 
                    mDx(DomainParameters::getDx(pPackage.xMax, pPackage.xMin, pPackage.N)),
					mDtau(DomainParameters::getDTau(pPackage.T, pPackage.M)),
					mEpsilon1(pPackage.epsilon1),
					mEpsilon2(pPackage.epsilon2) { }

// The input vector here is the Green's function transform G, which we know in closed form.
RVec Preprocessor::execute() const
{
	CVec tempOutput(mN);
	CVec inputVector(mN);

	shiftedIfft(inputVector, tempOutput);
	shiftedFft(tempOutput, inputVector);

	RVec outputVector(mN);

	for(size_t i = 0; i < mN; i++)
	{
		outputVector[i] = inputVector[i].real();
	}

	return outputVector;
}

void Preprocessor::shiftedFft(CVec& inputVec, CVec& outputVec) const
{
	// The fourier transform mTransformer implements a normalised fourier transform (divided by sqrt(N)),
	// summing from j = 0 to N-1.
	// We would like the FFT with no scaling term sqrt(N), beginning at -N/2 and summing up to N/2-.

	for(auto& element : inputVec)
	{
		element.imag(-element.imag());
	}
	mTransformer->ifft(inputVec, outputVec);
	double sqrtN = sqrt(mN);
	for(auto& element : outputVec)
	{
		element = -element/sqrtN;
	}
}

void Preprocessor::shiftedIfft(const CVec& inputVec, CVec& outputVec) const
{
	// The fourier transform mTransformer implements a normalised fourier transform (divided by sqrt(N)),
	// summing from j = 0 to N-1.
	// We would like the FFT with no scaling term sqrt(N), beginning at -N/2 and summing up to N/2-.

	mTransformer->fft(inputVec, outputVec);
	double sqrtN = sqrt(mN);
	for(auto& element : outputVec)
	{
		element.real(-element.real()*sqrtN);
		element.imag(element.imag()*sqrtN);
	}
}

void Preprocessor::calculateH(IFourierTransformer::ComplexVec& HOut, int lambda) const
{
	// resize Hout
	int mNLambda = mN * lambda;
	HOut.resize(mNLambda);
}
