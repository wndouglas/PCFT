#include "numerics/Preprocessor.hpp"

using namespace PCFT;
using namespace PCFT::numerics;

typedef IFourierTransformer::RealVec RVec;
typedef IFourierTransformer::ComplexVec CVec;

namespace
{
   const double taylor_n_bound = 1e-5;
   double sinc(double x)
   {
    	if (fabs(x) >= taylor_n_bound)
		{
			return sin(x) / x;
		}
     	else
		{
			double xSquared = x*x;
			return 1 - xSquared/6 + xSquared*xSquared/120;
		}
   }

}

Preprocessor::Preprocessor(std::unique_ptr<IFourierTransformer> transformer,
				DomainParameters pPackage) :
                    mTransformer(std::move(transformer)),
					mN(pPackage.N),
					mP(getP(pPackage.xMax, pPackage.xMin)),
					mXMin(pPackage.xMin),
					mXMax(pPackage.xMax),
                    mDx(getDx(pPackage.xMax, pPackage.xMin, pPackage.N)),
					mDtau(getDTau(pPackage.T, pPackage.M)),
					mEpsilon1(pPackage.epsilon1),
					mEpsilon2(pPackage.epsilon2),
					mLambda(1) { }

// The input vector here is the Green's function transform G, which we know in closed form.
RVec Preprocessor::execute(const GFunction& greensFunctionTransform) const
{
	CVec inputVec(mN);
	calculateH(inputVec, greensFunctionTransform);

	CVec tempOutput(mN);
	CVec tempInput(mN);
	shiftedIfft(inputVec, tempOutput);
	shiftedFft(tempOutput, tempInput);

	double l2Err = 0.0;
	for(int i = 0; i < mN; i++)
	{
		l2Err += pow(tempInput[i].real() - inputVec[i].real(), 2) + pow(tempInput[i].imag() - inputVec[i].imag(), 2);
	}
	return {sqrt(l2Err/mN)};
}

void Preprocessor::shiftedFft(CVec& inputVec, CVec& outputVec) const
{
	// The fourier transform mTransformer implements a normalised fourier transform (divided by sqrt(N)),
	// summing from j = 0 to N-1.
	// We would like the FFT with scaling term 1/N, beginning at -N/2 and summing up to N/2-1.
	mTransformer->fft(inputVec, outputVec);
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
	// We would like the FFT with no scaling term, beginning at -N/2 and summing up to N/2-1.

	mTransformer->ifft(inputVec, outputVec);
	double sqrtN = sqrt(mN);
	for(auto& element : outputVec)
	{
		element.real(-element.real()*sqrtN);
		element.imag(-element.imag()*sqrtN);
	}
}

void Preprocessor::calculateH(CVec& vector, const GFunction& greensFunctionTransform) const
{
	if(vector.size() != mN)
	{
		vector.resize(mN);
	}

	for(int i = 0; i < mN; i++)
	{
		int k = i - mN/2;
		double omega_k = k/mP;
		double sincOmegaSquared = sinc(PI*omega_k*mDx)*sinc(PI*omega_k*mDx);
		std::complex<double> GF = greensFunctionTransform(omega_k);
		vector[i] = sincOmegaSquared*GF;
	}
}
