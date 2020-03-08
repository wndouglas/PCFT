#include "numerics/Preprocessor.hpp"
#include <iostream>

using namespace PCFT;
using namespace PCFT::numerics;

typedef IFourierTransformer::RealVec RVec;
typedef IFourierTransformer::ComplexVec CVec;

namespace
{
	double sincSquared(double x)
	{
		x = std::fabs(x);
		if (x < 1e-3)
		{
			double xSquared = x * x;
			return 1 - xSquared / 3 + 2 * xSquared * xSquared / 45;
		}
		else
		{
			return std::sin(x) * std::sin(x) / x / x;
		}
	}
}


Preprocessor::Preprocessor(std::unique_ptr<IFourierTransformer> transformer,
				GFunction greensFunctionTransform,
				DomainParameters pPackage) :
                    mTransformer(std::move(transformer)),
					mGreensFunctionTransform(greensFunctionTransform),
					mN(pPackage.N), 
					mP(DomainParameters::getP(pPackage.xMax, pPackage.xMin)),
                    mDx(DomainParameters::getDx(pPackage.xMax, pPackage.xMin, pPackage.N)),
					mDtau(DomainParameters::getDTau(pPackage.T, pPackage.M)),
					mEpsilon1(pPackage.epsilon1),
					mEpsilon2(pPackage.epsilon2) { }

// The input vector here is the Green's function transform G, which we know in closed form.
RVec Preprocessor::execute() const
{
	RVec outputVector(mN);
	CVec gTilde(mN);
	calculateLittlegTilde(gTilde, 1);
	return outputVector;
}

void Preprocessor::shiftedFft(CVec& inputVec, CVec& outputVec) const
{
	// The fourier transform mTransformer implements a normalised fourier transform (divided by sqrt(N)),
	// summing from j = 0 to N-1.
	// We would like the FFT with no scaling term sqrt(N), beginning at -N/2 and summing up to N/2-.
	const size_t N = inputVec.size();
	CVec tempInput(N);
	for (int r = 0; r < N; r++)
	{
		double flipper = 1 - (r % 2) * 2;
		tempInput[r] = flipper * inputVec[r];
	}
	mTransformer->fft(inputVec, outputVec);
	double sqrtN = sqrt(N);
	for (int l = 0; l < N; l++)
	{
		double flipper = 1 - (l % 2) * 2;
		outputVec[l] *= flipper/sqrt(N);
	}
}

void Preprocessor::shiftedIfft(const CVec& inputVec, CVec& outputVec) const
{
	// The fourier transform mTransformer implements a normalised fourier transform (divided by sqrt(N)),
	// summing from j = 0 to N-1.
	// We would like the FFT with no scaling term sqrt(N), beginning at -N/2 and summing up to N/2-.
	const size_t N = inputVec.size();
	CVec tempInput(N);
	for (int r = 0; r < N; r++)
	{
		double flipper = 1 - (r % 2) * 2;
		tempInput[r] = flipper * inputVec[r];
	}
	mTransformer->ifft(tempInput, outputVec);
	double sqrtN = sqrt(N);

	for(int l = 0; l < N; l++)
	{
		double flipper = 1 - (l % 2) * 2;
		outputVec[l] *= flipper*sqrt(N);
	}
}

void Preprocessor::calculateH(IFourierTransformer::ComplexVec& HOut, int lambda) const
{
	// resize Hout
	int mNLambda = mN * lambda;
	// Calculate H
	for (int i = 0; i < mNLambda; i++)
	{
		int k = i - mNLambda / 2;
		double omegaK = k / mP;
		double xIn = PI * omegaK * mDx;

		//HOut[i] = sincSquared(xIn) * mGreensFunctionTransform(omegaK);
		HOut[i] = mGreensFunctionTransform(omegaK);
	}
}

void Preprocessor::calculateLittleh(IFourierTransformer::ComplexVec& hOut, int lambda) const
{
	calculateH(hOut, lambda);
	CVec tempOut(hOut.size());
	shiftedIfft(hOut, tempOut);

	int mNLambda = mN * lambda;
	for (int i = 0; i < mNLambda; i++)
	{
		hOut[i] = tempOut[i] / mP;
	}
}

void Preprocessor::calculateLittlegTilde(IFourierTransformer::ComplexVec& gOut, int lambda) const
{
	int mNLambda = mN * lambda;
	if (mN % 2 != 0 || mNLambda % 2 != 0)
	{
		throw std::runtime_error("Invalid input argument");
	}

	CVec hOut(mNLambda);
	calculateLittleh(hOut, lambda);

	for (int i = 0; i < mN; i++)
	{
		int j = i - mN / 2;
		int l = lambda * j + mNLambda/2;

		gOut[i] = hOut[l];
	}
}

