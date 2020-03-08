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
					mM(pPackage.M),
					mP(DomainParameters::getP(pPackage.xMax, pPackage.xMin)),
                    mDx(DomainParameters::getDx(pPackage.xMax, pPackage.xMin, pPackage.N)),
					mDtau(DomainParameters::getDTau(pPackage.T, pPackage.M)),
					mEpsilon1(pPackage.epsilon1),
					mEpsilon2(pPackage.epsilon2) { }

// The input vector here is the Green's function transform G, which we know in closed form.
CVec Preprocessor::execute() const
{
	CVec GTilde(mN);
	CVec gTildePrev(mN);
	CVec gTilde(mN);

	int lambda = 1;
	calculateLittlegTilde(gTildePrev, lambda);

	lambda = 2;
	double t1 = __DBL_MAX__;
	double t2 = t1;
	while(std::fabs(t1) >= mEpsilon1 ||
		std::fabs(t2) >= mEpsilon2)
	{
		calculateLittlegTilde(gTilde, lambda);
		calculateGTilde(GTilde, gTilde);

		t1 = test1(gTilde);
		t2 = test2(gTilde, gTildePrev);
		lambda *= 2;
		std::copy(gTilde.begin(), gTilde.end(), gTildePrev.begin());
	}

	return GTilde;
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
		HOut[i] = sincSquared(xIn) * mGreensFunctionTransform(omegaK);
	}
}

void Preprocessor::calculateLittleh(IFourierTransformer::ComplexVec& hOut, int lambda) const
{
	calculateH(hOut, lambda);
	mTransformer->shiftedIfft(hOut, hOut);

	int mNLambda = mN * lambda;
	for (int i = 0; i < mNLambda; i++)
	{
		hOut[i] /= mP;
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

void Preprocessor::calculateGTilde(IFourierTransformer::ComplexVec& GTilde, const IFourierTransformer::ComplexVec& gTilde) const
{
	GTilde.resize(gTilde.size());
	mTransformer->shiftedFft(gTilde, GTilde);


	for (int i = 0; i < GTilde.size(); i++)
	{
		GTilde[i] *= mP;
	}
}

double Preprocessor::test1(const IFourierTransformer::ComplexVec& gTilde) const
{
	double output = 0.0;
	for(int i = 0; i < mN; i++)
	{
		output += std::min(gTilde[i].real(), 0.0);
	}
	output *= -mM*mDx;
	return output;
}

double Preprocessor::test2(const IFourierTransformer::ComplexVec& gTilde, const IFourierTransformer::ComplexVec& gTildePrev) const
{
	double runningMax = 0.0;
	for(int i = 0; i < mN; i++)
	{
		double currMax = std::fabs(gTilde[i].real() - gTildePrev[i].real());
		if(currMax >= runningMax)
		{
			runningMax = currMax;
		}
	}
	return mDx*runningMax;
}
