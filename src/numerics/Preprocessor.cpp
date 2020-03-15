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


Preprocessor::Preprocessor(GFunction greensFunctionTransform) : mGreensFunctionTransform(greensFunctionTransform) { }

// The input vector here is the Green's function transform G, which we know in closed form.
CVec Preprocessor::execute(const IFourierTransformer* transformer, const DomainParameters& pPackage)
{
	const double dtau = DomainParameters::getDTau(pPackage.T, pPackage.M);
	const double mEpsilon1 = pPackage.epsilon1;
	const double mEpsilon2 = pPackage.epsilon2;

    double xMin = pPackage.xMin;
    double xMax = pPackage.xMax;
    double P = DomainParameters::getP(xMax, xMin);
    int N = pPackage.N;
	int M = pPackage.M;
	double dx = pPackage.getDx(xMax, xMin, N);

	CVec GTilde(N);
	CVec gTildePrev(N);
	CVec gTilde(N);

	int lambda = 1;
	calculateLittlegTilde(gTildePrev, transformer, lambda, N, P, dx);

	lambda = 2;
	double t1 = __DBL_MAX__;
	double t2 = t1;
	while(std::fabs(t1) >= mEpsilon1 ||
		std::fabs(t2) >= mEpsilon2)
	{
		calculateLittlegTilde(gTilde, transformer, lambda, N, P, dx);
		calculateGTilde(GTilde, gTilde, transformer, P);

		t1 = test1(gTilde, N, M, dx);
		t2 = test2(gTilde, gTildePrev, N, dx);
		lambda *= 2;
		std::copy(gTilde.begin(), gTilde.end(), gTildePrev.begin());
	}

	return GTilde;
}

void Preprocessor::calculateH(IFourierTransformer::ComplexVec& HOut,
							  const int lambda, const int N, const double P, const double dx) const
{
	// resize Hout
	int mNLambda = N * lambda;
	// Calculate H
	for (int i = 0; i < mNLambda; i++)
	{
		int k = i - mNLambda / 2;
		double omegaK = k / P;
		double xIn = PI * omegaK * dx;
		HOut[i] = sincSquared(xIn) * mGreensFunctionTransform(omegaK);
	}
}

void Preprocessor::calculateLittleh(IFourierTransformer::ComplexVec& hOut,
								    const IFourierTransformer* transformer,
									int lambda, const int N, const double P, const double dx) const
{
	calculateH(hOut, lambda, N, P, dx);
	transformer->shiftedIfft(hOut, hOut);

	int NLambda = N * lambda;
	for (int i = 0; i < NLambda; i++)
	{
		hOut[i] /= P;
	}
}

void Preprocessor::calculateLittlegTilde(IFourierTransformer::ComplexVec& gOut,
										 const IFourierTransformer* transformer,
										 int lambda, const int N, const double P, const double dx) const
{
	int NLambda = N * lambda;
	if (N % 2 != 0 || NLambda % 2 != 0)
	{
		throw std::runtime_error("Invalid input argument");
	}

	CVec hOut(NLambda);
	calculateLittleh(hOut, transformer, lambda, N, P, dx);

	for (int i = 0; i < N; i++)
	{
		int j = i - N / 2;
		int l = lambda * j + NLambda/2;

		gOut[i] = hOut[l];
	}
}

void Preprocessor::calculateGTilde(IFourierTransformer::ComplexVec& GTilde,
								   const IFourierTransformer::ComplexVec& gTilde,
								   const IFourierTransformer* transformer,
								   const double P) const
{
	GTilde.resize(gTilde.size());
	transformer->shiftedFft(gTilde, GTilde);


	for (int i = 0; i < GTilde.size(); i++)
	{
		GTilde[i] *= P;
	}
}

double Preprocessor::test1(const IFourierTransformer::ComplexVec& gTilde, 
						   const int N, const int M, const double dx) const
{
	double output = 0.0;
	for(int i = 0; i < N; i++)
	{
		output += std::min(gTilde[i].real(), 0.0);
	}
	output *= -M*dx;
	return output;
}

double Preprocessor::test2(const IFourierTransformer::ComplexVec& gTilde,
						   const IFourierTransformer::ComplexVec& gTildePrev, const int N, const double dx) const
{
	double runningMax = 0.0;
	for(int i = 0; i < N; i++)
	{
		double currMax = std::fabs(gTilde[i].real() - gTildePrev[i].real());
		if(currMax >= runningMax)
		{
			runningMax = currMax;
		}
	}
	return dx*runningMax;
}
