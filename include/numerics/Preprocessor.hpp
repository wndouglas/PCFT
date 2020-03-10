#ifndef PREPROCESSOR_HPP
#define PREPROCESSOR_HPP

#include <vector>
#include "numerics/IFourierTransformer.hpp"
#include "DomainParameters.hpp"
#include "GFunction.hpp"

namespace PCFT
{
	namespace numerics
	{
		class Preprocessor
		{
		public:
			Preprocessor(GFunction greensFunctionTransform);

			IFourierTransformer::ComplexVec execute(
				const std::unique_ptr<IFourierTransformer>& mTransformer,
				const DomainParameters& pPackage);

		private:
			GFunction mGreensFunctionTransform;

			void calculateH(IFourierTransformer::ComplexVec& HOut, const int lambda, const int N, const double P, const double dx) const;
			void calculateLittleh(IFourierTransformer::ComplexVec& hOut, const std::unique_ptr<IFourierTransformer>& transformer, int lambda, const int N, const double P, const double dx) const;
			void calculateLittlegTilde(IFourierTransformer::ComplexVec& gOut, const std::unique_ptr<IFourierTransformer>& transformer, int lambda, const int N, const double P, const double dx) const;
			void calculateGTilde(IFourierTransformer::ComplexVec& GTilde, const IFourierTransformer::ComplexVec& gTilde, const std::unique_ptr<IFourierTransformer>& transformer, const double P) const;
			double test1(const IFourierTransformer::ComplexVec& gTilde, const int N, const int M, const double dx) const;
			double test2(const IFourierTransformer::ComplexVec& gTilde, const IFourierTransformer::ComplexVec& gTildePrev, const int N, const double dx) const;
		};
	}
}

#endif // !PREPROCESSOR_HPP