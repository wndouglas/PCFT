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
			Preprocessor(
				std::unique_ptr<IFourierTransformer> transformer,
				GFunction greensFunctionTransform,
				DomainParameters pPackage);

			IFourierTransformer::ComplexVec execute() const;

		private:
			std::unique_ptr<IFourierTransformer> mTransformer;
			GFunction mGreensFunctionTransform;
			const int mN;
			const int mM;
			const double mP;
			const double mDx;
			const double mDtau;
			const double mEpsilon1;
			const double mEpsilon2;

			void calculateH(IFourierTransformer::ComplexVec& HOut, int lambda) const;
			void calculateLittleh(IFourierTransformer::ComplexVec& hOut, int lambda) const;
			void calculateLittlegTilde(IFourierTransformer::ComplexVec& gOut, int lambda) const;
			void calculateGTilde(IFourierTransformer::ComplexVec& GTilde, const IFourierTransformer::ComplexVec& gTilde) const;
			double test1(const IFourierTransformer::ComplexVec& gTilde) const;
			double test2(const IFourierTransformer::ComplexVec& gTilde, const IFourierTransformer::ComplexVec& gTildePrev) const;
		};
	}
}

#endif // !PREPROCESSOR_HPP