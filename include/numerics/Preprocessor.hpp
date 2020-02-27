#ifndef PREPROCESSOR_HPP
#define PREPROCESSOR_HPP

#include <vector>
#include "IFourierTransformer.hpp"

namespace PCFT
{
	namespace numerics
	{
		class Preprocessor
		{
		public:
			Preprocessor(
				std::unique_ptr<IFourierTransformer> transformer,
				const int N,
				const double dx, 
				const double dTau, 
				const double epsilon1,
				const double epsilon2);

			void execute(const std::vector<double>& inputVector, std::vector<double>& outputVector);

		private:
			std::unique_ptr<IFourierTransformer> mTransformer;
			const int mN;
			const double mDx;
			const double mDtau;
			const double mEpsilon1;
			const double mEpsilon2;
		};
	}
}

#endif // !PREPROCESSOR_HPP