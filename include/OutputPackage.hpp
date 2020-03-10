#ifndef OUTPUT_PACKAGE_HPP
#define OUTPUT_PACKAGE_HPP

#include <vector>

namespace PCFT
{
    struct OutputPackage
    {
        std::vector<double> mTimes;
        std::vector<double> mSpace;
        std::vector<std::vector<double> > mValueFunction;
        std::vector<std::vector<double> > mPolicy;
    };
}
#endif // !OUTPUT_PACKAGE_HPP