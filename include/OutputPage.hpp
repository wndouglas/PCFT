#ifndef OUTPUT_PAGE_HPP
#define OUTPUT_PAGE_HPP

#include <vector>

namespace PCFT
{
    struct OutputPage
    {
        std::vector<double> times;
        std::vector<double> state;
        std::vector<std::vector<double> > valueFunction;
        std::vector<std::vector<double> > policy;
    };
}
#endif // !OUTPUT_PAGE_HPP