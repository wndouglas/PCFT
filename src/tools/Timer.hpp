#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>

namespace PCFT
{
    namespace tools
    {
        template<typename type = float, typename period = std::milli>
        class Timer
        {
        public:
            using clock = std::chrono::high_resolution_clock;
            using duration = std::chrono::duration<type, period>;
            using time_point = std::chrono::time_point<clock, duration>;

            Timer() : time_(clock::now()) { }
            Timer(const Timer& that) = default;
            Timer(Timer&& temp) = default;
            ~Timer() = default;
            Timer& operator=(const Timer& that) = default;
            Timer& operator=(Timer&& temp) = default;

            duration tick()
            {
                time_point time = clock::now();
                duration delta = time - time_;
                time_ = time;

                return delta;
            }

        private:
            time_point time_;
        };
    }
}

#endif // !TIMER_HPP