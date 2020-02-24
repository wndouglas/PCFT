#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>

namespace PCFT
{
    namespace tools
    {
        class Timer
        {
        public:
            typedef double milliseconds;
            void start();
            void stop();
            milliseconds duration();

            Timer();

        private:
            std::chrono::steady_clock::time_point mStartTime;
            std::chrono::steady_clock::time_point mFinishTime;
            bool mClockRunning;
            bool mHasClockStopped;
        };
    }
}

#endif // !TIMER_HPP