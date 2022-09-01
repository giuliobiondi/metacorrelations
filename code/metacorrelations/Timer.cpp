#include "Timer.h"
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

namespace Time
{
    //ticks
	extern unsigned long long get_ticks_time()
	{
		unsigned long long val;
#ifdef _WIN32
		QueryPerformanceCounter((LARGE_INTEGER *)&val);
#else 
		timeval time_val;
		gettimeofday(&time_val, nullptr);
		val = (unsigned long long)time_val.tv_sec * (1000 * 1000) + (unsigned long long)time_val.tv_usec;
#endif
		return val;
	}
	//time
	extern double get_ms_time()
	{
		static double coe = 1.0 / 1000.0;
#ifdef _WIN32
		static unsigned long long freq = 0;
		if (freq == 0)
		{
			QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
			coe = 1000.0 / freq;
		}
#endif
		return get_ticks_time() * coe;
	}

    Timer::Timer(bool _start)
    {
        if(_start)
            start();
    }

    void Timer::start()
    {
        start_time = get_ms_time();
    }

    double Timer::stop()
    {
        end_time = get_ms_time();
        return result();
    }

    double Timer::result()
    {
        return end_time - start_time;
    }

    Timer::operator double()
    {
        return result();
    }
}