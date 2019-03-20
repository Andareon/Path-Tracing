#ifndef RAY_TRACING_TIMER_H
#define RAY_TRACING_TIMER_H

#include <chrono>
#include <iostream>
#include <map>

class TimerStorage {
    struct TimerSummary {
        long long minValue = std::numeric_limits<long long>::max();
        long long maxValue = 0;
        long long sum = 0;
        long long count = 0;
    };
private:
    TimerStorage() = default;

public:
    TimerStorage(const TimerStorage&) = delete;
    TimerStorage& operator=(TimerStorage&) = delete;
    std::map<std::string, TimerSummary> timers;
    static TimerStorage& get() {
        static TimerStorage instance;
        return instance;
    }

    void AddTimerValue(const std::string &timer_name, long long time) {
        timers[timer_name].sum += time;
        timers[timer_name].count++;
        timers[timer_name].maxValue = std::max(timers[timer_name].maxValue, time);
        timers[timer_name].minValue = std::min(timers[timer_name].minValue, time);
    }

    void PrintTimers() {
        for (auto &timer : timers) {
            std::cout << timer.first << " " << timer.second.minValue << " " << (float) timer.second.sum / timer.second.count <<
                         " " << timer.second.maxValue << std::endl;
        }
    }
};

template<typename T>
T get_current_time() {
    return std::chrono::duration_cast<T>(std::chrono::system_clock::now().time_since_epoch());
}

class Timer {
private:
    std::string timer_name_;
    std::chrono::microseconds start_time_;

public:
    explicit Timer(std::string name) :
        timer_name_(std::move(name)),
        start_time_(get_current_time<std::chrono::microseconds>())
    {}

    ~Timer() {
        std::chrono::microseconds end_time = get_current_time<std::chrono::microseconds>();
        TimerStorage::get().AddTimerValue(timer_name_, (end_time - start_time_).count());
    }
};

#endif //RAY_TRACING_TIMER_H
