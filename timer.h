#ifndef RAY_TRACING_TIMER_H
#define RAY_TRACING_TIMER_H

#include <chrono>
#include <iostream>
#include <map>

class TimerStorage {
    struct TimerSummary {
        int minValue;
        int maxValue;
        int sum;
        int count;
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

    void AddTimerValue(const std::string &timer_name, int time) {
        timers[timer_name].sum += time;
        timers[timer_name].count++;
        timers[timer_name].maxValue = std::max(timers[timer_name].maxValue, time);
        timers[timer_name].minValue = std::min(timers[timer_name].minValue, time);
    }

    void PrintTimers() {
        for (auto it = timers.begin(); it != timers.end(); it++) {
            std::cout << it->first << " " << it->second.maxValue << " " << (float)it->second.sum / it->second.count <<
                         " " << it->second.maxValue << std::endl;
        }
    }
};

class Timer {
private:
    std::string timer_name_;
    std::chrono::milliseconds start_time_;

public:
    Timer(std::string name)
              : timer_name_(std::move(name)) {
        start_time_ = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::system_clock::now().time_since_epoch());
    }

    ~Timer() {
        std::chrono::milliseconds end_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::system_clock::now().time_since_epoch());
        TimerStorage::get().AddTimerValue(timer_name_, (end_time - start_time_).count());
    }
};

#endif //RAY_TRACING_TIMER_H
