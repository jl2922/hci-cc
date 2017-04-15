#ifndef TIMER_H_
#define TIMER_H_

#include <chrono>
#include <cstdio>
#include <string>
#include <unordered_map>

class Timer {
 public:
  static char* str() {
    using namespace std::chrono;
    static Timer instance;
    auto now = high_resolution_clock::now();
    double time_diff =
        duration_cast<duration<double>>(now - instance.last).count();
    double time_total =
        duration_cast<duration<double>>(now - instance.begin).count();
    instance.last = now;
    static char buf[100];
    snprintf(buf, sizeof(buf), "[time: %.3f/%.3f]", time_total, time_diff);
    return buf;
  }
  Timer(const Timer&) = delete;
  void operator=(const Timer&) = delete;

 private:
  Timer() { last = begin = std::chrono::high_resolution_clock::now(); }
  std::chrono::high_resolution_clock::time_point begin;
  std::chrono::high_resolution_clock::time_point last;
};

#endif