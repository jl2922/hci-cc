#ifndef STATUS_H_
#define STATUS_H_

#include <chrono>
#include <cstdio>
#include <string>

class Status {
 public:
  static char* time() {
    using namespace std::chrono;
    static Status instance;
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
  Status(const Status&) = delete;
  void operator=(const Status&) = delete;

 private:
  Status() { last = begin = std::chrono::high_resolution_clock::now(); }
  std::chrono::high_resolution_clock::time_point begin;
  std::chrono::high_resolution_clock::time_point last;
};

#endif