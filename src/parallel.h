#ifndef HCI_PARALLEL_H_
#define HCI_PARALLEL_H_

#include <cstdarg>
#include <cstddef>
#include <list>
#include "mpi.h"

namespace hci {

class Parallel {
  public:
    static Parallel& getInstance() {
      static Parallel instance;
      return instance;
    }
    void init(int argc, char ** argv) {
      MPI_Init(&argc, &argv);
      MPI_Comm_size(MPI_COMM_WORLD, &n_tasks);
      MPI_Comm_rank(MPI_COMM_WORLD,&task_id);
      int len;
      MPI_Get_processor_name(hostname, &len);
      printf ("Task %d: Ready running on %s.\n", task_id, hostname);
      checkpoint("All %d task(s) ready.\n", n_tasks);
    }
    int get_n_tasks() { return n_tasks; }
    int get_task_id() { return task_id; }
    void checkpoint(const char* format, ...) {
      fflush(stdout);
      MPI_Barrier(MPI_COMM_WORLD);
      if (is_master()) {
        va_list args;
        va_start(args, format);
        vfprintf(stdout, format, args);
        fflush(stdout);
        va_end(args);
      }
    }
    template<class T>
    void send(int target_task_id, const T& data) { }
    template<class T>
    void broadcast(const T& data) { }
    template<class T>
    std::list<T> collect(const T& sample) {
      std::list<T> res;
      return res;
    }
    bool is_master() {
      return task_id == 0;
    }
    double sum(const double partial_sum) {
      double total_sum = 0.0;
      MPI_Reduce(&partial_sum, &total_sum,
          1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      return total_sum;
    }
    int sum(const int partial_sum) {
      int total_sum = 0.0;
      MPI_Reduce(&partial_sum, &total_sum,
          1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      return total_sum;
    }
    void finish() {
      MPI_Finalize();
    }
    Parallel(const Parallel&) = delete;
    void operator=(Parallel const&) = delete;
  private:
    Parallel() { }

    int n_tasks;
    int task_id;
    char hostname[MPI_MAX_PROCESSOR_NAME];
};

}

#endif