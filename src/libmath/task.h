// -*-C++-*-

#ifndef __TASK_H
#define __TASK_H

#ifdef TOAST_THREAD

#include <thread>
#include <mutex>
#include <condition_variable>

#define NUMTHREAD 2

typedef struct {
    int proc;
    int np;
    void *data;
} task_data;

MATHLIB void Task_Init (int nth);

class MATHLIB Task {
public:
    Task() {}

    static int nProcessor();
    // number of processors available

    static void SetThreadCount (int _nthread) { nthread = _nthread; }
    // set the number of threads to use for multiprocessing tasks.
    // Default is nProcessor()

    static int GetThreadCount () { return nthread; }
    // return current thread count setting

    static double GetThreadCPUTiming () { return ttime; }
    // returns user time [sec] spent by threads inside Multiprocess().
    // Note that the interpretation of this time is platform-dependent:
    // For Linux, this is the time spent by the master thread. For Sun and
    // SGI, it is the sum of all threads.

    static double GetThreadWallTiming () { return wtime; }
    // Wall clock (real) time spent inside Multiprocess()

    static void Multiprocess (void (*func)(task_data*), void *data, int np = 0);
    // run function 'func' in parallel. User data 'data' are passed to
    // each instance of 'func'. 'np' is the number of threads to create. If
    // np==0 then nthread is used

    static bool IsMultiprocessing() { return is_multiprocessing; }

    inline static void UserMutex_lock() { user_mutex.lock(); }

    inline static void UserMutex_unlock() { user_mutex.unlock(); }

private:
    static int nthread;
    static double ttime; // accumulated cpu time spent multiprocessing
    static double wtime; // accumulated wall clock time spent multiprocessing
    static std::mutex user_mutex;
    static bool is_multiprocessing;  // we are currently in Multiprocess
};



#endif // TOAST_THREAD
#endif // !__TASK_H
