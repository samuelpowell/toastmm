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

// ===========================================================================
// class ThreadPool

typedef struct tpool_work {
    void (*routine)(void*,int,int); // task
    void *arg;                      // task arguments
    int idx0, idx1;                 // task index range
    int id;                         // task no. in sequence
    int *counter;                   // pointer to sequence task counter
    struct tpool_work *next;
} tpool_work_t;

typedef struct tpool {
    int num_threads;                // number of threads
    int queue_size;                 // current queue size
    std::thread *threads;             // array of worker threads
    tpool_work_t *queue_head, *queue_tail;

    std::mutex queue_lock;
    std::condition_variable queue_not_empty;
    std::condition_variable thread_done;
} tpool_t;

class ThreadPool {
public:
    ThreadPool (int num_threads);
    // Construct a pool with `num_threads' threads

    void ProcessSequence (void (*routine)(void*,int,int), void *arg,
        int from, int to, int grain=1);
    // Calls `routine(arg,i)' for a sequence from <= i < to of indices
    // `grain' defines the granularity, i.e. the number of indices assigned
    // to each task
    // Function returns when complete sequence is processed

    inline void LockUserMutex() { user_lock.lock(); }
    inline void UnlockUserMutex() { user_lock.unlock(); }

private:
    tpool_t *tpool;              // pool properties
    std::mutex user_lock;   // user-space mutex
};


// ===========================================================================
// class ThreadPool2

typedef struct {
    void(*task)(int,void*);
    void *context;
    std::mutex mutex;  // general-purpose mutex to be used by workers
} THREAD_GLOBAL;

typedef struct {
    int nth;                      // thread index
    std::thread thread;             // thread handle
    std::mutex wakeup_mutex; // locked by worker during task processing
    std::mutex done_mutex;
    std::condition_variable wakeup_cond;
    std::condition_variable done_cond;
    bool wakeup;
    bool done;
    THREAD_GLOBAL *tg;           // pointer to global pool data
} THREAD_DATA;

class ThreadPool2 {
public:
    ThreadPool2 (int num_threads);
    ~ThreadPool2 ();
    static void Initialise (int num_threads);
    static ThreadPool2 *Pool() { return g_tpool2; }

    void MutexLock() { tg.mutex.lock(); }
    void MutexUnlock() { tg.mutex.unlock(); }
    inline int NumThread() const { return nthread; }

    void Invoke (void(*func)(int,void*), void *context);

private:
    int nthread;                   // number of threads
    THREAD_GLOBAL tg;              // global pool data
    THREAD_DATA *td;               // array of worker threads
    static ThreadPool2 *g_tpool2;

    static void *tpool_thread (void *context);
};

#ifndef __TASK_CC
#ifndef OLD_PARALLEL
void TPool_Init (int nt);
extern ThreadPool *g_tpool;
#endif
#endif

#endif // TOAST_THREAD
#endif // !__TASK_H
