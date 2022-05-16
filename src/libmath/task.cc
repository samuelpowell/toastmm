#ifdef TOAST_THREAD

#define __TASK_CC
#define MATHLIB_IMPLEMENTATION
//#define DEBUG_THREAD

#include <mutex>
#include <thread>

#include <sys/types.h>
#include <iostream>
#include <stdio.h>
#include "mathlib.h"
#include "task.h"
#include "timing.h"

using namespace std;

// initialisation of static members

int Task::nthread = Task::nProcessor();
double Task::ttime = 0.0;
double Task::wtime = 0.0;
bool Task::is_multiprocessing = false;

std::mutex Task::user_mutex;

MATHLIB void Task_Init (int nth)
{
    if (!nth) nth = Task::nProcessor();
    Task::SetThreadCount (nth);
    cout << "Toast using " << nth << " threads" << std::endl;
    // int ompnth = omp_get_max_threads();
    // cout << "Toast thread count: " << nth << " OMP thread max: " << ompnth << std::endl;
}

int Task::nProcessor ()
{
    int ntauto = std::thread::hardware_concurrency();
    return (ntauto < 1) ? 1 : ntauto;
}

void Task::Multiprocess (void (*func)(task_data*), void *data, int np)
{
    if (!np) np = nthread;
    int p;
    task_data *td = new task_data[np];
    std::thread *thread = new std::thread[np];

    is_multiprocessing = true;
    LOGOUT2("Multiprocess: Branch %d", np);
    double t0 = tic();
    double w, w0 = walltic();

    // create worker threads
    for (p = 1; p < np; p++) {
        td[p].proc = p;
        td[p].np   = np;
        td[p].data = data;
        thread[p] = std::thread((void*(*)(void*))func, (void*)(td+p));
    }

    // let the master thread do some work
    td[0].proc = 0;
    td[0].np   = np;
    td[0].data = data;
    (*func) (td);

    // wait for workers to finish
    for (p = 1; p < np; p++) {
        thread[p].join();
    }

    ttime += toc(t0);
    wtime += (w=walltoc(w0));
    LOGOUT2("Multiprocess: Join (t=%lf)", w);
    is_multiprocessing = false;
    
    delete []td;
    delete []thread;
}



#endif // TOAST_THREAD