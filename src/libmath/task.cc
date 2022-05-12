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
    cout << "Default thread count set to " << nth << endl;
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

// // ===========================================================================
// // class ThreadPool

// // this is the loop for worker threads waiting for tasks to arrive
// // in the queue

// void tpool_thread (tpool_t *tpool)
// {
//     tpool_work_t *my_workp;

//     for (;;) {
//         tpool->queue_lock.lock();
// 	while (!tpool->queue_size)
// 	    tpool->queue_not_empty.wait(tpool->queue_lock);
// 	my_workp = tpool->queue_head;
// 	tpool->queue_size--;
// 	if (tpool->queue_size == 0)
// 	    tpool->queue_head = tpool->queue_tail = NULL;
// 	else
// 	   tpool->queue_head = my_workp->next;

// #ifdef DEBUG_THREAD
// 	printf ("Waking thread for task %3d of sequence %x\n",
// 		my_workp->id, my_workp->counter);
// #endif

// 	tpool->queue_lock.unlock();
// 	(*(my_workp->routine))(my_workp->arg, my_workp->idx0, my_workp->idx1);
// 	tpool->queue_lock.lock();
// 	(*my_workp->counter)--;

// #ifdef DEBUG_THREAD
// 	printf ("Thread finished task   %3d of sequence %x (%d remaining)\n",
// 		my_workp->id, my_workp->counter, *my_workp->counter);
// #endif
// 	delete my_workp;
// 	tpool->thread_done.notify_all();
// 	tpool->queue_lock.unlock();
//     }
// }

// ThreadPool::ThreadPool (int num_threads)
// {
//     tpool = new tpool_t;
//     tpool->num_threads = num_threads;
//     tpool->queue_size = 0;
//     tpool->queue_head = NULL;
//     tpool->queue_tail = NULL;
//     tpool->threads = new std::thread[num_threads];
//     // if (pthread_mutex_init (&(tpool->queue_lock), NULL))
//     //     cerr << "ThreadPool: pthread_mutex_init failed\n";
//     // if (pthread_mutex_init (&user_lock, NULL))
//     //     cerr << "ThreadPool: pthread_mutex_init failed\n";
//     // if (pthread_cond_init (&(tpool->queue_not_empty), NULL))
//     //     cerr << "ThreadPool: pthread_cond_init failed\n";
//     // if (pthread_cond_init (&(tpool->thread_done), NULL))
//     //     cerr << "ThreadPool: pthread_cond_init failed\n";

//     // create worker threads
//     for (int i = 0; i < num_threads; i++) {
//         tpool->threads[i] = std::thread((void*(*)(void*))tpool_thread, (void*)tpool);
//     }

// #ifdef DEBUG_THREAD
//     cerr << "Created thread pool with " << num_threads << " threads\n";
// #endif
// }

// void ThreadPool::ProcessSequence (void (*routine)(void*,int,int), void *arg,
//     int from, int to, int grain)
// {
//     tpool_work_t *workp;
//     int n, *task_count;

//     tpool->queue_lock.lock();
//     if (grain <= 0) grain = 1;
//     n = (to-from+grain-1)/grain;
//     if (n <= 0) { // nothing to do
//         tpool->queue_lock.unlock();
// 	return;
//     } else {
//         task_count = new int; // counts tasks remaining for this sequence
// 	*task_count = n;
//     }

//     // add loops to queue
//     for (int i = 0; i < n; i++) {
//         workp = new tpool_work_t;
// 	workp->routine = routine;
// 	workp->arg     = arg;
// 	if ((workp->idx1 = (workp->idx0 = i*grain+from) + grain) > to)
// 	    workp->idx1 = to;
// 	workp->id      = i;
// 	workp->counter = task_count;
// 	workp->next    = NULL;
// 	if (tpool->queue_tail) tpool->queue_tail->next = workp;
// 	else                   tpool->queue_head = workp;
// 	tpool->queue_tail = workp;
//     }

//     tpool->queue_size += n;
//     if (tpool->queue_size == n) // wake up dormant threads
//         tpool->queue_not_empty.notify_all();

// #ifdef DEBUG_THREAD
//     cerr << "Added sequence " << task_count << " (" << n
// 	 << " tasks) to queue, queue size now " << tpool->queue_size << endl;
// #endif // DEBUG_THREAD

//     // wait until sequence counter is down to zero
//     while (*task_count)
//         tpool->thread_done.wait(tpool->queue_lock);

// #ifdef DEBUG_THREAD
//     cerr << "Sequence finished, queue size now " << tpool->queue_size << endl;
//     cerr << "Returning from ProcessSequence\n";
// #endif

//     delete task_count;
//     tpool->queue_lock.unlock();
// }

// #ifndef OLD_PARALLEL
// ThreadPool *g_tpool = 0;

// void TPool_Init (int nt)
// {
//     g_tpool = new ThreadPool (nt);
// }

// #endif



// // ===========================================================================
// // class ThreadPoo2

// ThreadPool2 *ThreadPool2::g_tpool2 = NULL;

// ThreadPool2::ThreadPool2 (int num_threads)
// {
//     nthread = num_threads;
//     td = new THREAD_DATA[nthread];

//     tg.task = NULL;
//     pthread_mutex_init (&tg.mutex, NULL);

//     for (int i = 0; i < nthread; i++) {
//         td[i].tg = &tg;
// 	td[i].nth = i;
// 	td[i].wakeup = false;
// 	td[i].done = false;
// 	pthread_mutex_init (&td[i].done_mutex, NULL);
// 	pthread_mutex_init (&td[i].wakeup_mutex, NULL);
// 	pthread_cond_init (&td[i].wakeup_cond, NULL);
// 	pthread_cond_init (&td[i].done_cond, NULL);
// 	pthread_mutex_lock (&td[i].done_mutex);

//         if (pthread_create (&td[i].thread, NULL,
// 	    (void*(*)(void*))tpool_thread, (void*)(td+i)))
// 	        cerr << "TreadPool: pthread_create failed" << endl;
//     }

// }

// ThreadPool2::~ThreadPool2 ()
// {
//     // should destroy threads and mutexes here
// }

// void ThreadPool2::Initialise (int num_threads)
// {
//     g_tpool2 = new ThreadPool2 (num_threads);
// }

// void ThreadPool2::Invoke (void(*func)(int,void*), void *context)
// {
//     int i;

//     // set the task function
//     tg.task = func;
//     tg.context = context;

//     // wake the threads
//     for (i = 0; i < nthread; i++) {
// 	td[i].wakeup = true;
// 	pthread_cond_signal (&td[i].wakeup_cond);
//     }
//     for (i = 0; i < nthread; i++) {
//         while (!td[i].done)
//   	    pthread_cond_wait (&td[i].done_cond, &td[i].done_mutex);

// 	td[i].done = false;
//     }
//     //cerr << "threads done" << endl;
// }

// void *ThreadPool2::tpool_thread (void *context)
// {
//     THREAD_DATA *td = (THREAD_DATA*)context;
//     THREAD_GLOBAL *tg = td->tg;
//     pthread_mutex_lock (&td->wakeup_mutex);
//     //cerr << "worker: locked wakeup" << endl;

//     // worker loop
//     for (;;) {
// 	// wait for a task
// 	while (!td->wakeup)
// 	    pthread_cond_wait (&td->wakeup_cond, &td->wakeup_mutex);
// 	td->wakeup = false;

// 	// process the task
// 	if (tg->task)
// 	    (*(tg->task))(td->nth, tg->context);

// 	// signal finished
// 	td->done = true;
// 	pthread_cond_signal (&td->done_cond);
//     }
//     return NULL;
// }

#endif // TOAST_THREAD