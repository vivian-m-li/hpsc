10/15

- we finally need omp.h and the timer function
- print myTH (my thread) instead of myPE
- omp_get_max_threads()
- omp_set_num_threads(4)
- omp_get_thread_num()
- #pragma omp parallel
  - decoration that launches a parallel region: no longer running only on Core0
- startTimer(t0), endTimer(str, t0, t1)
- #pragma omp parallel for
  - tell omp this is a for loop and split up the work
  - fork join process takes time so more cores is not always better
- this week's lab: use multiple threads to populate the matrix and solve the linear system
- #pragma omp parallel for firstprivate(FPvar)
  - now every thread gets the value of FPvar

- for timer: if it goes over a minute, then the seconds can roll over (in omp.h)