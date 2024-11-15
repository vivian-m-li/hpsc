- only thing that's new to linear system compared to lab9: minus quarter\*rhoCp on all diagonals in Acoef
- on the RHS: b[p] = quarter*rhoCP*phi[p]
  - due to rhoCp\*T_old term in the equation (see notes); phi[p] = T_old
- if restart is set to 1 then the code will read the restart file
- time marching loop:
  - need an MPI barrier at each timestep to make sure everything is in sync
  - form the linear system
- new routine called dump.h which we will have to write in lab
- when doing a dump file, make a decision on which processor writes the boundary nodes
  - multiple ways to do this; which value gets written depends on what your strategy is
  - reading process is different from the writing process
    - the start index of the file read while writing is different than the start index of the file read while reading (the processors have to read the boundary nodes but don't have to write all of them)
- in the dump file:

  - "header" information: inputs: dt, nCellx (total in the global grid), nCelly (total in the global grid)... we want to be able to restart with a different number of processes
  - scalar state variables: simulation time, plot count... how many plots have we written
  - phi[]

- phi = 0 is the "initial condition": initial temp distribution before boundary conditions are applied, then immediately BCs take over on the boundaries

- mySlurm.py: watches for a job and times it, kills it if too much time has passed
- keepRunning.py: restart script
  - executable
    - pathToExe: 'xclock'
  - remove os.system('clear') to put output in the lab
- restart strategy
  - each time a restart is written, we lose computation time. but each time we have to restart, we lose any computations after the last dump
  - advanced capability would query the runtime, compare it against remaining queue time, and remember dump-writing time, and then write at the very end... still a bit risky
- verification of a restart system
  - on this lab and moving forward: verification plots must show both sets of data on one plot. no more side-by-side plots. and also, be sure to use either lines for one and points for the other data set, or thick lines for one and thin lines for the other. do not just have one dataset disappear under the other dataset
- requirements/goals of a restart system
  - goal: have our science code seamlessrun the across resource limitations (time and number of PE available)
  - capabilities needed:
    1. must be able to launch a slurm/batch job
    2. must be able to monitor the job (is it still running? is it done?)
    3. track the restart file creation and status. know what the executable is doing and what the restart files represent, time-wise
    4. manage the entire run, just as you would
  - we will use a slurm "surrogate" that does a nice job of killing your code
  - another goal: debugging
    - imagine a big run, that runs for 5 days on a shared resource. on day 3 it crashes. you write a restart dump one timestep before the crash... or further back as necessary
- ps -elf
  - prints the jobs running on the computer
- xclock
- tee tty_out: writes as the command is writing
- dump.h: get restart capability working in writeRestart - gets called every time we write a gnuplot file and a dump
  -
