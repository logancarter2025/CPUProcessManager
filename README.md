
# CPU Process Manager

This project simulates an operating system given the set of processes pseudo-randomly generated. The overall focus will be on processes, assumed to be resident in memory, waiting to use the CPU.

This program implements four algorithms to determine the CPU usage. When the program is ran, all four algorithms are to be simulated in succession with the same set of processes. Each algorithm is summarized below.

First-come-first-served (FCFS)
The FCFS algorithm is a non-preemptive algorithm in which processes simply line up in the ready
queue, waiting to use the CPU. This is the baseline algorithm.

Shortest job first (SJF)
In SJF, processes are stored in the ready queue in order of priority based on their anticipated CPU burst times. More specifically, the process with the shortest predicted CPU burst time will be
selected as the next process executed by the CPU.

Shortest remaining time (SRT)
The SRT algorithm is a preemptive version of the SJF algorithm. In SRT, when a process arrives, if it has a predicted CPU burst time that is less than the remaining predicted time of the currently
running process, a preemption occurs. When such a preemption occurs, the currently running process is simply added to the ready queue.

Round robin (RR)
The RR algorithm is essentially the FCFS algorithm with time slice tslice. Each process is given tslice amount of time to complete its CPU burst. If the time slice expires, the process is preempted and added to the end of the ready queue.
If a process completes its CPU burst before a time slice expiration, the next process on the ready queue is context-switched in to use the CPU.
For the simulation, if a preemption occurs and there are no other processes on the ready queue,
a context switch is not performed. For example, given process G is using the CPU and the ready
queue is empty, if process G is preempted by a time slice expiration, do not context-switch process G back to the empty queue; instead, keep process G running with the CPU and do not count this as a context switch. In other words, when the time slice expires, the queue is checked to determine if a context switch should occur.

## Acknowledgements

 - Worked in collaboration with Anthony Basnight and Isaac Efrosman
## Usage/Examples

```c++
To compile:
g++ -g main.cpp Process.cpp
execute:
./a.out 3 1 1024 0.001 3000 4 0.75 256

argv[1]: Define n as the number of processes to simulate. Process IDs are assigned in
alphabetical order A through Z. Therefore, you will have at most 26 processes to simulate.

argv[2]: Define ncpu as the number of processes that are CPU-bound. For this project, we
will classify processes as I/O-bound or CPU-bound. The ncpu CPU-bound processes, when
generated, will have CPU burst times that are longer by a factor of 4 and will have I/O burst
times that are shorter by a factor of 8.

argv[3]: This command-line argument serves as the seed for the pseudorandom number sequence. To ensure
predictability and repeatability, I use srand48() with this given seed before simulating each scheduling
algorithm and drand48() to obtain the next value in the range [0.0, 1.0).

argv[4]: To determine interarrival times, we will use an exponential distribution; therefore, this
command-line argument is parameter λ. Remember that 1 λ will be the average random value generated,
e.g., if λ = 0.01, then the average should be appoximately 100. See the the formula shown in the code,
i.e., -log( r )  lambda, where log is the natural logarithm.

argv[5]: For the exponential distribution, this command-line argument represents the upper bound for
valid pseudo-random numbers. This threshold is used to avoid values far down the long tail of the
exponential distribution. As an example, if this is set to 3000, all generated values above 3000
should be skipped.

argv[6]: Define tcs as the time, in milliseconds, that it takes to perform a context switch. Specifically,
the first half of the context switch time (i.e., tcs 2) is the time required to remove the given process
from the CPU; the second half of the context switch time is the time required to bring the next process
in touse the CPU. Therefore, expect tcs to be a positive even integer.

argv[7]: For the SJF and SRT algorithms, since we cannot know the actual CPU burst times beforehand,
we will rely on estimates determined via exponential averaging. As such, this command-line argument
is the constant α. Note that the initial guess for each process is τ0 = 1 λ.

argv[8]: For the RR algorithm, define the time slice value, tslice, measured in milliseconds.



```


## Authors

- [Logan Carter](https://github.com/logancarter2025)
- [Anthony Basnight](https://github.com/anthony-basnight)
- [Isaac Efrosman](https://github.com/IsaacEf)



## License

[MIT](https://choosealicense.com/licenses/mit/)

