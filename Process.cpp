#include "Process.h"
#include <climits>

Process::Process(char l, bool cb, int arr_time, int n_bursts, std::vector<int> &cpu, std::vector<int> &io) {
    label = l;
    cpu_bound = cb;
    arrival_time = arr_time;
    num_bursts = n_bursts;
    cpu_bursts = cpu;
    io_bursts = io;
    cpu_availability = false;
    done_running = INT_MAX;
    done_waiting = INT_MAX;
    tau = 0;
    time_ran = 0;
    time_started = 0;
    time_started_waiting = 0;
    time_start_running = -1;
}

void Process::printProcess(std::ostream& output) {
    if (cpu_bound) {
        output << "CPU-bound process " << label << ": " << cpu_bursts.size() << " bursts" << std::endl;
    } else {
        output << "I/O-bound process " << label << ": " << cpu_bursts.size() << " bursts" << std::endl;
    }
}


bool fcfs_compare_arr(Process& p1, Process& p2) {
    if (p1.getArrivalTime() == p2.getArrivalTime()) {
        return p1.getLabel() < p2.getLabel();
    } else {
        return p1.getArrivalTime() < p2.getArrivalTime();
    }
}

bool fcfs_compare_wait(Process& p1, Process& p2) {
    if(p1.getDoneWaiting() == p2.getDoneWaiting()){
        return p1.getLabel() < p2.getLabel();
    }else{
        return p1.getDoneWaiting() < p2.getDoneWaiting();
    }
}

bool sjf_compare_tau(Process&p1, Process& p2){
    if(p1.getTau() == p2.getTau()){
        return p1.getLabel() < p2.getLabel();
    } else {
         return p1.getTau() < p2.getTau();
    }
}

bool srt_compare_time_remaining(Process &p1, Process &p2){
    if (p1.getTau() - p1.getTimeRan() == p2.getTau() - p2.getTimeRan()){
        return p1.getLabel() < p2.getLabel();
    } else {
        return p1.getTau() - p1.getTimeRan() < p2.getTau() - p2.getTimeRan();
    }
}

   