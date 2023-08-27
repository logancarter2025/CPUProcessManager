#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <math.h>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include "Process.h"

double next_exp(double& lambda, int& threshold) {
    double r = drand48();
    double x = -log(r)/lambda;
    if (x > threshold) {
        x = next_exp(lambda, threshold);
    }
    return x;
}

void fcfs(std::vector<Process> processes, int context_switch, std::ofstream &output, clock_t start) {

    //Variables used for analysis to output file
    int timeCpuRunning = 0;
    int cpuCS = 0;
    int ioCS = 0;
    int sumCpuTime_CPU = 0;
    int numCpuBursts_CPU = 0;
    int sumCpuTime_IO = 0;
    int numCpuBursts_IO = 0;
    int sumWaitTime_CPU = 0;
    int sumWaitTime_IO = 0;

    int time = 0;
    std::cout << "time " << time << "ms: Simulator started for FCFS [Q <empty>]" << std::endl;

    std::sort(processes.begin(), processes.end(), fcfs_compare_arr);
    int next_arrival = 0;

    std::list<Process> running;
    std::vector<Process> waiting;
    std::list<Process> ready;
    std::map<char, int> cpu_burst_index;
    std::map<char, int> io_burst_index;
    std::vector<Process> terminated;

    for (unsigned int i = 0; i < processes.size(); i++) {
        cpu_burst_index[processes[i].getLabel()] = 0;
        io_burst_index[processes[i].getLabel()] = 0;
        std::vector<int> temp_cpu_bursts = processes[i].getCpuBursts();
        std::vector<int> temp_io_bursts = processes[i].getIoBursts();
        
        if (processes[i].getCpuBound()) { 
            for (unsigned int j = 0; j < temp_cpu_bursts.size(); j++) { sumCpuTime_CPU += temp_cpu_bursts[j]; }
            numCpuBursts_CPU += temp_cpu_bursts.size();
        } else {
            for (unsigned int j = 0; j < temp_cpu_bursts.size(); j++) { sumCpuTime_IO += temp_cpu_bursts[j]; }
            numCpuBursts_IO += temp_cpu_bursts.size();
        }
    }
    clock_t end;
    
    while (processes.size() > terminated.size()) { //Remove process from processes once finished, not arrived
        // Figure out what happens next: arrival, finishing waiting state, or finishing running state
        end = clock();
        if (end - start > 100000){
            break;
        }
        // next_arrival is the index of the next process to arrive from processes

        /* the time when the next process arrives */
        int next_arrival_time = INT_MAX;
        if ((unsigned int)next_arrival < processes.size()) {
            next_arrival_time = processes[next_arrival].getArrivalTime();
        }

        /* the time when the burst using the CPU finishes and is ready to leave the running queue */
        int done_running_time = INT_MAX;
        if (running.size() > 0) {
            done_running_time = running.front().getDoneRunning();
        }

        /* the time when the burst blocking on IO finishes and is ready to leave the waiting queue */
        int done_blocking_time = INT_MAX;
        if (waiting.size() > 0) {
            done_blocking_time = waiting.front().getDoneWaiting();
        }

        if (next_arrival_time <= done_running_time && next_arrival_time <= done_blocking_time) {
            /* NEW ARRIVAL */
            Process temp_process = processes[next_arrival];
            
            time = temp_process.getArrivalTime();
            temp_process.setTimeStarted(time);
            
            temp_process.setTimeStartedWaiting(time);                        
            ready.push_back(temp_process);

            next_arrival += 1;
            if (time < 10000){
                std::cout << "time " << time << "ms: Process " << temp_process.getLabel()
                      << " arrived; added to ready queue [Q";
                if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                else {
                    for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                        std::cout << " " << it->getLabel();
                    }
                    std::cout << "]" << std::endl;
                }
            }
        } else if (done_running_time < next_arrival_time && done_running_time <= done_blocking_time) {
            /* DONE WITH RUNNING STATE */
            
            Process temp_process = running.front();
            running.pop_front();
                        
            time = temp_process.getDoneRunning();

            //Updates how much time CPU has been in use
            temp_process.setDoneRunning(INT_MAX);
            timeCpuRunning += time - temp_process.getTimeStarted();

            if (temp_process.getCpuBursts().size() == (unsigned int)cpu_burst_index[temp_process.getLabel()]) {
                std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " terminated [Q";
                if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                else {
                    for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                        std::cout << " " << it->getLabel();
                    }
                    std::cout << "]" << std::endl;
                }
                terminated.push_back(temp_process);

                time += (context_switch/2);

            } else {

                if (time < 10000){
                    std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " completed a CPU burst; "
                            << temp_process.getCpuBursts().size() - cpu_burst_index[temp_process.getLabel()]
                            << " bursts to go [Q";
                    if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                    else {
                        for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                            std::cout << " " << it->getLabel();
                        }
                        std::cout << "]" << std::endl;   
                    }
                }

                int temp_burst_len = temp_process.getIOBurst(io_burst_index[temp_process.getLabel()]);
                temp_process.setDoneWaiting(time + temp_burst_len + (context_switch / 2));
                
                if (time < 10000){
                    std::cout << "time " << time << "ms: Process " << temp_process.getLabel()
                          << " switching out of CPU; blocking on I/O until time "
                          << time + temp_burst_len + (context_switch / 2) << "ms [Q";
                    if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                    else {
                        for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                            std::cout << " " << it->getLabel();
                        }
                        std::cout << "]" << std::endl;
                    }
                }
                io_burst_index[temp_process.getLabel()] += 1;
                waiting.push_back(temp_process);
                std::sort(waiting.begin(), waiting.end(), fcfs_compare_wait);
                time += (context_switch / 2);
            }
            
            if (temp_process.getCpuBound()) {
                cpuCS++;
            } else {
                ioCS++;
            }

        } else {
            /* DONE WITH WAITING STATE */
            Process temp_process = waiting.front();
            waiting.erase(waiting.begin());

            time = temp_process.getDoneWaiting();
            temp_process.setDoneWaiting(INT_MAX);
            
            temp_process.setTimeStartedWaiting(time);
            ready.push_back(temp_process);

            if (time < 10000) {
                std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " completed I/O; added to ready queue [Q";
                if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                else {
                    for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                        std::cout << " " << it->getLabel();
                    }
                    std::cout << "]" << std::endl;
                }
            }          
        }

        if (running.size() == 0 && ready.size() > 0) {
            /* MOVE NEXT PROCESS FROM READY STATE TO RUNNING STATE */
            
            Process temp_process = ready.front();
            ready.pop_front();

            if (temp_process.getCpuBound()) {
                sumWaitTime_CPU += time - temp_process.getTimeStartedWaiting();
            } else {
                sumWaitTime_IO += time - temp_process.getTimeStartedWaiting();
            }
        
            time += (context_switch / 2);

            int temp_burst_len = temp_process.getCpuBurst(cpu_burst_index[temp_process.getLabel()]);
            temp_process.setDoneRunning(time + temp_burst_len);
            temp_process.setTimeStarted(time);
            
            if (time < 10000) {
                std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " started using the CPU for "
                        << temp_burst_len << "ms burst [Q";
                if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                else {
                    for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                        std::cout << " " << it->getLabel();
                    }
                    std::cout << "]" << std::endl;
                }
            }

            cpu_burst_index[temp_process.getLabel()] += 1;
            running.push_back(temp_process);
        }
    }
    std::cout << "time " << time << "ms: Simulator ended for FCFS [Q <empty>]" << std::endl;


    double cpu_utilization = 0;
    if (time != 0){
        cpu_utilization = ceil(1000 * ((double)timeCpuRunning / time) * 100) / 1000;
    }
    
    double cpuBurstAvg_CPU = 0;
    if (numCpuBursts_CPU != 0){
        cpuBurstAvg_CPU = ceil(1000 * (double)sumCpuTime_CPU / numCpuBursts_CPU) / 1000;
    }
    
    double cpuBurstAvg_IO = 0;
    if (numCpuBursts_IO != 0){
        cpuBurstAvg_IO = ceil(1000 * (double)sumCpuTime_IO / numCpuBursts_IO) / 1000;
    }

    double cpuBurstAvg_SUM = 0;
    if (numCpuBursts_CPU + numCpuBursts_IO != 0){
        cpuBurstAvg_SUM = ceil(1000.0 * (double)(sumCpuTime_CPU + sumCpuTime_IO) / (numCpuBursts_CPU + numCpuBursts_IO)) / 1000;
    }
    
    double IOBurstAvg_CPU = 0;
    if (numCpuBursts_CPU != 0){
        IOBurstAvg_CPU = ceil(1000 * (double)sumWaitTime_CPU / numCpuBursts_CPU ) / 1000;
    }

    double IOBurstAvg_IO = 0;
    if (numCpuBursts_IO != 0){
        IOBurstAvg_IO = ceil(1000 * (double)sumWaitTime_IO / numCpuBursts_IO ) / 1000;
    }
    double IOBurstAvg_SUM = 0;
    if (numCpuBursts_CPU + numCpuBursts_IO != 0){
        IOBurstAvg_SUM = ceil(1000.0 * (double)(sumWaitTime_CPU + sumWaitTime_IO) / (numCpuBursts_CPU + numCpuBursts_IO)) / 1000;
    }

 
    output << "Algorithm FCFS\n";
    output << "-- CPU utilization: " << std::fixed << std::setprecision(3) << cpu_utilization << "%\n";
    output << "-- average CPU burst time: " << cpuBurstAvg_SUM << " ms (" << cpuBurstAvg_CPU << " ms/" << cpuBurstAvg_IO << " ms)\n";
    output << "-- average wait time: " << IOBurstAvg_SUM << " ms (" << IOBurstAvg_CPU << " ms/" << IOBurstAvg_IO << " ms)\n";
    output << "-- average turnaround time: " << cpuBurstAvg_SUM + IOBurstAvg_SUM + context_switch << " ms (" << cpuBurstAvg_CPU + IOBurstAvg_CPU + context_switch 
           << " ms/" << cpuBurstAvg_IO + IOBurstAvg_IO + context_switch << " ms)\n";
    output << "-- number of context switches: " << (cpuCS + ioCS) << " (" << cpuCS << "/" << ioCS << ")\n";
    output << "-- number of preemptions: 0 (0/0)\n\n";

}

void sjf(std::vector<Process> processes, int context_switch, double lambda, double alpha, std::ofstream &output, clock_t start){
    int time = 0;
    std::cout << "time " << time << "ms: Simulator started for SJF [Q <empty>]" << std::endl;

    //sort processes by their arrival time
    std::sort(processes.begin(), processes.end(), fcfs_compare_arr);
    int next_arrival = 0;
    
    int timeCpuRunning = 0; 
    int cpuCS = 0;
    int ioCS = 0;
    int sumCpuTime_CPU = 0;
    int numCpuBursts_CPU = 0;
    int sumCpuTime_IO = 0;
    int numCpuBursts_IO = 0;
    int sumWaitTime_CPU = 0;
    int sumWaitTime_IO = 0;


    
    std::list<Process> running;
    std::vector<Process> waiting;
    std::vector<Process> ready;
    
    //Maps process id ('A', 'B', ... 'Z') to how many cpu bursts have been completed
    std::map<char, int> cpu_burst_index;

    //Maps process id ('A', 'B', ... 'Z') to how many io bursts have been completed
    std::map<char, int> io_burst_index;

    //Vector of processes which have finished
    std::vector<Process> terminated;

    //Initialize tau values, and maps. At beginning of simulation, no bursts have been completed
    for (unsigned int i = 0; i < processes.size(); i++) {
        cpu_burst_index[processes[i].getLabel()] = 0;
        io_burst_index[processes[i].getLabel()] = 0;
        processes[i].setTau(ceil(1 / lambda));
        
        std::vector<int> temp_cpu_bursts = processes[i].getCpuBursts();
        std::vector<int> temp_io_bursts = processes[i].getIoBursts();
        
        if (processes[i].getCpuBound()) { 
            for (unsigned int j = 0; j < temp_cpu_bursts.size(); j++) { sumCpuTime_CPU += temp_cpu_bursts[j]; }
            numCpuBursts_CPU += temp_cpu_bursts.size();
        } else {
            for (unsigned int j = 0; j < temp_cpu_bursts.size(); j++) { sumCpuTime_IO += temp_cpu_bursts[j]; }
            numCpuBursts_IO += temp_cpu_bursts.size();
        }
    }

    clock_t end;
    while (processes.size() > terminated.size()) { //Remove process from processes once finished, not arrived

        end = clock();
        if (end - start > 100000){
            break;
        }
        
        // Figure out what happens next: arrival, finishing waiting state, or finishing running state

        // next_arrival is the index of the next process to arrive from processes

        /* the time when the next process arrives */
        int next_arrival_time = INT_MAX;
        if ((unsigned int)next_arrival < processes.size()) {
            next_arrival_time = processes[next_arrival].getArrivalTime();
        }

        /* the time when the burst using the CPU finishes and is ready to leave the running queue */
        int done_running_time = INT_MAX;
        if (running.size() > 0) {
            done_running_time = running.front().getDoneRunning();
        }

        /* the time when the burst blocking on IO finishes and is ready to leave the waiting queue */
        int done_blocking_time = INT_MAX;
        if (waiting.size() > 0) {
            done_blocking_time = waiting.front().getDoneWaiting();
        }

        if (next_arrival_time <= done_running_time && next_arrival_time <= done_blocking_time) {
            /* NEW ARRIVAL */
            Process temp_process = processes[next_arrival];
        
            time = temp_process.getArrivalTime();
            temp_process.setTimeStartedWaiting(time);                        

            ready.push_back(temp_process);
            std::sort(ready.begin(), ready.end(), sjf_compare_tau);
            next_arrival += 1;

            if (time < 10000) { 
                std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " (tau " << (int)ceil(temp_process.getTau()) << "ms) arrived; added to ready queue [Q"; 
                if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                else {
                    for (std::vector<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                        std::cout << " " << it->getLabel();
                    }
                    std::cout << "]" << std::endl;
                }
            }
        } else if (done_running_time < next_arrival_time && done_running_time <= done_blocking_time) {
            /* DONE WITH RUNNING STATE */
            Process temp_process = running.front();
            running.pop_front();
            

            time = temp_process.getDoneRunning();
            temp_process.setDoneRunning(INT_MAX);
            timeCpuRunning += time - temp_process.getTimeStarted();

            if (temp_process.getCpuBursts().size() == (unsigned int)cpu_burst_index[temp_process.getLabel()]) {
                
                
                std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " terminated [Q";
                if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                else {
                    for (std::vector<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                        std::cout << " " << it->getLabel();
                    }
                    std::cout << "]" << std::endl;
                }
                
                time += (context_switch / 2);

                terminated.push_back(temp_process);
            } else {
                if (time < 10000) {
                    std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " (tau " << (int)ceil(temp_process.getTau()) << "ms) completed a CPU burst; "
                            << temp_process.getCpuBursts().size() - cpu_burst_index[temp_process.getLabel()] << " bursts to go [Q";
                    if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                    else {
                        for (std::vector<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                            std::cout << " " << it->getLabel();
                        }
                        std::cout << "]" << std::endl;
                    }
                }

                double old_tau = temp_process.getTau();
                
                temp_process.setTau(ceil(alpha * temp_process.getCpuBurst(cpu_burst_index[temp_process.getLabel()] - 1) + ( 1.0 - alpha ) * (int)old_tau));

                if (time < 10000) { 
                    std::cout << "time " << time <<  "ms: Recalculating tau for process " << temp_process.getLabel() << ": old tau " << (int)old_tau <<
                              "ms ==> new tau " << (int)ceil(temp_process.getTau()) << "ms [Q"; 
                    if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                    else {
                        for (std::vector<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                            std::cout << " " << it->getLabel();
                        }
                        std::cout << "]" << std::endl;
                    }
                }


                int temp_burst_len = temp_process.getIOBurst(io_burst_index[temp_process.getLabel()]);
                temp_process.setDoneWaiting(time + temp_burst_len + (context_switch / 2));
                if (time < 10000) {
                    std::cout << "time " << time << "ms: Process " << temp_process.getLabel()
                            << " switching out of CPU; blocking on I/O until time "
                            << time + temp_burst_len + (context_switch / 2) << "ms [Q";
                    if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                    else {
                        for (std::vector<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                            std::cout << " " << it->getLabel();
                        }
                        std::cout << "]" << std::endl;
                    }
                }


                if (temp_process.getCpuBound()) {
                    cpuCS++;
                } else {
                    ioCS++;
                }

                io_burst_index[temp_process.getLabel()] += 1;
                waiting.push_back(temp_process);
                std::sort(waiting.begin(), waiting.end(), fcfs_compare_wait);
                time += (context_switch / 2);

            }
            
        } else {
            /* DONE WITH WAITING STATE */
            Process temp_process = waiting.front();
            waiting.erase(waiting.begin());

            time = temp_process.getDoneWaiting();
            temp_process.setDoneWaiting(INT_MAX);
            temp_process.setTimeStartedWaiting(time);
            
            ready.push_back(temp_process);
            std::sort(ready.begin(), ready.end(), sjf_compare_tau);

            if (time < 10000) {
                std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " (tau " << (int)ceil(temp_process.getTau()) << "ms) completed I/O; added to ready queue [Q";
                if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                else {
                    for (std::vector<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                        std::cout << " " << it->getLabel();
                    }
                    std::cout << "]" << std::endl;
                }
            }
        }

        if (running.size() == 0 && ready.size() > 0) {
            /* MOVE NEXT PROCESS FROM READY STATE TO RUNNING STATE */
            Process temp_process = ready.front();
            ready.erase(ready.begin());

            if (temp_process.getCpuBound()) {
                sumWaitTime_CPU += time - temp_process.getTimeStartedWaiting();
            } else {
                sumWaitTime_IO += time - temp_process.getTimeStartedWaiting();
            }

            time += (context_switch / 2);

            int temp_burst_len = temp_process.getCpuBurst(cpu_burst_index[temp_process.getLabel()]);
            temp_process.setDoneRunning(time + temp_burst_len);
            if (time < 10000) {
                std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " (tau " << (int)ceil(temp_process.getTau()) << "ms) started using the CPU for "
                        << temp_burst_len << "ms burst [Q";
                if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                else {
                    for (std::vector<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                        std::cout << " " << it->getLabel();
                    }
                    std::cout << "]" << std::endl;
                }
            }

            temp_process.setTimeStarted(time);
            cpu_burst_index[temp_process.getLabel()] += 1;
            running.push_back(temp_process);
        }
    }
    std::cout << "time " << time << "ms: Simulator ended for SJF [Q <empty>]" << std::endl;

    double cpu_utilization = 0;
    if (time != 0){
        cpu_utilization = ceil(1000 * ((double)timeCpuRunning / time) * 100) / 1000;
    }
    
    double cpuBurstAvg_CPU = 0;
    if (numCpuBursts_CPU != 0){
        cpuBurstAvg_CPU = ceil(1000 * (double)sumCpuTime_CPU / numCpuBursts_CPU) / 1000;
    }
    
    double cpuBurstAvg_IO = 0;
    if (numCpuBursts_IO != 0){
        cpuBurstAvg_IO = ceil(1000 * (double)sumCpuTime_IO / numCpuBursts_IO) / 1000;
    }
    
    double cpuBurstAvg_SUM = 0;
    if (numCpuBursts_CPU + numCpuBursts_IO != 0){
        cpuBurstAvg_SUM = ceil(1000.0 * (double)(sumCpuTime_CPU + sumCpuTime_IO) / (numCpuBursts_CPU + numCpuBursts_IO)) / 1000;
    }
     
    double IOBurstAvg_CPU = 0;
    if(numCpuBursts_CPU != 0){
        IOBurstAvg_CPU = ceil(1000 * (double)sumWaitTime_CPU / numCpuBursts_CPU ) / 1000;
    }

    double IOBurstAvg_IO = 0;
    if(numCpuBursts_IO != 0){
        IOBurstAvg_IO = ceil(1000 * (double)sumWaitTime_IO / numCpuBursts_IO ) / 1000;
    }

    double IOBurstAvg_SUM = 0;
    if(numCpuBursts_CPU + numCpuBursts_IO != 0){
        IOBurstAvg_SUM = ceil(1000.0 * (double)(sumWaitTime_CPU + sumWaitTime_IO) / (numCpuBursts_CPU + numCpuBursts_IO)) / 1000;
    }
    
    double turnaround = 0;
    if(numCpuBursts_CPU + numCpuBursts_IO != 0){
        turnaround = (double)(sumWaitTime_CPU + sumWaitTime_IO) / (numCpuBursts_CPU + numCpuBursts_IO) 
                        + (double)(sumCpuTime_CPU + sumCpuTime_IO) / (numCpuBursts_CPU + numCpuBursts_IO) + context_switch;
    }
    
    output << "Algorithm SJF\n";
    output << "-- CPU utilization: " << std::fixed << std::setprecision(3) << cpu_utilization << "%\n";
    output << "-- average CPU burst time: " << cpuBurstAvg_SUM << " ms (" << cpuBurstAvg_CPU << " ms/" << cpuBurstAvg_IO << " ms)\n";
    output << "-- average wait time: " << IOBurstAvg_SUM << " ms (" << IOBurstAvg_CPU << " ms/" << IOBurstAvg_IO << " ms)\n";
    output << "-- average turnaround time: " << ceil(1000 * turnaround) / 1000 << " ms (" << cpuBurstAvg_CPU + IOBurstAvg_CPU + context_switch 
           << " ms/" << cpuBurstAvg_IO + IOBurstAvg_IO + context_switch << " ms)\n";
    output << "-- number of context switches: " << (numCpuBursts_CPU + numCpuBursts_IO) << " (" << numCpuBursts_CPU << "/" << numCpuBursts_IO << ")\n";
    output << "-- number of preemptions: 0 (0/0)\n\n";
   

}

void srt(std::vector<Process> processes, int context_switch, double lambda, double alpha, std::ofstream &output, clock_t start){

    int timeCpuRunning = 0; 

    int sumCpuTime_CPU = 0;
    int numCpuBursts_CPU = 0;
    int sumCpuTime_IO = 0;
    int numCpuBursts_IO = 0;
    int sumWaitTime_CPU = 0;
    int sumWaitTime_IO = 0;
    int num_preemptions_CPU = 0;
    int num_preemptions_IO = 0;

    int time = 0;
    std::cout << "time " << time << "ms: Simulator started for SRT [Q <empty>]" << std::endl;

    //sort processes by their arrival time
    std::sort(processes.begin(), processes.end(), fcfs_compare_arr);
    int next_arrival = 0;

    std::list<Process> running;
    std::vector<Process> waiting;
    std::vector<Process> ready;
    
    //Maps process id ('A', 'B', ... 'Z') to how many cpu bursts have been completed
    std::map<char, int> cpu_burst_index;

    //Maps process id ('A', 'B', ... 'Z') to how many io bursts have been completed
    std::map<char, int> io_burst_index;

    //Vector of processes which have finished
    std::vector<Process> terminated;

    //Initialize tau values, and maps. At beginning of simulation, no bursts have been completed
    for (unsigned int i = 0; i < processes.size(); i++) {
        cpu_burst_index[processes[i].getLabel()] = 0;
        io_burst_index[processes[i].getLabel()] = 0;
        processes[i].setTau(ceil(1 / lambda));

        std::vector<int> temp_cpu_bursts = processes[i].getCpuBursts();
        std::vector<int> temp_io_bursts = processes[i].getIoBursts();
        
        if (processes[i].getCpuBound()) { 
            for (unsigned int j = 0; j < temp_cpu_bursts.size(); j++) { sumCpuTime_CPU += temp_cpu_bursts[j]; }
            numCpuBursts_CPU += temp_cpu_bursts.size();
        } else {
            for (unsigned int j = 0; j < temp_cpu_bursts.size(); j++) { sumCpuTime_IO += temp_cpu_bursts[j]; }
            numCpuBursts_IO += temp_cpu_bursts.size();
        }  
    }

    int temp_total_cpu_time = sumCpuTime_CPU;
    
    clock_t end;
    unsigned int counter = 0;
    while (processes.size() > counter) { //Remove process from processes once finished, not arrived
        end = clock();
        if (end - start > 100000){
            break;
        }
        


        // Figure out what happens next: arrival, finishing waiting state, or finishing running state

        // next_arrival is the index of the next process to arrive from processes

        /* the time when the next process arrives */
        int next_arrival_time = INT_MAX;
        if ((unsigned int)next_arrival < processes.size()) {
            next_arrival_time = processes[next_arrival].getArrivalTime();
        }

        /* the time when the burst using the CPU finishes and is ready to leave the running queue */
        int done_running_time = INT_MAX;
        if (running.size() > 0) {
            done_running_time = running.front().getDoneRunning();
        }

        /* the time when the burst blocking on IO finishes and is ready to leave the waiting queue */
        int done_blocking_time = INT_MAX;
        if (waiting.size() > 0) {
            done_blocking_time = waiting.front().getDoneWaiting();
        }

        if (next_arrival_time <= done_running_time && next_arrival_time <= done_blocking_time && next_arrival_time != INT_MAX) {
            /* NEW ARRIVAL */
            /* PREEMPTING MAY OCCUR */

            Process temp_process = processes[next_arrival];

            time = temp_process.getArrivalTime();

            temp_process.setTimeStarted(time);

            temp_process.setTimeStartedWaiting(time);
            ready.push_back(temp_process);
            std::sort(ready.begin(), ready.end(), srt_compare_time_remaining);
            next_arrival += 1;

            if (time < 10000) { 
                std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " (tau " << (int)ceil(temp_process.getTau()) << "ms) arrived; added to ready queue [Q"; 
                if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                else {
                    for (std::vector<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                        std::cout << " " << it->getLabel();
                    }
                    std::cout << "]" << std::endl;
                }
            }
        } else if (done_running_time < next_arrival_time && done_running_time <= done_blocking_time) {
            /* DONE WITH RUNNING STATE */
            /* no preempting */
            Process temp_process = running.front();
            running.pop_front();

            time = temp_process.getDoneRunning();
            temp_process.setDoneRunning(INT_MAX);
            timeCpuRunning += time - temp_process.getTimeStarted();


            if (temp_process.getCpuBursts().size() == (unsigned int)cpu_burst_index[temp_process.getLabel()]) {
                /* PROCESS TERMINATES */
                std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " terminated [Q";
                if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                else {
                    for (std::vector<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                        std::cout << " " << it->getLabel();
                    }
                    std::cout << "]" << std::endl;
                }
                counter++;
                //terminated.push_back(temp_process);
                time += (context_switch / 2);

            } else {
                if (time < 10000) {
                    std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " (tau " << (int)ceil(temp_process.getTau()) << "ms) completed a CPU burst; "
                            << temp_process.getCpuBursts().size() - cpu_burst_index[temp_process.getLabel()] << " bursts to go [Q";
                    if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                    else 
                    {
                        for (std::vector<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                            std::cout << " " << it->getLabel();
                        }
                        std::cout << "]" << std::endl;
                    }
                }
                
                double old_tau = temp_process.getTau();
                
                temp_process.setTau(ceil(alpha * temp_process.getCpuBurst(cpu_burst_index[temp_process.getLabel()] - 1) + ( 1.0 - alpha ) * (int)old_tau));
                
                //Setting time ran to zero
                temp_process.updateTimeRan(-1 * temp_process.getTimeRan());

                if (time < 10000) { 
                    std::cout << "time " << time <<  "ms: Recalculating tau for process " << temp_process.getLabel() << ": old tau " << (int)old_tau <<
                              "ms ==> new tau " << (int)ceil(temp_process.getTau()) << "ms [Q"; 
                    if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                    else {
                        for (std::vector<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                            std::cout << " " << it->getLabel();
                        }
                        std::cout << "]" << std::endl;
                    }
                }


                int temp_burst_len = temp_process.getIOBurst(io_burst_index[temp_process.getLabel()]);
                temp_process.setDoneWaiting(time + temp_burst_len + (context_switch / 2));
                if (time < 10000) {
                    std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " switching out of CPU; blocking on I/O until time "
                            << time + temp_burst_len + (context_switch / 2) << "ms [Q";
                    if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                    else {
                        for (std::vector<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                            std::cout << " " << it->getLabel();
                        }
                        std::cout << "]" << std::endl;
                    }
                }
                io_burst_index[temp_process.getLabel()] += 1;
                waiting.push_back(temp_process);
                std::sort(waiting.begin(), waiting.end(), fcfs_compare_wait);
                time += (context_switch / 2);
            }
        } else if (done_blocking_time != INT_MAX) {
            /* DONE WITH WAITING STATE */
            /* PREEMPTING MAY OCCUR */

                           


            /* take the process that is exiting the waiting queue and store it in temp_process */
            Process temp_process = waiting.front();
            waiting.erase(waiting.begin());

            /* update the time to when the new process leaves the waiting queue */
            time = temp_process.getDoneWaiting();
            temp_process.setDoneWaiting(INT_MAX);
            temp_process.setTimeStartedWaiting(time);


            if (running.size() == 0) {
                /* if nothing is in the running queue, add the process to the ready queue and sort it */
                ready.push_back(temp_process);
                std::sort(ready.begin(), ready.end(), srt_compare_time_remaining);
                
                if (time < 10000) {
                    std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " (tau " << (int)ceil(temp_process.getTau()) << "ms) completed I/O; added to ready queue [Q";
                    if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                    else {
                        for (std::vector<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                            std::cout << " " << it->getLabel();
                        }
                        std::cout << "]" << std::endl;
                    }
                }

            } else {
                /* there is a process in the running queue, so check if we need to preempt it */
                
                Process preempted_process = running.front();
                /* update the time ran for the preempted_process */
                preempted_process.updateTimeRan(time - preempted_process.getTimeStarted());


                
                if (preempted_process.getTau() - preempted_process.getTimeRan() > temp_process.getTau() - temp_process.getTimeRan()) {
                    /* preempt the current process in the running queue and replace it with the new process */
                    temp_process.setTimeStarted(time);
                    
                    timeCpuRunning += time - preempted_process.getTimeStarted();

                    ready.push_back(temp_process);
                    std::sort(ready.begin(), ready.end(), srt_compare_time_remaining);

                    if (time < 10000) {
                        std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " (tau " << (int)ceil(temp_process.getTau()) << "ms) completed I/O; preempting " << preempted_process.getLabel() << " [Q";
                        if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                        else {
                            for (std::vector<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                                std::cout << " " << it->getLabel();
                            }
                            std::cout << "]" << std::endl;
                        }
                    }
                    
                    running.erase(running.begin());
                    
                    temp_process.setTimeStarted(time);
                    preempted_process.setTimeStartedWaiting(time);
                    ready.push_back(preempted_process);
                    std::sort(ready.begin(), ready.end(), srt_compare_time_remaining);
                    cpu_burst_index[preempted_process.getLabel()] -= 1;
                    
                    if (preempted_process.getCpuBound()) {
                        num_preemptions_CPU++;
                    } else {
                        num_preemptions_IO++;
                    }
    
                    time += (context_switch / 2);
                    
                } else {
                    /* no preemption necessary */
                    temp_process.setTimeStartedWaiting(time);
                    ready.push_back(temp_process);
                    std::sort(ready.begin(), ready.end(), srt_compare_time_remaining);

                    if (time < 10000) {
                        std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " (tau " << (int)ceil(temp_process.getTau()) << "ms) completed I/O; added to ready queue [Q";
                        if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                        else {
                            for (std::vector<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                                std::cout << " " << it->getLabel();
                            }
                            std::cout << "]" << std::endl;
                        }
                    }
                }
            }
        }

        if (running.size() == 0 && ready.size() > 0) {
            /* MOVE NEXT PROCESS FROM READY STATE TO RUNNING STATE */
            Process temp_process = ready.front();
            ready.erase(ready.begin());

            /* fully completed last cpu burst */
            if (temp_process.getTimeRan() == 0) {
                if (temp_process.getCpuBound()) {
                    sumWaitTime_CPU += time - temp_process.getTimeStartedWaiting();
                } else {
                    sumWaitTime_IO += time - temp_process.getTimeStartedWaiting();
                }
                time += (context_switch / 2);
                
                int temp_burst_len = temp_process.getCpuBurst(cpu_burst_index[temp_process.getLabel()]);
                temp_process.setDoneRunning(time + temp_burst_len);
                temp_process.setTimeStarted(time);
                if (time < 10000) {
                    std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " (tau " << (int)ceil(temp_process.getTau()) << "ms) started using the CPU for "
                            << temp_burst_len << "ms burst [Q";
                    if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                    else {
                        for (std::vector<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                            std::cout << " " << it->getLabel();
                        }
                        std::cout << "]" << std::endl;
                    }
                }

                cpu_burst_index[temp_process.getLabel()] += 1;

                running.push_back(temp_process);
            } else {
                /* finish previous cpu burst */

                if (temp_process.getCpuBound()) {
                    sumWaitTime_CPU += time - temp_process.getTimeStartedWaiting() - (context_switch / 2);
                    // sumCpuTime_CPU += context_switch;
                    temp_total_cpu_time += context_switch;
                } else {
                    sumWaitTime_IO += time - temp_process.getTimeStartedWaiting();
                }

                time += (context_switch / 2);
                int temp_burst_len = temp_process.getCpuBurst(cpu_burst_index[temp_process.getLabel()]);
                temp_process.setDoneRunning(time + temp_burst_len - temp_process.getTimeRan());
                temp_process.setTimeStarted(time);
                if (time < 10000) {
                    std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " (tau " << (int)ceil(temp_process.getTau()) << "ms) started using the CPU for remaining "
                            << temp_burst_len - temp_process.getTimeRan() << "ms of " << temp_burst_len << "ms burst [Q";
                    if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                    else {
                        for (std::vector<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                            std::cout << " " << it->getLabel();
                        }
                        std::cout << "]" << std::endl;
                    }
                }

                cpu_burst_index[temp_process.getLabel()] += 1;
                running.push_back(temp_process);
            }
        }
    }
    std::cout << "time " << time << "ms: Simulator ended for SRT [Q <empty>]" << std::endl;

    double cpu_utilization = 0;
    if(time != 0){
        cpu_utilization = ceil(1000 * ((double)timeCpuRunning / time) * 100) / 1000;
    }
    
    double cpuBurstAvg_CPU = 0;
    if(numCpuBursts_CPU != 0){
        cpuBurstAvg_CPU = ceil(1000 * (double)sumCpuTime_CPU / numCpuBursts_CPU) / 1000;
    }

    double cpuBurstAvg_IO = 0;
    if(numCpuBursts_IO != 0){
        cpuBurstAvg_IO = ceil(1000 * (double)sumCpuTime_IO / numCpuBursts_IO) / 1000;
    }

    double cpuBurstAvg_SUM = 0;
    if(numCpuBursts_CPU + numCpuBursts_IO != 0){
        cpuBurstAvg_SUM = ceil(1000.0 * (double)(sumCpuTime_CPU + sumCpuTime_IO) / (numCpuBursts_CPU + numCpuBursts_IO)) / 1000;
    }

    double IOBurstAvg_CPU = 0;
    if(numCpuBursts_CPU != 0){
        IOBurstAvg_CPU = ceil(1000 * (double)sumWaitTime_CPU / numCpuBursts_CPU) / 1000;
    }

    double IOBurstAvg_IO = 0;
    if(numCpuBursts_IO != 0){
        IOBurstAvg_IO = ceil(1000 * (double)sumWaitTime_IO / numCpuBursts_IO) / 1000;
    }

    double IOBurstAvg_SUM = 0;
    if(numCpuBursts_CPU + numCpuBursts_IO != 0){
        IOBurstAvg_SUM = ceil(1000.0 * (double)(sumWaitTime_CPU + sumWaitTime_IO) / (numCpuBursts_CPU + numCpuBursts_IO)) / 1000;
    }

    double turnaround = 0;
    if (numCpuBursts_CPU + numCpuBursts_IO != 0) {
        turnaround = (double)(sumWaitTime_CPU + sumWaitTime_IO) / (numCpuBursts_CPU + numCpuBursts_IO) 
                        + (double)(temp_total_cpu_time + sumCpuTime_IO) / (numCpuBursts_CPU + numCpuBursts_IO) + context_switch;
    }

    double temp_avg = 0;
    if (numCpuBursts_CPU != 0) {
        temp_avg = ceil(1000 * (double)temp_total_cpu_time / numCpuBursts_CPU) / 1000;
    }
    
    output << "Algorithm SRT\n";
    output << "-- CPU utilization: " << std::fixed << std::setprecision(3) << cpu_utilization << "%\n";
    output << "-- average CPU burst time: " << cpuBurstAvg_SUM << " ms (" << cpuBurstAvg_CPU << " ms/" << cpuBurstAvg_IO << " ms)\n";
    output << "-- average wait time: " << IOBurstAvg_SUM << " ms (" << IOBurstAvg_CPU << " ms/" << IOBurstAvg_IO << " ms)\n";
    output << "-- average turnaround time: " << ceil(1000 * turnaround) / 1000 << " ms (" << temp_avg + IOBurstAvg_CPU + context_switch 
           << " ms/" << cpuBurstAvg_IO + IOBurstAvg_IO + context_switch << " ms)\n";
    output << "-- number of context switches: " << numCpuBursts_CPU + numCpuBursts_IO + num_preemptions_IO + num_preemptions_CPU << " (" << numCpuBursts_CPU + num_preemptions_CPU << "/" << numCpuBursts_IO + num_preemptions_IO << ")\n";
    output << "-- number of preemptions: " << num_preemptions_IO + num_preemptions_CPU << " (" << num_preemptions_CPU << "/" << num_preemptions_IO << ")\n\n";
}

void rr(std::vector<Process> processes, int context_switch, int time_slice, std::ofstream &output, clock_t start) {

    int timeCpuRunning = 0; 
    int sumCpuTime_CPU = 0;
    int numCpuBursts_CPU = 0;
    int sumCpuTime_IO = 0;
    int numCpuBursts_IO = 0;
    int sumWaitTime_CPU = 0;
    int sumWaitTime_IO = 0;
    int num_preemptions_CPU = 0;
    int num_preemptions_IO = 0;
    clock_t end;

    int time = 0;
    std::cout << "time " << time << "ms: Simulator started for RR [Q <empty>]" << std::endl;

    std::sort(processes.begin(), processes.end(), fcfs_compare_arr);
    int next_arrival = 0;
    int magic_number = 1;

    std::list<Process> running;
    std::vector<Process> waiting;
    std::list<Process> ready;
    std::map<char, int> cpu_burst_index;
    std::map<char, int> io_burst_index;
    std::vector<Process> terminated;
    

    for (unsigned int i = 0; i < processes.size(); i++) {
        cpu_burst_index[processes[i].getLabel()] = 0;
        io_burst_index[processes[i].getLabel()] = 0;

        std::vector<int> temp_cpu_bursts = processes[i].getCpuBursts();
        std::vector<int> temp_io_bursts = processes[i].getIoBursts();
        
        if (processes[i].getCpuBound()) { 
            for (unsigned int j = 0; j < temp_cpu_bursts.size(); j++) { sumCpuTime_CPU += temp_cpu_bursts[j]; }
            numCpuBursts_CPU += temp_cpu_bursts.size();
        } else {
            for (unsigned int j = 0; j < temp_cpu_bursts.size(); j++) { sumCpuTime_IO += temp_cpu_bursts[j]; }
            numCpuBursts_IO += temp_cpu_bursts.size();
        }
    }

    int temp_total_cpu_time = sumCpuTime_CPU;
    int temp_total_cpu_time_io = sumCpuTime_IO;

    while (processes.size() > terminated.size()) {
        end = clock();
        if (end - start > 100000){
            break;
        }
        
         //Remove process from processes once finished, not arrived
        // Figure out what happens next: arrival, finishing waiting state, or finishing running state

        // next_arrival is the index of the next process to arrive from processes

        /* the time when the next process arrives */
        int next_arrival_time = INT_MAX;
        if ((unsigned int)next_arrival < processes.size()) {
            next_arrival_time = processes[next_arrival].getArrivalTime();
        }

        /* the time when the burst using the CPU finishes and is ready to leave the running queue */
        int done_running_time = INT_MAX;
        if (running.size() > 0) {
            done_running_time = running.front().getDoneRunning();
        }

        /* the time when the burst blocking on IO finishes and is ready to leave the waiting queue */
        int done_blocking_time = INT_MAX;
        if (waiting.size() > 0) {
            done_blocking_time = waiting.front().getDoneWaiting();
        }
        
        int next_time_slice = INT_MAX;
        if (running.size() > 0) {
            Process preempted_process = running.front();
            next_time_slice = (magic_number * time_slice) + preempted_process.getTimeStarted();
        }

        
        

        if (next_time_slice < next_arrival_time && next_time_slice < done_running_time && next_time_slice < done_blocking_time) {
            /* TIME SLICE EXPIRED */

            if (ready.size() > 0) {
                /* preempt the current running process and send it to the back of the ready queue */
                Process preempted_process = running.front();
                preempted_process.updateTimeRan(magic_number * time_slice);
                running.erase(running.begin());
                
                int temp_burst_len = preempted_process.getCpuBurst(cpu_burst_index[preempted_process.getLabel()] - 1);
                time = preempted_process.getTimeStarted() + (magic_number * time_slice);
                timeCpuRunning += time - preempted_process.getTimeStarted();
                                    
                if (time < 10000) {
                    std::cout << "time " << time << "ms: Time slice expired; preempting process " << preempted_process.getLabel() << " with " << temp_burst_len - preempted_process.getTimeRan() << "ms remaining [Q"; 
                    if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                    else {
                        for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                            std::cout << " " << it->getLabel();
                        }
                        std::cout << "]" << std::endl;
                    }
                }
                
                time += (context_switch / 2);
                preempted_process.setTimeStartedWaiting(time);
                cpu_burst_index[preempted_process.getLabel()] -= 1;
                ready.push_back(preempted_process);
                magic_number = 1;
                
                if (preempted_process.getCpuBound()) {
                    num_preemptions_CPU++;
                } else {
                    num_preemptions_IO++;
                }

            } else {
                time += time_slice;
                Process running_process = running.front();
                running_process.updateTimeRan(time_slice);
                magic_number++;
                
                if (time < 10000) {
                    std::cout << "time " << time << "ms: Time slice expired; no preemption because ready queue is empty [Q <empty>]" << std::endl;
                }
            }

        } else if (next_arrival_time <= done_running_time && next_arrival_time <= done_blocking_time && next_arrival_time != INT_MAX) {
            /* NEW ARRIVAL */
            Process temp_process = processes[next_arrival];
            time = temp_process.getArrivalTime();
            temp_process.setTimeStartedWaiting(time);
            ready.push_back(temp_process);
            next_arrival += 1;
            if (time < 10000) {
                std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " arrived; added to ready queue [Q";
                if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                else {
                    for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                        std::cout << " " << it->getLabel();
                    }
                    std::cout << "]" << std::endl;
                }
            }
            
        } else if (done_running_time < next_arrival_time && done_running_time <= done_blocking_time && running.size() > 0) {
            /* DONE WITH RUNNING STATE */
            Process temp_process = running.front();
            running.pop_front();

            time = temp_process.getDoneRunning();
            temp_process.setDoneRunning(INT_MAX);
            temp_process.updateTimeRan(-1 * temp_process.getTimeRan());
            timeCpuRunning += time - temp_process.getTimeStarted();
            magic_number = 1;

            if (temp_process.getCpuBursts().size() == (unsigned int)cpu_burst_index[temp_process.getLabel()]) {
                /* PROCESS TERMINATED */
                std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " terminated [Q";
                if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                else {
                    for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                        std::cout << " " << it->getLabel();
                    }
                    std::cout << "]" << std::endl;
                }
                terminated.push_back(temp_process);

                time += (context_switch / 2);

            } else {

                if (time < 10000) {
                    std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " completed a CPU burst; "
                            << temp_process.getCpuBursts().size() - cpu_burst_index[temp_process.getLabel()]
                            << " bursts to go [Q";
                    if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                    else {
                        for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                            std::cout << " " << it->getLabel();
                        }
                        std::cout << "]" << std::endl;   
                    }
                }


                int temp_burst_len = temp_process.getIOBurst(io_burst_index[temp_process.getLabel()]);
                temp_process.setDoneWaiting(time + temp_burst_len + (context_switch / 2));
                
                if (time < 10000) {
                    std::cout << "time " << time << "ms: Process " << temp_process.getLabel()
                          << " switching out of CPU; blocking on I/O until time "
                          << time + temp_burst_len + (context_switch / 2) << "ms [Q";
                    if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                    else {
                        for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                            std::cout << " " << it->getLabel();
                        }
                        std::cout << "]" << std::endl;
                    }
                }
                io_burst_index[temp_process.getLabel()] += 1;
                waiting.push_back(temp_process);
                std::sort(waiting.begin(), waiting.end(), fcfs_compare_wait);
                time += (context_switch / 2);
            }
        } else {
            /* DONE WITH WAITING STATE */
            Process temp_process = waiting.front();
            waiting.erase(waiting.begin());

            time = temp_process.getDoneWaiting();
            temp_process.setDoneWaiting(INT_MAX);
            temp_process.setTimeStartedWaiting(time);
            ready.push_back(temp_process);

            if (time < 10000) {
                std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " completed I/O; added to ready queue [Q";
                if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                else {
                    for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                        std::cout << " " << it->getLabel();
                    }
                    std::cout << "]" << std::endl;
                }
            }
            
        }

        if (running.size() == 0 && ready.size() > 0) {
            /* MOVE NEXT PROCESS FROM READY STATE TO RUNNING STATE */
            Process temp_process = ready.front();
            ready.pop_front();

            if (temp_process.getCpuBound()) {
                sumWaitTime_CPU += time - temp_process.getTimeStartedWaiting();
            } else {
                sumWaitTime_IO += time - temp_process.getTimeStartedWaiting();
            }

            time += (context_switch / 2);

            if (temp_process.getTimeRan() == 0) {
                /* the process never got cut off */
            
                int temp_burst_len = temp_process.getCpuBurst(cpu_burst_index[temp_process.getLabel()]);
                temp_process.setDoneRunning(time + temp_burst_len);
                temp_process.setTimeStarted(time);
                if (time < 10000) {
                    std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " started using the CPU for "
                            << temp_burst_len << "ms burst [Q";
                    if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                    else {
                        for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                            std::cout << " " << it->getLabel();
                        }
                        std::cout << "]" << std::endl;
                    }
                }

                cpu_burst_index[temp_process.getLabel()] += 1;
                
                running.push_back(temp_process);

            } else {
                /* finish the cpu burst that already started */

                if (temp_process.getCpuBound()) {
                    temp_total_cpu_time += context_switch;
                } else {
                    temp_total_cpu_time_io += context_switch;
                }

                int temp_burst_len = temp_process.getCpuBurst(cpu_burst_index[temp_process.getLabel()]);
                temp_process.setDoneRunning(time + temp_burst_len - temp_process.getTimeRan());
                temp_process.setTimeStarted(time);
                if (time < 10000) {
                    std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " started using the CPU for remaining "
                            << temp_burst_len - temp_process.getTimeRan() << "ms of " << temp_burst_len << "ms burst [Q";
                    if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                    else {
                        for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                            std::cout << " " << it->getLabel();
                        }
                        std::cout << "]" << std::endl;
                    }
                }

                cpu_burst_index[temp_process.getLabel()] += 1;
                running.push_back(temp_process);
            }
        }
    }

    std::cout << "time " << time << "ms: Simulator ended for RR [Q <empty>]" << std::endl;

    double cpu_utilization = 0;
    if(time != 0){
        cpu_utilization = ceil(1000 * ((double)timeCpuRunning / time) * 100) / 1000;
    }
    
    double cpuBurstAvg_CPU = 0;
    if(numCpuBursts_CPU != 0){
        cpuBurstAvg_CPU = ceil(1000 * (double)sumCpuTime_CPU / numCpuBursts_CPU) / 1000;
    }

    double cpuBurstAvg_IO = 0;
    if(numCpuBursts_IO != 0){
        cpuBurstAvg_IO = ceil(1000 * (double)sumCpuTime_IO / numCpuBursts_IO) / 1000;
    }

    double cpuBurstAvg_SUM = 0;
    if(numCpuBursts_CPU + numCpuBursts_IO != 0){
        cpuBurstAvg_SUM = ceil(1000.0 * (double)(sumCpuTime_CPU + sumCpuTime_IO) / (numCpuBursts_CPU + numCpuBursts_IO)) / 1000;
    }

    double IOBurstAvg_CPU = 0;
    if(numCpuBursts_CPU != 0){
        IOBurstAvg_CPU = ceil(1000 * (double)sumWaitTime_CPU / numCpuBursts_CPU) / 1000;
    }

    double IOBurstAvg_IO = 0;
    if(numCpuBursts_IO != 0){
        IOBurstAvg_IO = ceil(1000 * (double)sumWaitTime_IO / numCpuBursts_IO) / 1000;
    }

    double IOBurstAvg_SUM = 0;
    if(numCpuBursts_CPU + numCpuBursts_IO != 0){
        IOBurstAvg_SUM = ceil(1000.0 * (double)(sumWaitTime_CPU + sumWaitTime_IO) / (numCpuBursts_CPU + numCpuBursts_IO)) / 1000;
    }

    double turnaround = 0;
    if (numCpuBursts_CPU + numCpuBursts_IO != 0) {
        turnaround = (double)(sumWaitTime_CPU + sumWaitTime_IO) / (numCpuBursts_CPU + numCpuBursts_IO) + (double)(temp_total_cpu_time + temp_total_cpu_time_io) / (numCpuBursts_CPU + numCpuBursts_IO) + context_switch;
    }

    double temp_avg_cpu = 0;
    if (numCpuBursts_CPU != 0) {
        temp_avg_cpu = ceil(1000 * (double)temp_total_cpu_time / numCpuBursts_CPU) / 1000;
    }

    double temp_avg_io = 0;
    if (numCpuBursts_IO != 0) {
        temp_avg_io = ceil(1000 * (double)temp_total_cpu_time_io / numCpuBursts_IO) / 1000;
    }
    
    output << "Algorithm RR\n";
    output << "-- CPU utilization: " << std::fixed << std::setprecision(3) << cpu_utilization << "%\n";
    output << "-- average CPU burst time: " << cpuBurstAvg_SUM << " ms (" << cpuBurstAvg_CPU << " ms/" << cpuBurstAvg_IO << " ms)\n";
    output << "-- average wait time: " << IOBurstAvg_SUM << " ms (" << IOBurstAvg_CPU << " ms/" << IOBurstAvg_IO << " ms)\n";
    output << "-- average turnaround time: " << ceil(1000 * turnaround) / 1000 << " ms (" << temp_avg_cpu + IOBurstAvg_CPU + context_switch 
           << " ms/" << temp_avg_io + IOBurstAvg_IO + context_switch << " ms)\n";
    output << "-- number of context switches: " << numCpuBursts_CPU + numCpuBursts_IO + num_preemptions_IO + num_preemptions_CPU << " (" << numCpuBursts_CPU + num_preemptions_CPU << "/" << numCpuBursts_IO + num_preemptions_IO << ")\n";
    output << "-- number of preemptions: " << num_preemptions_IO + num_preemptions_CPU << " (" << num_preemptions_CPU << "/" << num_preemptions_IO << ")\n";

}

int main(int argc, char** argv) {
    if (argc != 9) {
        std::cerr << char(91) << "Usage" << char(93) << " " << *(argv) << " <n> <n_cpu> <seed> <lambda> <threshold>"
                  << std::endl;
        return EXIT_FAILURE;
    }

   
    int n = atoi(*(argv + 1));
    int n_cpu = atoi(*(argv + 2));
    long seed = atol(*(argv + 3));
    double lambda = atof(*(argv + 4));
    int threshold = atoi(*(argv + 5));
    int context_switch = atoi(*(argv + 6));
    double alpha = atof(*(argv + 7));
    int time_slice = atoi(*(argv + 8));


    if (n > 26) {
        std::cerr << "[Usage] too many processes were specified (max: 26)" << std::endl;
        return EXIT_FAILURE;
    }
    int n_io = n - n_cpu;
    if (n_io < 0) {
        std::cerr << "[Usage] there can't be more cpu processes than total processes"
                  << std::endl;
        return EXIT_FAILURE;
    }
    if (context_switch < 0 || context_switch % 2 == 1) {
        std::cerr << "[Usage] context switch must be a positive even number" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "<<< PROJECT PART I -- process set (n=" << n << ") with " << n_cpu;
    if (n_cpu == 1) { std::cout << " CPU-bound process >>>" << std::endl; }
    else { std::cout << " CPU-bound processes >>>" << std::endl; }

    std::vector<Process> processes;

    srand48(seed);
    for (int i = 0; i < n; i++) {
        int arrival_time = floor(next_exp(lambda, threshold));
        int num_bursts = ceil(drand48() * 64);
        if (n - i <= n_cpu) { // CPU-bound process
            std::cout << "CPU-bound process " << char(65 + i) << ": arrival time " << arrival_time << "ms; "
                      << num_bursts << " CPU bursts" << std::endl;
        } else { // IO-bound process
            std::cout << "I/O-bound process " << char(65 + i) << ": arrival time " << arrival_time << "ms; "
                      << num_bursts << " CPU bursts" << std::endl;
        }

        std::vector<int> cpu_bursts;
        std::vector<int> io_bursts;
        for (int j = 0; j < num_bursts; j++) {
            int cpu_time = ceil(next_exp(lambda, threshold));
            if (n - i <= n_cpu) { cpu_time *= 4; } // CPU-bound process
            cpu_bursts.push_back(cpu_time);
            if (j < num_bursts - 1) {
                int io_time = ceil(next_exp(lambda, threshold)) * 10;
                if (n - i <= n_cpu) { io_time /= 8; } // CPU-bound process
                io_bursts.push_back(io_time);
            }
        }
        processes.push_back(Process(char(65 + i), n - i <= n_cpu, arrival_time, num_bursts, cpu_bursts, io_bursts));
    }

    std::cout << std::endl;
    std::cout << "<<< PROJECT PART II -- t_cs=" << context_switch << "ms; alpha=" << std::fixed << std::setprecision(2) << alpha << "; t_slice=" << time_slice
              << "ms >>>" << std::endl;


    std::ofstream output;
    output.open("simout.txt");
    

    fcfs(processes, context_switch, output, start);
    std::cout << std::endl;
    sjf(processes, context_switch, lambda, alpha, output, start);
    std::cout << std::endl;
    srt(processes, context_switch, lambda, alpha, output, start);
    std::cout << std::endl;
    rr(processes, context_switch, time_slice, output, start);

    output.close();


    return EXIT_SUCCESS;
}