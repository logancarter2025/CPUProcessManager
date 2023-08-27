#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <math.h>
#include <fstream>
#include <algorithm>
#include <time.h>
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
    int timeCpuRunning = 0; //line 105
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
    std::list<Process> waiting;
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
        if ((double)(end - start)/CLOCKS_PER_SEC > 0.25){
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
                waiting.sort(fcfs_compare_wait);
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
    //time += (context_switch / 2);
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
    std::list<Process> waiting;
    std::list<Process> ready;
    
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
        if ((double)(end - start)/CLOCKS_PER_SEC > 0.25){
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
            ready.sort(sjf_compare_tau);
            next_arrival += 1;

            if (time < 10000) { 
                std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " (tau " << (int)ceil(temp_process.getTau()) << "ms) arrived; added to ready queue [Q"; 
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
                
                time += (context_switch / 2);

                terminated.push_back(temp_process);
            } else {
                if (time < 10000) {
                    std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " (tau " << (int)ceil(temp_process.getTau()) << "ms) completed a CPU burst; "
                            << temp_process.getCpuBursts().size() - cpu_burst_index[temp_process.getLabel()] << " bursts to go [Q";
                    if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                    else {
                        for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                            std::cout << " " << it->getLabel();
                        }
                        std::cout << "]" << std::endl;
                    }
                }

                double old_tau = temp_process.getTau();
                /* update tau for temp_process */
                // cpu->tau = ceil( alpha * cpu->burst + ( 1.0 - alpha ) * cpu->tau );
                
                temp_process.setTau(ceil(alpha * temp_process.getCpuBurst(cpu_burst_index[temp_process.getLabel()] - 1) + ( 1.0 - alpha ) * (int)old_tau));

                if (time < 10000) { 
                    std::cout << "time " << time <<  "ms: Recalculating tau for process " << temp_process.getLabel() << ": old tau " << (int)old_tau <<
                              "ms ==> new tau " << (int)ceil(temp_process.getTau()) << "ms [Q"; 
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


                if (temp_process.getCpuBound()) {
                    cpuCS++;
                } else {
                    ioCS++;
                }

                io_burst_index[temp_process.getLabel()] += 1;
                waiting.push_back(temp_process);
                waiting.sort(fcfs_compare_wait);
                time += (context_switch / 2);

            }
            
        } else {
            /* DONE WITH WAITING STATE */
            Process temp_process = waiting.front();
            waiting.pop_front();

            time = temp_process.getDoneWaiting();
            temp_process.setDoneWaiting(INT_MAX);
            temp_process.setTimeStartedWaiting(time);
            
            ready.push_back(temp_process);
            ready.sort(sjf_compare_tau);

            if (time < 10000) {
                std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " (tau " << (int)ceil(temp_process.getTau()) << "ms) completed I/O; added to ready queue [Q";
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
            if (time < 10000) {
                std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " (tau " << (int)ceil(temp_process.getTau()) << "ms) started using the CPU for "
                        << temp_burst_len << "ms burst [Q";
                if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                else {
                    for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
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
    std::list<Process> waiting;
    std::list<Process> ready;
    
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
        if ((double)(end - start)/CLOCKS_PER_SEC > 0.25){
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
            ready.sort(srt_compare_time_remaining);
            next_arrival += 1;

            if (time < 10000) { 
                std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " (tau " << (int)ceil(temp_process.getTau()) << "ms) arrived; added to ready queue [Q"; 
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
                    for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
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
                        for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                            std::cout << " " << it->getLabel();
                        }
                        std::cout << "]" << std::endl;
                    }
                }
                
                double old_tau = temp_process.getTau();
                /* update tau for temp_process */
                // cpu->tau = ceil( alpha * cpu->burst + ( 1.0 - alpha ) * cpu->tau );
                
                temp_process.setTau(ceil(alpha * temp_process.getCpuBurst(cpu_burst_index[temp_process.getLabel()] - 1) + ( 1.0 - alpha ) * (int)old_tau));

                //sets to zero
                
                temp_process.updateTimeRan(-1 * temp_process.getTimeRan());

                if (time < 10000) { 
                    std::cout << "time " << time <<  "ms: Recalculating tau for process " << temp_process.getLabel() << ": old tau " << (int)old_tau <<
                              "ms ==> new tau " << (int)ceil(temp_process.getTau()) << "ms [Q"; 
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
                    std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " switching out of CPU; blocking on I/O until time "
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
                waiting.sort(fcfs_compare_wait);
                time += (context_switch / 2);
            }
        } else if (done_blocking_time != INT_MAX) {
            /* DONE WITH WAITING STATE */
            /* PREEMPTING MAY OCCUR */

            /* steps: 
             *  
             *  1: check if anything is in the running queue 
             *    2: if so, check if the process leaving the waiting queue will finish its cpu burst first
             *      3: if so, preempt the current process in the running queue and replace it with the new process
             *         -- the tricky part here is updating the times correctly
             *      4: if the new process will take longer than the currently running process, put the new one in the ready queue and sort
             *  5: if nothing is in the running queue, add the process to the ready queue and sort it
             *

             */
                           


            /* take the process that is exiting the waiting queue and store it in temp_process */
            Process temp_process = waiting.front();
            waiting.pop_front();

            /* update the time to when the new process leaves the waiting queue */
            time = temp_process.getDoneWaiting();
            temp_process.setDoneWaiting(INT_MAX);
            temp_process.setTimeStartedWaiting(time);


            if (running.size() == 0) {
                /* if nothing is in the running queue, add the process to the ready queue and sort it */
                ready.push_back(temp_process);
                ready.sort(srt_compare_time_remaining);
                
                if (time < 10000) {
                    std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " (tau " << (int)ceil(temp_process.getTau()) << "ms) completed I/O; added to ready queue [Q";
                    if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                    else {
                        for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
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
                    ready.sort(srt_compare_time_remaining);

                    if (time < 10000) {
                        std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " (tau " << (int)ceil(temp_process.getTau()) << "ms) completed I/O; preempting " << preempted_process.getLabel() << " [Q";
                        if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                        else {
                            for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                                std::cout << " " << it->getLabel();
                            }
                            std::cout << "]" << std::endl;
                        }
                    }
                    
                    running.pop_front();
                    
                    temp_process.setTimeStarted(time);
                    preempted_process.setTimeStartedWaiting(time);
                    ready.push_back(preempted_process);
                    ready.sort(srt_compare_time_remaining);
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
                    ready.sort(srt_compare_time_remaining);

                    if (time < 10000) {
                        std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " (tau " << (int)ceil(temp_process.getTau()) << "ms) completed I/O; added to ready queue [Q";
                        if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                        else {
                            for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
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
            ready.pop_front();

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
                        for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
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


void rr(std::vector<Process> processes, int context_switch, int time_slice, std::ofstream &output) {
    
    /* FOR ANALYTICS */
    int num_cpu_bursts = 0;
    int num_io_bursts = 0;
    int sum_cpu_time = 0;
    int sum_io_time = 0;
    int sumWaitTime_CPU = 0;
    int sumWaitTime_IO = 0;

    int time = 0;
    std::cout << "time " << time << "ms: Simulator started for RR [Q <empty>]" << std::endl;

    std::sort(processes.begin(), processes.end(), fcfs_compare_arr);
    
    /* initialize three queues: running is for processes actively doing a CPU burst, waiting is for
       processes actively blocking on IO (doing an IO burst), and ready is for processes in line to 
       do a CPU burst */
    std::list<Process> running;
    std::list<Process> waiting;
    std::list<Process> ready;
    unsigned int num_terminated = 0;
    unsigned int next_arrival = 0;
    int magic_number = 1;

    /* initialize two maps that keep track of the respective indices each Process is on */
    std::map<char, int> cpu_burst_index;
    std::map<char, int> io_burst_index;

    /* set initial cpu burst and io burst indices to 0 and calculate the total cpu time and number of bursts */
    for (unsigned int i = 0; i < processes.size(); i++) {
        cpu_burst_index[processes[i].getLabel()] = 0;
        io_burst_index[processes[i].getLabel()] = 0;
        if (processes[i].getCpuBound()) {
            for (unsigned int j = 0; j < processes[i].getCpuBursts().size(); j++) {
                sum_cpu_time += processes[i].getCpuBursts()[j];
            }
            num_cpu_bursts += processes[i].getCpuBursts().size();
        } else {
            for (unsigned int j = 0; j < processes[i].getIoBursts().size(); j++) {
                sum_io_time += processes[i].getIoBursts()[j];
            }
            num_io_bursts += processes[i].getIoBursts().size();
        }
    }

    /* calculate which event occurs next: 
         1. CPU burst completion/time_slice expiring            (running -> waiting/terminated++)
         2. process starts using the CPU    (ready -> running)
         3. I/O burst completion            (waiting -> ready)
         4. New arrival                     (add process to ready) */

    while (processes.size() > num_terminated) {
        bool ts_expired = false;
        int cpu_burst_completion = INT_MAX;
        if (running.size() > 0) {
            Process running_process = running.front();
            int burst_completion = running_process.getCpuBurst(cpu_burst_index[running_process.getLabel()]) + running_process.getTimeStarted();
            int time_slice_expire = running_process.getTimeStarted() + (magic_number * time_slice);
            cpu_burst_completion = std::min(burst_completion, time_slice_expire);
            if (cpu_burst_completion == time_slice_expire) { ts_expired = true; }
        }

        int cpu_burst_starts = INT_MAX;
        if (running.size() == 0 && ready.size() > 0) { 
            cpu_burst_starts = time + (context_switch / 2); 
        }
        
        int io_burst_completion = INT_MAX;
        if (waiting.size() > 0) {
            Process waiting_process = waiting.front();
            io_burst_completion = waiting_process.getIOBurst(io_burst_index[waiting_process.getLabel()]) + time;
        }
        
        int next_arrival_time = INT_MAX;
        if (next_arrival < processes.size()){
            next_arrival_time = processes[next_arrival].getArrivalTime();
        }

        
        
        if (cpu_burst_completion <= cpu_burst_starts && cpu_burst_completion <= io_burst_completion && cpu_burst_completion <= next_arrival_time) {
            /* CPU BURST COMPLETION */
            
            time = cpu_burst_completion;
            if (ts_expired) {
                /* TIME SLICE EXPIRED */
                Process preempted_process = running.front();
                preempted_process.updateTimeRan(time_slice);

                if (ready.size() > 0) {
                    /* PREEMPTION OCCURS */
                    running.pop_front();

                    int temp_burst_len = preempted_process.getCpuBurst(cpu_burst_index[preempted_process.getLabel()]);
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

                } else {
                    //Nothing in ready queue
                    std::cout << "  process: " << preempted_process.getLabel() << ", time ran: " << preempted_process.getTimeRan() << ", est len: " << preempted_process.getCpuBurst(cpu_burst_index[preempted_process.getLabel()]) << std::endl;
                    magic_number++;
                    if (time < 10000) {
                        std::cout << "time " << time << "ms: Time slice expired; no preemption because ready queue is empty [Q <empty>]" << std::endl;
                    }
                    
                }

            } else {
                /* CPU BURST COMPLETED */
                Process finished_process = running.front();
                running.pop_front();
                finished_process.updateTimeRan(-1 * finished_process.getTimeRan());
                
                cpu_burst_index[finished_process.getLabel()]++;
                if (cpu_burst_index[finished_process.getLabel()] == finished_process.getCpuBursts().size()) {
                    /* PROCESS TERMINATED */
                    std::cout << "time " << time << "ms: Process " << finished_process.getLabel() << " terminated [Q";
                    if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                    else {
                        for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                            std::cout << " " << it->getLabel();
                        }
                        std::cout << "]" << std::endl;
                    }

                    num_terminated++;

                } else {
                    /* PROCESS GETS ADDED TO WAITING QUEUE */
                    int io_burst_len = io_burst_index[finished_process.getLabel()];
                    finished_process.setTimeStartedWaiting(time + io_burst_len + (context_switch / 2));
                    waiting.push_back(finished_process);
                    waiting.sort(fcfs_compare_wait);
                    if (time < 10000) {
                        std::cout << "time " << time << "ms: Process " << finished_process.getLabel()
                            << " switching out of CPU; blocking on I/O until time "
                            << time + io_burst_len + (context_switch / 2) << "ms [Q";
                        if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                        else {
                            for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                                std::cout << " " << it->getLabel();
                            }
                            std::cout << "]" << std::endl;
                        }
                    }
                }

                time += (context_switch / 2);
            }

        } else if (cpu_burst_starts <= io_burst_completion && cpu_burst_starts <= next_arrival_time) {
            /*PROCESS STARTS USING CPU*/
    
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
            magic_number = 1;


        } else if (io_burst_completion <= cpu_burst_completion && io_burst_completion <= next_arrival_time) {
            /*IO BURST COMPLETION*/

            Process temp_process = waiting.front();
            waiting.pop_front();

            time = temp_process.getDoneWaiting();
            temp_process.setDoneWaiting(INT_MAX);
            temp_process.setTimeStartedWaiting(time);
            
            ready.push_back(temp_process);

            if (time < 10000) {
                std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << "  completed I/O; added to ready queue [Q";
                if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
                else {
                    for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
                        std::cout << " " << it->getLabel();
                    }
                    std::cout << "]" << std::endl;
                }
            }
            
            io_burst_index[temp_process.getLabel()]++;
            
        } else { 
            /*NEW ARRIVAL*/

            Process temp_process = processes[next_arrival];
            
            time = temp_process.getArrivalTime();
            temp_process.setTimeStarted(time);
            
            temp_process.setTimeStartedWaiting(time);                        
            ready.push_back(temp_process);

            next_arrival++;
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
  
        }
        
    }
;
}



// void rr(std::vector<Process> processes, int context_switch, int time_slice, std::ofstream &output, clock_t start) {

//     int timeCpuRunning = 0; 
//     int sumCpuTime_CPU = 0;
//     int numCpuBursts_CPU = 0;
//     int sumCpuTime_IO = 0;
//     int numCpuBursts_IO = 0;
//     int sumWaitTime_CPU = 0;
//     int sumWaitTime_IO = 0;
//     int num_preemptions_CPU = 0;
//     int num_preemptions_IO = 0;
//     clock_t end;

//     int time = 0;
//     std::cout << "time " << time << "ms: Simulator started for RR [Q <empty>]" << std::endl;

//     std::sort(processes.begin(), processes.end(), fcfs_compare_arr);
//     int next_arrival = 0;
//     int magic_number = 1;

//     std::list<Process> running;
//     std::list<Process> waiting;
//     std::list<Process> ready;
//     std::map<char, int> cpu_burst_index;
//     std::map<char, int> io_burst_index;
//     std::vector<Process> terminated;
    

//     for (unsigned int i = 0; i < processes.size(); i++) {
//         cpu_burst_index[processes[i].getLabel()] = 0;
//         io_burst_index[processes[i].getLabel()] = 0;

//         std::vector<int> temp_cpu_bursts = processes[i].getCpuBursts();
//         std::vector<int> temp_io_bursts = processes[i].getIoBursts();
        
//         if (processes[i].getCpuBound()) { 
//             for (unsigned int j = 0; j < temp_cpu_bursts.size(); j++) { sumCpuTime_CPU += temp_cpu_bursts[j]; }
//             numCpuBursts_CPU += temp_cpu_bursts.size();
//         } else {
//             for (unsigned int j = 0; j < temp_cpu_bursts.size(); j++) { sumCpuTime_IO += temp_cpu_bursts[j]; }
//             numCpuBursts_IO += temp_cpu_bursts.size();
//         }
//     }

//     int temp_total_cpu_time = sumCpuTime_CPU;
//     int temp_total_cpu_time_io = sumCpuTime_IO;

//     bool arrival_inside_cs = false;

//     while (processes.size() > terminated.size()) {
//         end = clock();
//         if ((double)(end - start)/CLOCKS_PER_SEC > 0.25){
//             break;
//         }
        
//          //Remove process from processes once finished, not arrived
//         // Figure out what happens next: arrival, finishing waiting state, or finishing running state

//         // next_arrival is the index of the next process to arrive from processes

//         /* the time when the next process arrives */
//         int next_arrival_time = INT_MAX;
//         if ((unsigned int)next_arrival < processes.size()) {
//             next_arrival_time = processes[next_arrival].getArrivalTime();
//         }

//         /* the time when the burst using the CPU finishes and is ready to leave the running queue */
//         int done_running_time = INT_MAX;
//         if (running.size() > 0) {
//             done_running_time = running.front().getDoneRunning();
//         }

//         /* the time when the burst blocking on IO finishes and is ready to leave the waiting queue */
//         int done_blocking_time = INT_MAX;
//         if (waiting.size() > 0) {
//             done_blocking_time = waiting.front().getDoneWaiting();
//         }
        
//         int next_time_slice = INT_MAX;
//         if (running.size() > 0) {
//             Process preempted_process = running.front();
//             next_time_slice = (magic_number * time_slice) + preempted_process.getTimeStarted();
//         }

//         std::list<int> events;
//         events.push_back(next_arrival_time);
//         events.push_back(done_running_time);
//         events.push_back(done_blocking_time);
//         events.push_back(next_time_slice);

//         int next_event = *std::min_element(events.begin(), events.end());
//         // int next_event = std::min_element({next_arrival_time, done_running_time, done_blocking_time, next_time_slice});
        
//         if (next_event == next_arrival_time && next_arrival_time != INT_MAX) {
//             // if (arrival_inside_cs) {
//                 /* NEW ARRIVAL */
//                 Process temp_process = processes[next_arrival];
//                 time = temp_process.getArrivalTime();
//                 temp_process.setTimeStartedWaiting(time);
//                 ready.push_back(temp_process);
//                 next_arrival += 1;
//                 if (time < 10000) {
//                     std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " arrived; added to ready queue [Q";
//                     if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
//                     else {
//                         for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
//                             std::cout << " " << it->getLabel();
//                         }
//                         std::cout << "]" << std::endl;
//                     }
//                 }
//             // }
//         }

//         if (running.size() == 0 && ready.size() > 0) {
//             /* MOVE NEXT PROCESS FROM READY STATE TO RUNNING STATE */
//             Process temp_process = ready.front();
//             ready.pop_front();

//             if (temp_process.getCpuBound()) {
//                 sumWaitTime_CPU += time - temp_process.getTimeStartedWaiting();
//             } else {
//                 sumWaitTime_IO += time - temp_process.getTimeStartedWaiting();
//             }

//             time += (context_switch / 2);

//             if (temp_process.getTimeRan() == 0) {
//                 /* the process never got cut off */
            
//                 int temp_burst_len = temp_process.getCpuBurst(cpu_burst_index[temp_process.getLabel()]);
//                 temp_process.setDoneRunning(time + temp_burst_len);
//                 temp_process.setTimeStarted(time);
//                 if (time < 10000) {
//                     std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " started using the CPU for "
//                             << temp_burst_len << "ms burst [Q";
//                     if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
//                     else {
//                         for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
//                             std::cout << " " << it->getLabel();
//                         }
//                         std::cout << "]" << std::endl;
//                     }
//                 }

//                 cpu_burst_index[temp_process.getLabel()] += 1;
                
//                 running.push_back(temp_process);

//             } else {
//                 /* finish the cpu burst that already started */

//                 if (temp_process.getCpuBound()) {
//                     temp_total_cpu_time += context_switch;
//                 } else {
//                     temp_total_cpu_time_io += context_switch;
//                 }

//                 int temp_burst_len = temp_process.getCpuBurst(cpu_burst_index[temp_process.getLabel()]);
//                 temp_process.setDoneRunning(time + temp_burst_len - temp_process.getTimeRan());
//                 temp_process.setTimeStarted(time);
//                 if (time < 10000) {
//                     std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " started using the CPU for remaining "
//                             << temp_burst_len - temp_process.getTimeRan() << "ms of " << temp_burst_len << "ms burst [Q";
//                     if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
//                     else {
//                         for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
//                             std::cout << " " << it->getLabel();
//                         }
//                         std::cout << "]" << std::endl;
//                     }
//                 }

//                 cpu_burst_index[temp_process.getLabel()] += 1;
//                 running.push_back(temp_process);
//             }
//         }

//         if (next_event == done_running_time) {
//             /* DONE WITH RUNNING STATE */
//             Process temp_process = running.front();
//             running.pop_front();

//             time = temp_process.getDoneRunning();
//             temp_process.setDoneRunning(INT_MAX);
//             temp_process.updateTimeRan(-1 * temp_process.getTimeRan());
//             timeCpuRunning += time - temp_process.getTimeStarted();
//             magic_number = 1;

//             if (temp_process.getCpuBursts().size() == (unsigned int)cpu_burst_index[temp_process.getLabel()]) {
//                 /* PROCESS TERMINATED */
//                 std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " terminated [Q";
//                 if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
//                 else {
//                     for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
//                         std::cout << " " << it->getLabel();
//                     }
//                     std::cout << "]" << std::endl;
//                 }
//                 terminated.push_back(temp_process);

//                 //MAKE SURE THIS BELONGS HERE
//                 time += (context_switch / 2);

//             } else {

//                 if (time < 10000) {
//                     std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " completed a CPU burst; "
//                             << temp_process.getCpuBursts().size() - cpu_burst_index[temp_process.getLabel()]
//                             << " bursts to go [Q";
//                     if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
//                     else {
//                         for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
//                             std::cout << " " << it->getLabel();
//                         }
//                         std::cout << "]" << std::endl;   
//                     }
//                 }


//                 int temp_burst_len = temp_process.getIOBurst(io_burst_index[temp_process.getLabel()]);
//                 temp_process.setDoneWaiting(time + temp_burst_len + (context_switch / 2));
                
//                 if (time < 10000) {
//                     std::cout << "time " << time << "ms: Process " << temp_process.getLabel()
//                           << " switching out of CPU; blocking on I/O until time "
//                           << time + temp_burst_len + (context_switch / 2) << "ms [Q";
//                     if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
//                     else {
//                         for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
//                             std::cout << " " << it->getLabel();
//                         }
//                         std::cout << "]" << std::endl;
//                     }
//                 }
//                 io_burst_index[temp_process.getLabel()] += 1;
//                 waiting.push_back(temp_process);
//                 waiting.sort(fcfs_compare_wait);
//                 //std::sort(waiting.begin(), waiting.end(), fcfs_compare_wait);
//                 time += (context_switch / 2);
//             }
//         }

//         if (next_event == next_time_slice) {
//             /* TIME SLICE EXPIRED */

//             if (ready.size() > 0) {
//                 /* preempt the current running process and send it to the back of the ready queue */
//                 Process preempted_process = running.front();
//                 preempted_process.updateTimeRan(magic_number * time_slice);
//                 running.erase(running.begin());
                
//                 int temp_burst_len = preempted_process.getCpuBurst(cpu_burst_index[preempted_process.getLabel()] - 1);
//                 time = preempted_process.getTimeStarted() + (magic_number * time_slice);
//                 timeCpuRunning += time - preempted_process.getTimeStarted();
                                    
//                 if (time < 10000) {
//                     std::cout << "time " << time << "ms: Time slice expired; preempting process " << preempted_process.getLabel() << " with " << temp_burst_len - preempted_process.getTimeRan() << "ms remaining [Q"; 
//                     if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
//                     else {
//                         for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
//                             std::cout << " " << it->getLabel();
//                         }
//                         std::cout << "]" << std::endl;
//                     }
//                 }
                
//                 time += (context_switch / 2);
//                 preempted_process.setTimeStartedWaiting(time);
//                 cpu_burst_index[preempted_process.getLabel()] -= 1;
//                 ready.push_back(preempted_process);
//                 magic_number = 1;
                
//                 if (preempted_process.getCpuBound()) {
//                     num_preemptions_CPU++;
//                 } else {
//                     num_preemptions_IO++;
//                 }

//             } else {
//                 time += time_slice;
//                 Process running_process = running.front();
//                 running_process.updateTimeRan(time_slice);
//                 magic_number++;
                
//                 if (time < 10000) {
//                     std::cout << "time " << time << "ms: Time slice expired; no preemption because ready queue is empty [Q <empty>]" << std::endl;
//                 }
//             }

//         } 
        
//         if (next_event == done_blocking_time) {
//             /* DONE WITH WAITING STATE */
//             Process temp_process = waiting.front();
//             waiting.pop_front();
//             //waiting.erase(waiting.begin());

//             time = temp_process.getDoneWaiting();
//             temp_process.setDoneWaiting(INT_MAX);
//             temp_process.setTimeStartedWaiting(time);
//             ready.push_back(temp_process);

//             if (time < 10000) {
//                 std::cout << "time " << time << "ms: Process " << temp_process.getLabel() << " completed I/O; added to ready queue [Q";
//                 if (ready.size() == 0) { std::cout << " <empty>]" << std::endl; }
//                 else {
//                     for (std::list<Process>::iterator it = ready.begin(); it != ready.end(); it++) {
//                         std::cout << " " << it->getLabel();
//                     }
//                     std::cout << "]" << std::endl;
//                 }
//             }
            
//         }

//         // if (time > 690 && time < 700) { std::cout << "  time: " << time << ", next_arrival in: " << time - next_arrival_time << ", context_switch: " << context_switch << std::endl; }
        
//         if (time + context_switch > next_arrival_time) { arrival_inside_cs = true; continue; }

//         arrival_inside_cs = false;

        
//     }

//     std::cout << "time " << time << "ms: Simulator ended for RR [Q <empty>]" << std::endl;

//     double cpu_utilization = 0;
//     if(time != 0){
//         cpu_utilization = ceil(1000 * ((double)timeCpuRunning / time) * 100) / 1000;
//     }
    
//     double cpuBurstAvg_CPU = 0;
//     if(numCpuBursts_CPU != 0){
//         cpuBurstAvg_CPU = ceil(1000 * (double)sumCpuTime_CPU / numCpuBursts_CPU) / 1000;
//     }

//     double cpuBurstAvg_IO = 0;
//     if(numCpuBursts_IO != 0){
//         cpuBurstAvg_IO = ceil(1000 * (double)sumCpuTime_IO / numCpuBursts_IO) / 1000;
//     }

//     double cpuBurstAvg_SUM = 0;
//     if(numCpuBursts_CPU + numCpuBursts_IO != 0){
//         cpuBurstAvg_SUM = ceil(1000.0 * (double)(sumCpuTime_CPU + sumCpuTime_IO) / (numCpuBursts_CPU + numCpuBursts_IO)) / 1000;
//     }

//     double IOBurstAvg_CPU = 0;
//     if(numCpuBursts_CPU != 0){
//         IOBurstAvg_CPU = ceil(1000 * (double)sumWaitTime_CPU / numCpuBursts_CPU) / 1000;
//     }

//     double IOBurstAvg_IO = 0;
//     if(numCpuBursts_IO != 0){
//         IOBurstAvg_IO = ceil(1000 * (double)sumWaitTime_IO / numCpuBursts_IO) / 1000;
//     }

//     double IOBurstAvg_SUM = 0;
//     if(numCpuBursts_CPU + numCpuBursts_IO != 0){
//         IOBurstAvg_SUM = ceil(1000.0 * (double)(sumWaitTime_CPU + sumWaitTime_IO) / (numCpuBursts_CPU + numCpuBursts_IO)) / 1000;
//     }

//     double turnaround = 0;
//     if (numCpuBursts_CPU + numCpuBursts_IO != 0) {
//         turnaround = (double)(sumWaitTime_CPU + sumWaitTime_IO) / (numCpuBursts_CPU + numCpuBursts_IO) + (double)(temp_total_cpu_time + temp_total_cpu_time_io) / (numCpuBursts_CPU + numCpuBursts_IO) + context_switch;
//     }

//     double temp_avg_cpu = 0;
//     if (numCpuBursts_CPU != 0) {
//         temp_avg_cpu = ceil(1000 * (double)temp_total_cpu_time / numCpuBursts_CPU) / 1000;
//     }

//     double temp_avg_io = 0;
//     if (numCpuBursts_IO != 0) {
//         temp_avg_io = ceil(1000 * (double)temp_total_cpu_time_io / numCpuBursts_IO) / 1000;
//     }
    
//     output << "Algorithm RR\n";
//     output << "-- CPU utilization: " << std::fixed << std::setprecision(3) << cpu_utilization << "%\n";
//     output << "-- average CPU burst time: " << cpuBurstAvg_SUM << " ms (" << cpuBurstAvg_CPU << " ms/" << cpuBurstAvg_IO << " ms)\n";
//     output << "-- average wait time: " << IOBurstAvg_SUM << " ms (" << IOBurstAvg_CPU << " ms/" << IOBurstAvg_IO << " ms)\n";
//     output << "-- average turnaround time: " << ceil(1000 * turnaround) / 1000 << " ms (" << temp_avg_cpu + IOBurstAvg_CPU + context_switch 
//            << " ms/" << temp_avg_io + IOBurstAvg_IO + context_switch << " ms)\n";
//     output << "-- number of context switches: " << numCpuBursts_CPU + numCpuBursts_IO + num_preemptions_IO + num_preemptions_CPU << " (" << numCpuBursts_CPU + num_preemptions_CPU << "/" << numCpuBursts_IO + num_preemptions_IO << ")\n";
//     output << "-- number of preemptions: " << num_preemptions_IO + num_preemptions_CPU << " (" << num_preemptions_CPU << "/" << num_preemptions_IO << ")\n";

// }

int main(int argc, char** argv) {
    if (argc != 9) {
        std::cerr << char(91) << "Usage" << char(93) << " " << *(argv) << " <n> <n_cpu> <seed> <lambda> <threshold>"
                  << std::endl;
        return EXIT_FAILURE;
    }

    /* ==== command line args for part I ==== */
    int n = atoi(*(argv + 1));
    int n_cpu = atoi(*(argv + 2));
    long seed = atol(*(argv + 3));
    double lambda = atof(*(argv + 4));
    int threshold = atoi(*(argv + 5));

    /* ==== command line args for part II === */
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
            //    std::cout << "--> CPU burst " << cpu_time << "ms --> I/O burst " << io_time << "ms" << std::endl;
            } else {
            //    std::cout << "--> CPU burst " << cpu_time << "ms" << std::endl;
            }
        }
        processes.push_back(Process(char(65 + i), n - i <= n_cpu, arrival_time, num_bursts, cpu_bursts, io_bursts));
    }

    std::cout << std::endl;
    std::cout << "<<< PROJECT PART II -- t_cs=" << context_switch << "ms; alpha=" << std::fixed << std::setprecision(2) << alpha << "; t_slice=" << time_slice
              << "ms >>>" << std::endl;


    std::ofstream output;
    output.open("simout.txt");

    clock_t start;
    start = clock();

    

    fcfs(processes, context_switch, output, start);
    std::cout << std::endl;
    sjf(processes, context_switch, lambda, alpha, output, start);
    std::cout << std::endl;
    srt(processes, context_switch, lambda, alpha, output, start);
    std::cout << std::endl;
    rr(processes, context_switch, time_slice, output);

    output.close();


    //clock_t end;
    //end = clock();
    //std::cout << start << ", " << end << " diff = " << end - start << std::endl;

    //clock_t final;
    //final = clock();
    //std::cout << (double)(final - start)/CLOCKS_PER_SEC << std::endl;


    return EXIT_SUCCESS;
}