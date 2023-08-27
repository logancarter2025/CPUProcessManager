#ifndef OPSYSPROJ1_PROCESS_H
#define OPSYSPROJ1_PROCESS_H

#include <iostream>
#include <vector>
#include <climits>

class Process {
public:
    Process(char l, bool cb, int arr_time, int n_bursts, std::vector<int> &cpu, std::vector<int> &io);
    
    void printProcess(std::ostream& output);

    /* ======= get() functions ======= */

    char getLabel() {
        return label;
    }

    double getTau(){
        return tau;
    }

    bool getCpuBound() {
        return cpu_bound;
    }

    int getArrivalTime() {
        return arrival_time;
    }

    int getNumBursts() {
        return num_bursts;
    }

    std::vector<int>& getCpuBursts() {
        return cpu_bursts;
    }

    std::vector<int>& getIoBursts() {
        return io_bursts;
    }

    bool getCpuAvailability() {
        return cpu_availability;
    }

    bool hasArrived(int& time) {
        return time >= arrival_time;
    }

    int getCpuBurst(int index) {
        return cpu_bursts[index];
    }

    int getIOBurst(int index) {
        return io_bursts[index];
    }

    int getDoneRunning() {
        return done_running;
    }

    int getDoneWaiting() {
        return done_waiting;
    }

    int getTimeStarted() {
        return time_started;
    }

    int getTimeRan(){
        return time_ran;
    }
    
    int getTimeStartedWaiting() {
        return time_started_waiting;
    }

    int getTimeStartRunning() {
        return time_start_running;
    }
    
        
    /* ======= set() functions ======= */
    
    void setCpuAvailability(bool state) {
        cpu_availability = state;
    }

    void setDoneRunning(int time) {
        done_running = time;
    }

    void setDoneWaiting(int time) {
        done_waiting = time;
    }

    void setTau(double t){
        tau = t;
    }

    //t is current time minus time started
    void updateTimeRan(int t) {
        time_ran += t;
    }

    void setTimeStarted(int t) {
        time_started = t;
    }

    void setTimeStartedWaiting(int t) {
        time_started_waiting = t;
    }

    void setTimeStartRunning(int t) {
        time_start_running = t;
    }

    

    


private:
    char label;
    bool cpu_bound;
    int arrival_time;
    int num_bursts;
    std::vector<int> cpu_bursts;
    std::vector<int> io_bursts;
    bool cpu_availability;

    int done_running;
    int done_waiting;
    double tau;
    int time_ran; //Total time in the CPU
    int time_started; //Time entered into current CPU Burst
    int time_started_waiting;

    int time_start_running; 

};


bool fcfs_compare_arr(Process& p1, Process& p2);
bool fcfs_compare_wait(Process& p1, Process& p2);
bool sjf_compare_tau(Process& p1, Process& p2);
bool srt_compare_time_remaining(Process &p1, Process &p2);


#endif //OPSYSPROJ1_PROCESS_H
