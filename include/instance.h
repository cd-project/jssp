#ifndef JSSP_INSTANCE_H
#define JSSP_INSTANCE_H

#include <vector>
#include <string>
#include <tuple>
using namespace std;
class Job {
public:
    int Index;
    vector<tuple<int, int>> MachineRuntimes;
    explicit Job(int jobIndex, vector<tuple<int,int>> machineRuntimes);
};
class Instance {
public:
    int NumberOfJobs;
    int NumberOfMachines;
    vector<Job> Jobs;
    vector<vector<int>> ProcessingTime;

    explicit Instance(string filePath, string spec);
    void PrintInstanceInfo();

};
#endif //JSSP_INSTANCE_H
