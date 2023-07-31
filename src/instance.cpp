//
// Created by cuong on 19/07/2023.
//
#include <utility>
#include <fstream>
#include <iostream>
#include "../include/instance.h"

vector<string> SplitStringWithDelimiter(string s, string delimiter) {
    vector<string> returnValue;
    string::size_type start = 0;
    string::size_type end = s.find(delimiter);

    while(end != string::npos) {
        returnValue.push_back(s.substr(start, end-start));
        start = end + 1;
        end = s.find(delimiter, start);
    }

    returnValue.push_back(s.substr(start));
    return returnValue;
}



Job::Job(int jobIndex, vector<tuple<int,int>> machineRuntimes) {
    this->Index = jobIndex;
    this->MachineRuntimes = std::move(machineRuntimes);
}

// Initializer for Taillard's specification.
Instance::Instance(string filePath, string spec) {
    ifstream dataFile(filePath);
    string line;
    if (spec == "Taillard") {
        bool startReadOverallStat = true;
        int numTimesLine = 0;
        int numMachsLine = 0;
        vector<vector<int>> timeArr;
        vector<vector<int>> machineArr;
        int cnt = 0;
        while(getline(dataFile, line)) {
            if (startReadOverallStat) {
                auto overallStatStr = SplitStringWithDelimiter(line, "\t");

                cout << overallStatStr[0] << endl;
                cout << overallStatStr[1] << endl;

                NumberOfJobs = stoi(overallStatStr[0]);
                NumberOfMachines = stoi(overallStatStr[1]);

                startReadOverallStat = false;
                continue;
            }
            if (numTimesLine < NumberOfJobs) {
                auto times = SplitStringWithDelimiter(line, "\t");
                vector<int> row;
                for (int i = 0; i < times.size(); i++) {
                    row.push_back(stoi(times[i]));
                }
                timeArr.push_back(row);
                numTimesLine++;
                continue;
            }


            if (numMachsLine < NumberOfJobs) {
                auto machineSeq = SplitStringWithDelimiter(line, "\t");
                vector<int> row;
                for (int i = 0; i < machineSeq.size(); i++) {
                    row.push_back(stoi(machineSeq[i]));
                }
                machineArr.push_back(row);
                numMachsLine++;
                continue;
            }
        }
        // Taillard's spec: Corresponding position for machine - operate time on that machine

        for (int i = 0; i < machineArr.size(); i++) {
            int jobIndex = i;
            vector<tuple<int, int>> machineRuntimes;
            for (int j = 0; j < machineArr[i].size(); j++) {
                auto machineIndex = machineArr[i][j]-1;
                auto timeForThisMachine = timeArr[i][j];
                tuple<int, int> machineAndTime;
                machineAndTime = make_tuple(machineIndex, timeForThisMachine);
                machineRuntimes.push_back(machineAndTime);
            }
            auto job = Job(jobIndex, machineRuntimes);
            Jobs.push_back(job);
        }
        vector<vector<int>> procTime(NumberOfMachines, vector<int>(NumberOfJobs));

        // Processing time Pij
        for (int j = 0; j < Jobs.size(); j++) {
            auto machineProcTimeInfo = Jobs[j].MachineRuntimes;
            for (const auto& tuple:machineProcTimeInfo) {
                auto i = get<0>(tuple);
                auto time = get<1>(tuple);
                procTime[i][j] = time;
            }
        }

        ProcessingTime = procTime;
    } else {
        bool startReadOverallStat = true;
        bool startReadMach = false;
        int jobIndex = 0;
        while(getline(dataFile, line)) {
            if (startReadOverallStat) {
                auto overallStatStr = SplitStringWithDelimiter(line, "\t");

                cout << overallStatStr[0] << endl;
                cout << overallStatStr[1] << endl;

                NumberOfJobs = stoi(overallStatStr[0]);
                NumberOfMachines = stoi(overallStatStr[1]);

                startReadOverallStat = false;
                startReadMach = true;
                continue;
            }
            if (startReadMach) {
                auto inputArray = SplitStringWithDelimiter(line, "\t");
                vector<tuple<int, int>> tuples;
                for (size_t i = 0; i < inputArray.size(); i += 2) {
                    // Check if there are at least two elements left in the array
                    if (i + 1 < inputArray.size()) {
                        // Create a tuple from the current two elements
                        tuple<int, int> currentTuple = make_tuple(stoi(inputArray[i]), int(stoi(inputArray[i + 1])));
                        // Add the tuple to the vector
                        tuples.push_back(currentTuple);
                    }
                }
                auto job = Job(jobIndex, tuples);
                jobIndex++;
                Jobs.push_back(job);
            }

        }
        vector<vector<int>> procTime(NumberOfMachines, vector<int>(NumberOfJobs));

        // Processing time Pij
        for (int j = 0; j < Jobs.size(); j++) {
            auto machineProcTimeInfo = Jobs[j].MachineRuntimes;
            for (const auto& tuple:machineProcTimeInfo) {
                auto i = get<0>(tuple);
                auto time = get<1>(tuple);
                procTime[i][j] = time;
            }
        }

        ProcessingTime = procTime;
    }

}

void Instance::PrintInstanceInfo() {
    cout << "Number of jobs: " << NumberOfJobs << endl;
    cout << "Number of machines: " << NumberOfMachines << endl;

    cout << "Jobs info:" << endl;
    for (int i = 0; i < Jobs.size(); i++) {
        cout << Jobs[i].Index << " [ ";
        for (const auto& tuple : Jobs[i].MachineRuntimes) {
            cout << "(" << get<0>(tuple) << ", " << get<1>(tuple) << ") ";
        }
        cout << "]" << endl;
    }

    cout << "Processing P[i][j] info:" << endl;
    for (int i = 0; i < ProcessingTime.size(); i++) {
        for (int j = 0; j < ProcessingTime[i].size(); j++) {
            cout << ProcessingTime[i][j] << " ";
        }
        cout << endl;
    }
}
