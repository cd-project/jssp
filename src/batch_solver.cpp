//
// Created by cuong on 25/07/2023.
//
#include "../include/batch_solver.h"
#include "../include/solver.h"
BatchSolver::BatchSolver() = default;

void BatchSolver::CPLEXBatch(string folderPath, string outputPath, int formula, double timeLimit, string spec) {
    vector<string> fileList;
    for (const auto& entry : filesystem::directory_iterator(folderPath)) {
        if (entry.is_regular_file()) {
            fileList.push_back(entry.path().filename().string());
        }
    }
    string opF;
    switch (formula) {
        case 1:
            opF = "cplex_disjunctive_batch_" + outputPath;
            break;
        case 2:
            opF = "cplex_liao_disjunctive_batch_" + outputPath;
            break;
        case 3:
            opF = "cplex_time_indexed_batch_" + outputPath;
            break;
        case 4:
            opF = "cplex_rank_based_batch_" + outputPath;
            break;
        default:
            return;
    }
    ofstream file(opF);
    struct InstanceResult {
        string InstanceName;
        int LB;
        int UB;
        double RunningTime;
        bool IsOptimal;
        double Gap;
    };
    vector<InstanceResult> data;
    if (file.is_open()) {
        file << "Instance name,Lower bound,Upper bound,Time,Optimal,Gap\n";
        for (int i = 0; i < fileList.size(); i++) {
            auto fName = folderPath + "/" + fileList[i];
            cout << "Current instance name: " << fName << endl;
            auto inst = Instance(fName, spec);
            auto solver = Solver(inst);
            tuple<int, int, double, bool, double> res;
            InstanceResult iR;
            switch (formula) {
                case 1: // Disjunctive
                    res = solver.CPLEXDisjunctiveSolver(timeLimit);
                    iR = (InstanceResult{fileList[i], get<0>(res), get<1>(res), get<2>(res), get<3>(res), get<4>(res)});
                    break;
                case 2:
                    res = solver.CPLEXLiaoDisjunctiveSolver(timeLimit);
                    iR = (InstanceResult{fileList[i], get<0>(res), get<1>(res), get<2>(res), get<3>(res), get<4>(res)});
                    break;
                case 3:
                    res = solver.CPLEXTimeIndexedSolver(timeLimit);
                    iR = (InstanceResult{fileList[i], get<0>(res), get<1>(res), get<2>(res), get<3>(res), get<4>(res)});
                    break;
                case 4:
                    res = solver.CPLEXRankBasedSolver(timeLimit);
                    iR = (InstanceResult{fileList[i], get<0>(res), get<1>(res), get<2>(res), get<3>(res), get<4>(res)});
                    break;
            }
            file << iR.InstanceName << ","
                 << iR.LB << ","
                 << iR.UB << ","
                 << iR.RunningTime << ","
                 << iR.IsOptimal << ","
                 << iR.Gap << "\n";
        }
        file.close();
        std::cout << "CSV file created successfully." << std::endl;
    } else {
        std::cerr << "Error creating CSV file." << std::endl;
    }
}

vector<double> BatchSolver::GurobiBatch(std::string folderPath, std::string outputPath, int formula, double timeLimit, string spec, int v, bool useVTimeIndexed) {
    vector<double> runtimes;
    vector<string> fileList;
    for (const auto& entry : filesystem::directory_iterator(folderPath)) {
        if (entry.is_regular_file()) {
            fileList.push_back(entry.path().filename().string());
        }
    }
    string opF;
    cout << formula << endl;
    switch (formula) {
        case 1:
            opF = "gurobi_disjunctive_batch_" + outputPath;
            break;
        case 2:
            opF = "gurobi_liao_disjunctive_batch_" + outputPath;
            break;
        case 3:
            opF = "gurobi_time_indexed_batch_" + outputPath;
            break;
        case 4:
            opF = "gurobi_rank_based_batch_" + outputPath;
            break;
        default:
            return runtimes;
    }
    cout << opF << endl;
    ofstream file(opF);
    ofstream oF("OUT.txt");
    struct InstanceResult {
        string InstanceName;
        int LB;
        int UB;
        double RunningTime;
        bool IsOptimal;
        double Gap;
    };
    vector<InstanceResult> data;
    if (file.is_open()) {
        file << "Instance name,Lower bound,Upper bound,Time,Optimal,Gap\n";
        for (int i = 0; i < fileList.size(); i++) {
            auto fName = folderPath + "/" + fileList[i];
            cout << "Current instance name: " << fName << endl;
            auto inst = Instance(fName, spec);

            auto solver = Solver(inst);
            tuple<int, int, double, bool, double> res;
            InstanceResult iR;
            switch (formula) {
                case 1: // Disjunctive
                    res = solver.GurobiDisjunctiveSolver(timeLimit);
                    iR =(InstanceResult{fileList[i], get<0>(res), get<1>(res), get<2>(res), get<3>(res), get<4>(res)});
                    break;

                case 2:
                    res = solver.GurobiLiaoDisjunctiveSolver(timeLimit);
                    iR =(InstanceResult{fileList[i], get<0>(res), get<1>(res), get<2>(res), get<3>(res), get<4>(res)});
                    break;
                case 3:
                    res = solver.GurobiTimeIndexedSolver(timeLimit, v, useVTimeIndexed);
                    iR =(InstanceResult{fileList[i], get<0>(res), get<1>(res), get<2>(res), get<3>(res), get<4>(res)});
                    break;
                case 4:
                    res = solver.GurobiRankBasedSolver(timeLimit);
                    iR =(InstanceResult{fileList[i], get<0>(res), get<1>(res), get<2>(res), get<3>(res), get<4>(res)});
                    break;
            }
            file << iR.InstanceName << ","
                 << iR.LB << ","
                 << iR.UB << ","
                 << iR.RunningTime << ","
                 << iR.IsOptimal << ","
                 << iR.Gap << "\n";
            runtimes.push_back(iR.RunningTime);
        }

        file.close();
        std::cout << "CSV file created successfully." << std::endl;
    } else {
        std::cerr << "Error creating CSV file." << std::endl;
    }
    return runtimes;
}

