#include <iostream>
#include "../include/instance.h"
#include "../include/batch_solver.h"
#include <filesystem>
int main(int   argc,
         char *argv[]) {
    string folderPath = "/home/22800033/jssp/jsp_8_8_small";
    string outputPath = "v_test_report";
    double timeLimit = 3600.0;
    int formula = 3; // Liao's Disjunctive
    string spec = "O";
    auto batchSolver = BatchSolver();
    bool useVTimeIndexed = true;
    ofstream opf("output_test_v.csv");
//    batchSolver.CPLEXBatch(folderPath, outputPath, formula, timeLimit, spec);
    for (int v = 71; v <= 347; v++) {
        double avg = 0;
        for (int i = 0; i < 10; i++) {
            auto time = batchSolver.GurobiBatch(folderPath, outputPath, formula, timeLimit, spec, v, true);
            avg += time[0];
        }
        avg /= 10;
        cout << v << ", " << avg/100 << endl;
        opf << v << "," << avg/100 << "\n";
    }

}

