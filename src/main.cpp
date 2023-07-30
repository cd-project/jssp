#include <iostream>
#include "../include/instance.h"
#include "../include/batch_solver.h"
#include <filesystem>
int main(int   argc,
         char *argv[]) {
    string folderPath = "/home/cuong/CLionProjects/JSSP/jsp50_3";
    string outputPath = "jsp_50_3_all";
    double timeLimit = 300.0;
    int formula = 3; // Liao's Disjunctive
    string spec = "O";
    auto batchSolver = BatchSolver();
    bool useVTimeIndexed = true;
//    ofstream opf("output_test_v.txt");
//    batchSolver.CPLEXBatch(folderPath, outputPath, formula, timeLimit, spec);
//    for (int v = 71; v <= 150; v++) {
//        double avg = 0;
//        for (int i = 0; i < 1; i++) {
//            avg += time[0];
//        }
//        avg /= 1;
//        cout << v << ", " << avg/100 << endl;
//        opf << v << "," << avg/100 << "\n";
//    }
    auto time = batchSolver.GurobiBatch(folderPath, outputPath, formula, timeLimit, spec, 150, false);

}

