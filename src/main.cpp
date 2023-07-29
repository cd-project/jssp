#include <iostream>
#include "../include/instance.h"
#include "../include/batch_solver.h"
#include <filesystem>
int main(int   argc,
         char *argv[]) {
    string folderPath = "/home/cuong/CLionProjects/JSSP/jsp_3_3_small";
    string outputPath = "ta_batch";
    double timeLimit = 3600.0;
    int formula = 3 ; // Liao's Disjunctive
    string spec = "O";
    auto batchSolver = BatchSolver();
    batchSolver.CPLEXBatch(folderPath, outputPath, formula, timeLimit, spec);

//    batchSolver.GurobiBatch(folderPath, outputPath, formula, timeLimit, spec);
}
