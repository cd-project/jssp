//
// Created by cuong on 25/07/2023.
//

#ifndef JSSP_BATCH_SOLVER_H
#define JSSP_BATCH_SOLVER_H

#include <filesystem>

#include "ilcplex/ilocplex.h"
#include "gurobi_c++.h"
#include "instance.h"

class BatchSolver {
public:
    BatchSolver();
    void CPLEXBatch(string folderPath, string outputPath, int formula, double tiLim, string spec);
    vector<double> GurobiBatch(string folderPath, string outputPath, int formula, double tiLim, string spec, int v, bool useVTimeIndexed);
};
#endif //JSSP_BATCH_SOLVER_H

