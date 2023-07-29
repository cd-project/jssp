//
// Created by cuong on 19/07/2023.
//

#ifndef JSSP_SOLVER_H
#define JSSP_SOLVER_H

#include <filesystem>

#include "ilcplex/ilocplex.h"
#include "gurobi_c++.h"
#include "instance.h"

class Solver {
public:
    Instance instance;

    explicit Solver(Instance ins): instance{std::move(ins)} {}

    // CPLEX solver
    tuple<int, int, double, bool, double> CPLEXDisjunctiveSolver(double timeLimit);
    tuple<int, int, double, bool, double> CPLEXLiaoDisjunctiveSolver(double timeLimit);
    tuple<int, int, double, bool, double> CPLEXTimeIndexedSolver(double timeLimit);
    tuple<int, int, double, bool, double> CPLEXRankBasedSolver(double timeLimit);

    // Gurobi solver
    tuple<int, int, double, bool, double> GurobiDisjunctiveSolver(double timeLimit);
    tuple<int, int, double, bool, double> GurobiLiaoDisjunctiveSolver(double timeLimit);
    tuple<int, int, double, bool, double> GurobiTimeIndexedSolver(double timeLimit);
    tuple<int, int, double, bool, double> GurobiRankBasedSolver(double timeLimit);

};
#endif //JSSP_SOLVER_H
