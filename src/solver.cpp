//
// Created by cuong on 19/07/2023.
//
#include "../include/solver.h"
#include <chrono>
#include <filesystem>

using namespace std;


tuple<int, int, double, bool, double> Solver::CPLEXDisjunctiveSolver(double timeLimit) {
    IloEnv env;
    IloModel model(env);
    IloCplex cplex(model);
    cplex.setParam(IloCplex::TiLim, timeLimit);
    auto J = instance.NumberOfJobs;
    auto M = instance.NumberOfMachines;
    auto jobs = instance.Jobs;
    auto p = instance.ProcessingTime;
    cout << "J is: " << J << ", " << "M is: " << M << endl;
    // big V variable
    int V = 0;

    // V = sum(pij)
    for (int i = 0; i < p.size(); i++) {
        for (int j = 0; j < p[i].size(); j++) {
            V += p[i][j];
        }
    }
    cout << "Big V number: " << V << endl;

    typedef IloArray<IloNumVarArray> NumVarMat;
    typedef IloArray<NumVarMat> NumVar3Mat;
    typedef IloArray<IloNumArray>NumResMat;
    typedef IloArray<NumResMat> NumRes3Mat;

    // x[i][j]
    IloArray<IloNumVarArray> x(env, M);
    for (int i = 0; i < M; i++) {
        x[i] = IloNumVarArray(env, J, 0, IloInfinity, ILOINT);
    }

    IloArray<IloArray<IloNumVarArray>> z(env, M);
    for (int i = 0; i < M; ++i) {
        z[i] = IloArray<IloNumVarArray>(env, J);
        for (int j = 0; j < J; ++j) {
            z[i][j] = IloNumVarArray(env, J);
            for (int k = 0; k < J; k++) {
                if (j < k) {
                    z[i][j][k] = IloNumVar(env, 0, 1, ILOINT);
              }
            }
        }
    }


    // Constraint (2): x[i][j] >= 0 forall j in J, i in M
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < J; j++) {
            model.add(x[i][j] >= 0);
        }
    }

    // Constraint (3): x[currentMachine, currentJob] >= x[lastMachine, currentJob] + processing_time[lastMachine, currentJob]
    for (int j = 0; j < jobs.size(); j++) {
        auto currentJob = jobs[j].Index;
        for (int i = 1; i < M; i++) {
            auto currentMachine = get<0>(jobs[j].MachineRuntimes[i]);
            auto lastMachine = get<0>(jobs[j].MachineRuntimes[i-1]);
            auto lastMachineProcessingTime = get<1>(jobs[j].MachineRuntimes[i-1]);
            model.add(x[currentMachine][currentJob] >= x[lastMachine][currentJob] + lastMachineProcessingTime);
            cout << "Constraint added: " << "x[" << currentMachine << "]["  << currentJob << "] >= x[" << lastMachine << "][" << currentJob << "] + " << lastMachineProcessingTime << endl;
        }
    }

    // Constraint (4)/(5): x[i][j] >= x[i][k] + p[i][k] - V * z[i][j][k]
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < J; j++) {
            for (int k = 0; k < J; k++) {
                if (j < k) {
                    model.add(x[i][j] >= x[i][k] + p[i][k] - V * z[i][j][k]);
                    model.add(x[i][k] >= x[i][j] + p[i][j] - V * (1 - z[i][j][k]));
                }
            }
        }
    }

//    for (int i = 0; i < M; i++) {
//        for (int j = 0; j < J; j++) {
//            for (int k = 0; k < J; k++) {
//                model.add(0 <= z[i][j][k] <= 1);
//            }
//        }
//    }

    IloNumVar C(env);
    for (int j = 0; j < J; j++) {
        auto lastMachineOfJob = get<0>(jobs[j].MachineRuntimes[jobs[j].MachineRuntimes.size()-1]);
        auto lastMachineProcessingTime = get<1>(jobs[j].MachineRuntimes[jobs[j].MachineRuntimes.size()-1]);
        model.add(C >= x[lastMachineOfJob][j] + lastMachineProcessingTime);
        cout << "Constraint added: " << "C >= x[" << lastMachineOfJob << "][" << j << "] + " << lastMachineProcessingTime << endl;
    }
    model.add(IloMinimize(env, C));
    auto startTime = std::chrono::high_resolution_clock::now();
    cplex.solve();
    auto endTime = std::chrono::high_resolution_clock::now();

    cplex.exportModel("model.lp");
    cout <<"Objective value: " << cplex.getObjValue() << endl;
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    // Print the execution time
    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;
    bool isOptimal = false;
    if (cplex.getStatus() == IloAlgorithm::Optimal) {
        isOptimal = true;
    }
    tuple<int, int, double, bool, double> result;
    result = make_tuple(cplex.getBestObjValue(), cplex.getObjValue(), duration.count()/1000, isOptimal, cplex.getMIPRelativeGap());

    cplex.end();
    model.end();
    env.end();

    return result;
}

tuple<int, int, double, bool, double> Solver::CPLEXLiaoDisjunctiveSolver(double timeLimit) {
    IloEnv env;
    IloModel model(env);
    IloCplex cplex(model);
    cplex.setParam(IloCplex::TiLim, timeLimit);

    auto J = instance.NumberOfJobs;
    auto M = instance.NumberOfMachines;
    auto jobs = instance.Jobs;
    auto p = instance.ProcessingTime;
//    cout << "J is: " << J << ", " << "M is: " << M << endl;
    // big V variable
    int V = 0;

    // V = sum(pij)
    for (int i = 0; i < p.size(); i++) {
        for (int j = 0; j < p[i].size(); j++) {
            V += p[i][j];
        }
    }
//    cout << "Big V number: " << V << endl;

    typedef IloArray<IloNumVarArray> NumVarMat;
    typedef IloArray<NumVarMat> NumVar3Mat;
    typedef IloArray<IloNumArray>NumResMat;
    typedef IloArray<NumResMat> NumRes3Mat;

    // x[i][j]
    IloArray<IloNumVarArray> x(env, M);
    for (int i = 0; i < M; i++) {
        x[i] = IloNumVarArray(env, J, 0, IloInfinity, ILOINT);
    }

    IloArray<IloArray<IloNumVarArray>> z(env, M);
    for (int i = 0; i < M; ++i) {
        z[i] = IloArray<IloNumVarArray>(env, J);
        for (int j = 0; j < J; ++j) {
            z[i][j] = IloNumVarArray(env, J, 0, 1, ILOINT);
        }
    }

    IloArray<IloArray<IloNumVarArray>> q(env, M);
    for (int i = 0; i < M; ++i) {
        q[i] = IloArray<IloNumVarArray>(env, J);
        for (int j = 0; j < J; ++j) {
            q[i][j] = IloNumVarArray(env, J, 0, IloInfinity, ILOINT);
        }
    }
    // Constraint (2): x[i][j] >= 0 forall j in J, i in M
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < J; j++) {
            model.add(x[i][j] >= 0);
        }
    }

    // Constraint (3): x[currentMachine, currentJob] >= x[lastMachine, currentJob] + processing_time[lastMachine, currentJob]
    for (int j = 0; j < jobs.size(); j++) {
        auto currentJob = jobs[j].Index;
        cout << "On job: " << currentJob << endl;
        for (int i = 1; i < M; i++) {
            auto currentMachine = get<0>(jobs[j].MachineRuntimes[i]);
            auto lastMachine = get<0>(jobs[j].MachineRuntimes[i-1]);
            auto lastMachineProcessingTime = get<1>(jobs[j].MachineRuntimes[i-1]);
            model.add(x[currentMachine][currentJob] >= x[lastMachine][currentJob] + lastMachineProcessingTime);
//            cout << "Constraint added: " << "x[" << currentMachine << "]["  << currentJob << "] >= x[" << lastMachine << "][" << currentJob << "] + " << lastMachineProcessingTime << endl;
        }
    }

    // Constraint (4): x[i][j] >= x[i][k] + p[i][k] - V * z[i][j][k]
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < J; j++) {
            for (int k = 0; k < J; k++) {
                if (j < k) {
                    model.add(q[i][j][k] == V * z[i][j][k] + (x[i][j] - x[i][k] - p[i][k]));
                    model.add(q[i][j][k] <= V - p[i][j] - p[i][k]);
                }
            }
        }
    }

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < J; j++) {
            for (int k = 0; k < J; k++) {
                model.add(0 <= z[i][j][k] <= 1);
            }
        }
    }

    IloNumVar C(env);
    for (int j = 0; j < J; j++) {
        auto lastMachineOfJob = get<0>(jobs[j].MachineRuntimes[jobs[j].MachineRuntimes.size()-1]);
        auto lastMachineProcessingTime = get<1>(jobs[j].MachineRuntimes[jobs[j].MachineRuntimes.size()-1]);
        cout << "Last machine of job number " << j << " is: " << lastMachineOfJob << endl;
        cout << "Last processing time of job number " << j << " is: " << lastMachineProcessingTime << endl;
        model.add(C >= x[lastMachineOfJob][j] + lastMachineProcessingTime);
        cout << "Constraint added: " << "C >= x[" << lastMachineOfJob << "][" << j << "] + " << lastMachineProcessingTime << endl;
    }
    model.add(IloMinimize(env, C));
    auto startTime = std::chrono::high_resolution_clock::now();
    cplex.solve();
    auto endTime = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    // Print the execution time
    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;
    bool isOptimal = false;
    if (cplex.getStatus() == IloAlgorithm::Optimal) {
        isOptimal = true;
    }
    tuple<int, int, double, bool, double> result;
    result = make_tuple(cplex.getBestObjValue(), cplex.getObjValue(), duration.count()/1000, isOptimal, cplex.getMIPRelativeGap());

    cplex.end();
    model.end();
    env.end();

    return result;
}

tuple<int, int, double, bool, double> Solver::CPLEXRankBasedSolver(double timeLimit) {
    IloEnv env;
    IloModel model(env);
    IloCplex cplex(model);
    cplex.setParam(IloCplex::TiLim, timeLimit);

    auto J = instance.NumberOfJobs;
    auto M = instance.NumberOfMachines;
    auto jobs = instance.Jobs;
    auto p = instance.ProcessingTime;

    int V = 0;
    // V = sum(pij)
    for (int i = 0; i < p.size(); i++) {
        for (int j = 0; j < p[i].size(); j++) {
            V += p[i][j];
        }
    }
    // Variable x[i][j][k]
    IloArray<IloArray<IloNumVarArray>> x(env, M);
    for (int i = 0; i < M; i++) {
        x[i] = IloArray<IloNumVarArray>(env, J);
        for (int j = 0; j < J; j++) {
            x[i][j] = IloNumVarArray(env, J, 0, 1, ILOINT);
        }
    }

    // Variable h[i][k]
    IloArray<IloNumVarArray> h(env, M);
    for(int i = 0; i < M; i++) {
        h[i] = IloNumVarArray(env, J, 0, IloInfinity, ILOINT);
    }

    // Variable r[i][j][k]
    // Variable r[i][j][k]
    std::vector<std::vector<std::vector<int>>> r(M, std::vector<std::vector<int>>(J, std::vector<int>(M)));

    for (int j = 0; j < jobs.size(); j++) {
        auto jobIndex = j;
        auto machRuntimes = jobs[j].MachineRuntimes;
        for (int idx = 0; idx < machRuntimes.size(); idx++) {
            auto kOp = idx;
            auto machIdx = get<0>(machRuntimes[idx]);
            r[machIdx][j][kOp] = 1;
            cout << "On mach: " << machIdx << ", job: " << jobIndex << ", op: " << kOp << endl;
        }
    }
   // Constraint 17
   for (int i = 0; i < M; i++) {
       for (int k = 0; k < J; k++) {
           IloExpr constr17(env);
           for (int j = 0; j < J; j++) {
               constr17 += x[i][j][k];
           }
           model.add(constr17 == 1);
           constr17.end();
       }
   }

   // Constraint 18
   for (int i = 0; i < M; i++) {
       for (int j = 0; j < J; j++) {
           IloExpr constr18(env);
           for (int k = 0; k < J; k++) {
               constr18 += x[i][j][k];
           }
           model.add(constr18 == 1);
       }
   }

   // Constraint 19
   for (int k = 0; k < J-1; k++) {
       for (int i = 0; i < M; i++) {
           IloExpr sumJ(env);
           for (int j = 0; j < J; j++) {
               sumJ += p[i][j] * x[i][j][k];
           }
           model.add(h[i][k] + sumJ <= h[i][k+1]);
           sumJ.end();
       }
   }

   // Constraint 20
//   for (int i = 0; i < M; i++) {
//       IloExpr rhBase(env);
//       IloExpr rp(env);
//       IloExpr rxBase(env);
//       IloExpr rxPlus(env);
//       IloExpr rhPlus(env);
//       for (int j = 0; j < J; j++) {
//           for (int l = 0; l < M - 1; l++) {
//               for (int k = 0; k < J; k++) {
//                   rhBase += r[i][j][l] * h[i][k];
//                   rxBase += r[i][j][l] * x[i][j][k];
//               }
//               rp += r[i][j][l] * p[i][j];
//               // kAp: k with apostrophe k'
//               for (int kAp = 0; kAp < J; kAp++) {
//                   rxPlus += r[i][j][l+1] * x[i][j][kAp];
//                   rhPlus += r[i][j][l+1] * h[i][kAp];
//               }
//           }
//
//       }
//       model.add(rhBase + rp <= V * (1 - rxBase) + V * (1 - rxPlus) + rhPlus);
//       rhBase.end();
//       rp.end();
//       rxBase.end();
//       rxPlus.end();
//       rhPlus.end();
//   }
    for (int j = 0; j < J; j++) {
        for (int k = 0; k < J; k++) {
            for (int kAp = 0; kAp < J; kAp++) {
                for (int l = 0; l < M-1; l++) {
                    IloExpr sum1(env);
                    IloExpr sum2(env);
                    IloExpr sum3(env);
                    IloExpr sum4(env);
                    IloExpr sum5(env);
                    for (int i = 0; i < M; i++) {
                        sum1 += r[i][j][l] * h[i][k];
                        sum2 += r[i][j][l] * p[i][j];
                        sum3 += r[i][j][l] * x[i][j][k];
                        sum4 += r[i][j][l+1] * x[i][j][kAp];
                        sum5 += r[i][j][l+1] * h[i][kAp];
                    }
                    model.add(sum1 + sum2 <= V * (1 - sum3) + V * (1 - sum4) + sum5);
                }
            }
        }
    }

    // Constraint 21
    IloNumVar C(env);
    for (int i = 0; i < M; i++) {
        IloExpr sumPx(env);
        for (int j = 0; j < J; j++) {
            sumPx += p[i][j] * x[i][j][J-1];
        }
        model.add(C >= h[i][J-1] + sumPx);
    }


   model.add(IloMinimize(env, C));
   auto startTime = std::chrono::high_resolution_clock::now();
   cplex.solve();
   auto endTime = std::chrono::high_resolution_clock::now();

   auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

   // Print the execution time
   std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;
   bool isOptimal = false;
   if (cplex.getStatus() == IloAlgorithm::Optimal) {
       isOptimal = true;
   }
   tuple<int, int, double, bool, double> result;
   result = make_tuple(cplex.getBestObjValue(), cplex.getObjValue(), duration.count()/1000, isOptimal, cplex.getMIPRelativeGap());

   cplex.end();
   model.end();
   env.end();

   return result;
}

tuple<int, int, double, bool, double> Solver::CPLEXTimeIndexedSolver(double timeLimit) {
    IloEnv env;
    IloModel model(env);
    IloCplex cplex(model);
    cplex.setParam(IloCplex::TiLim, timeLimit);

    auto J = instance.NumberOfJobs;
    auto M = instance.NumberOfMachines;
    auto jobs = instance.Jobs;
    auto p = instance.ProcessingTime;

    int V = 0;
    // V = sum(pij)
    for (int i = 0; i < p.size(); i++) {
        for (int j = 0; j < p[i].size(); j++) {
            V += p[i][j];
        }
    }
//
//    V = 100;
    // x[i][j][t]
    IloArray<IloArray<IloNumVarArray>> x(env, M);
    for (int i = 0; i < M; i++) {
        x[i] = IloArray<IloNumVarArray>(env, J);
        for (int j = 0; j < J; j++) {
            x[i][j] = IloNumVarArray(env, V, 0, 1, ILOINT);
        }
    }

    // Constraint 11
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < J; j++) {
            IloExpr constr11(env);
            for (int t  = 0; t < V; t++) {
                constr11 += x[i][j][t];
            }
            model.add(constr11 == 1);
            constr11.end();
;        }
    }

    IloNumVar C(env);
    // Constraint 12
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < J; j++) {
            IloExpr sumT(env);
            for (int t = 0; t < V; t++) {
                sumT += (t + p[i][j]) * x[i][j][t];
            }
            model.add(C >= sumT);
            sumT.end();
        }
    }

    // Constraint 13
    for (int i = 0; i < M; i++) {
        for (int t = 0; t < V; t++) {
            IloExpr sumjt(env);
            for (int j = 0; j < J; j++) {
                for (int tAp = t - p[i][j] + 1; tAp <= t; tAp++) {
                    if (tAp < 0) {
                        cout << t << " " << p[i][j]  << " " << tAp << endl;
                        continue;
                    } else {
                        sumjt += x[i][j][tAp];
                    }
                }
            }
            model.add(sumjt <= 1);
            sumjt.end();
        }
    }


    // Constraint 14
    for (int j = 0; j < J; j++) {
        for (int i  = 1; i < M; i++) {
            IloExpr lhs(env);
            IloExpr rhs(env);
            for (int t = 0; t < V; t++) {
                auto machAtHMinus1 = get<0>(jobs[j].MachineRuntimes[i-1]);
                lhs += (t + p[machAtHMinus1][j]) * x[machAtHMinus1][j][t];
                rhs += t * x[get<0>(jobs[j].MachineRuntimes[i])][j][t];
            }
            model.add(lhs <= rhs);
            lhs.end();
            rhs.end();
        }
    }

    model.add(IloMinimize(env, C));
    auto startTime = std::chrono::high_resolution_clock::now();
    cplex.solve();
    cplex.exportModel("model_ti.lp");
    auto endTime = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    // Print the execution time
    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;
    bool isOptimal = false;
    if (cplex.getStatus() == IloAlgorithm::Optimal) {
        isOptimal = true;
    }
    tuple<int, int, double, bool, double> result;
    result = make_tuple(cplex.getBestObjValue(), cplex.getObjValue(), duration.count()/1000, isOptimal, cplex.getMIPRelativeGap());

    cplex.end();
    model.end();
    env.end();

    return result;
}

tuple<int, int, double, bool, double> Solver::GurobiDisjunctiveSolver(double timeLimit) {
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);

    auto J = instance.NumberOfJobs;
    auto M = instance.NumberOfMachines;
    auto jobs = instance.Jobs;
    auto p = instance.ProcessingTime;

    int V = 0;

    // V = sum(pij)
    for (int i = 0; i < p.size(); i++) {
        for (int j = 0; j < p[i].size(); j++) {
            V += p[i][j];
        }
    }
    // Variable x
    auto** x = new GRBVar * [M];
    for (int i= 0; i < M; i++) {
        x[i] = reinterpret_cast<GRBVar *>(new GRBVar *[J]);
        for (int j = 0; j < J; j++) {
            x[i][j] = model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER);
        }
    }

    auto***z = new GRBVar** [M];
    for (int i = 0; i < M; i++) {
        z[i] = reinterpret_cast<GRBVar **>(new GRBVar **[J]);
        for(int j = 0; j < J; j++) {
            z[i][j] = reinterpret_cast<GRBVar *>(new GRBVar *[J]);
            for (int k = 0; k < J; k++) {
                z[i][j][k] = model.addVar(0, 1, 0, GRB_INTEGER);
            }
        }
    }

    // Constraint (2): x[i][j] >= 0 forall j in J, i in M
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < J; j++) {
            model.addConstr(x[i][j], GRB_GREATER_EQUAL, 0);
        }
    }

    // Constraint (3): x[currentMachine, currentJob] >= x[lastMachine, currentJob] + processing_time[lastMachine, currentJob]
    for (int j = 0; j < jobs.size(); j++) {
        auto currentJob = jobs[j].Index;
        cout << "On job: " << currentJob << endl;
        for (int i = 1; i < M; i++) {
            auto currentMachine = get<0>(jobs[j].MachineRuntimes[i]);
            auto lastMachine = get<0>(jobs[j].MachineRuntimes[i-1]);
            auto lastMachineProcessingTime = get<1>(jobs[j].MachineRuntimes[i-1]);
            model.addConstr(x[currentMachine][currentJob], GRB_GREATER_EQUAL, x[lastMachine][currentJob] + lastMachineProcessingTime);
//            cout << "Constraint added: " << "x[" << currentMachine << "]["  << currentJob << "] >= x[" << lastMachine << "][" << currentJob << "] + " << lastMachineProcessingTime << endl;
        }
    }
    // Constraint (4)/(5): x[i][j] >= x[i][k] + p[i][k] - V * z[i][j][k]
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < J; j++) {
            for (int k = 0; k < J; k++) {
                if (j < k) {
                    model.addConstr(x[i][j], GRB_GREATER_EQUAL, x[i][k] + p[i][k] - V * z[i][j][k]);
                    model.addConstr(x[i][k], GRB_GREATER_EQUAL, x[i][j] + p[i][j] - V * (1 - z[i][j][k]));
                }
            }
        }
    }

    GRBVar C = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER);
    GRBLinExpr objective = GRBLinExpr(C);
    for (int j = 0; j < J; j++) {
        auto lastMachineOfJob = get<0>(jobs[j].MachineRuntimes[jobs[j].MachineRuntimes.size()-1]);
        auto lastMachineProcessingTime = get<1>(jobs[j].MachineRuntimes[jobs[j].MachineRuntimes.size()-1]);
        model.addConstr(C, GRB_GREATER_EQUAL, x[lastMachineOfJob][j] + lastMachineProcessingTime);
    }

    model.setObjective(objective, GRB_MINIMIZE);
    model.set(GRB_DoubleParam_TimeLimit, timeLimit);
    model.update();
    auto startTime = std::chrono::high_resolution_clock::now();
    model.optimize();
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    bool optimal = false;
    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        optimal = true;
    }
    tuple<int, int, double, bool, double> result;
    result = make_tuple(model.get(GRB_DoubleAttr_ObjBound), model.get(GRB_DoubleAttr_ObjVal), duration.count(), optimal, model.get(GRB_DoubleAttr_MIPGap));

    return result;
}

tuple<int, int, double, bool, double> Solver::GurobiLiaoDisjunctiveSolver(double timeLimit) {
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);

    auto J = instance.NumberOfJobs;
    auto M = instance.NumberOfMachines;
    auto jobs = instance.Jobs;
    auto p = instance.ProcessingTime;

    int V = 0;

    // V = sum(pij)
    for (int i = 0; i < p.size(); i++) {
        for (int j = 0; j < p[i].size(); j++) {
            V += p[i][j];
        }
    }
    // Variable x
    auto** x = new GRBVar * [M];
    for (int i= 0; i < M; i++) {
        x[i] = reinterpret_cast<GRBVar *>(new GRBVar *[J]);
        for (int j = 0; j < J; j++) {
            x[i][j] = model.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER);
        }
    }

    auto***z = new GRBVar** [M];
    for (int i = 0; i < M; i++) {
        z[i] = reinterpret_cast<GRBVar **>(new GRBVar **[J]);
        for(int j = 0; j < J; j++) {
            z[i][j] = reinterpret_cast<GRBVar *>(new GRBVar *[J]);
            for (int k = 0; k < J; k++) {
                z[i][j][k] = model.addVar(0, 1, 0, GRB_INTEGER);
            }
        }
    }

    auto***q = new GRBVar** [M];
    for (int i = 0; i < M; i++) {
        q[i] = reinterpret_cast<GRBVar **>(new GRBVar **[J]);
        for(int j = 0; j < J; j++) {
            q[i][j] = reinterpret_cast<GRBVar *>(new GRBVar *[J]);
            for (int k = 0; k < J; k++) {
                q[i][j][k] = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER);
            }
        }
    }

    // Constraint (2): x[i][j] >= 0 forall j in J, i in M
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < J; j++) {
            model.addConstr(x[i][j], GRB_GREATER_EQUAL, 0);
        }
    }

    // Constraint (3): x[currentMachine, currentJob] >= x[lastMachine, currentJob] + processing_time[lastMachine, currentJob]
    for (int j = 0; j < jobs.size(); j++) {
        auto currentJob = jobs[j].Index;
        cout << "On job: " << currentJob << endl;
        for (int i = 1; i < M; i++) {
            auto currentMachine = get<0>(jobs[j].MachineRuntimes[i]);
            auto lastMachine = get<0>(jobs[j].MachineRuntimes[i-1]);
            auto lastMachineProcessingTime = get<1>(jobs[j].MachineRuntimes[i-1]);
            model.addConstr(x[currentMachine][currentJob], GRB_GREATER_EQUAL, x[lastMachine][currentJob] + lastMachineProcessingTime);
//            cout << "Constraint added: " << "x[" << currentMachine << "]["  << currentJob << "] >= x[" << lastMachine << "][" << currentJob << "] + " << lastMachineProcessingTime << endl;
        }
    }
    // Constraint (4)/(5): x[i][j] >= x[i][k] + p[i][k] - V * z[i][j][k]
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < J; j++) {
            for (int k = 0; k < J; k++) {
                if (j < k) {
                    model.addConstr(q[i][j][k], GRB_EQUAL, V * z[i][j][k] + (x[i][j] - x[i][k] - p[i][k]));
                    model.addConstr(q[i][j][k], GRB_LESS_EQUAL, V - p[i][j] - p[i][k]);
                }
            }
        }
    }

    GRBVar C = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER);
    GRBLinExpr objective = GRBLinExpr(C);
    for (int j = 0; j < J; j++) {
        auto lastMachineOfJob = get<0>(jobs[j].MachineRuntimes[jobs[j].MachineRuntimes.size()-1]);
        auto lastMachineProcessingTime = get<1>(jobs[j].MachineRuntimes[jobs[j].MachineRuntimes.size()-1]);
        model.addConstr(C, GRB_GREATER_EQUAL, x[lastMachineOfJob][j] + lastMachineProcessingTime);
    }

    model.set(GRB_DoubleParam_TimeLimit, timeLimit);
    model.setObjective(objective, GRB_MINIMIZE);
    model.update();
    auto startTime = std::chrono::high_resolution_clock::now();
    model.optimize();
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    bool optimal = false;
    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        optimal = true;
    }
    tuple<int, int, double, bool, double> result;
    result = make_tuple(model.get(GRB_DoubleAttr_ObjBound), model.get(GRB_DoubleAttr_ObjVal), duration.count()/1000, optimal, model.get(GRB_DoubleAttr_MIPGap));

    return result;
}

tuple<int, int, double, bool, double> Solver::GurobiTimeIndexedSolver(double timeLimit, int bigV, bool useV) {
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);
    env.set(GRB_IntParam_LogToConsole, 0);
    auto J = instance.NumberOfJobs;
    auto M = instance.NumberOfMachines;
    auto jobs = instance.Jobs;
    auto p = instance.ProcessingTime;

    int V = 0;
    if (useV) {
        V = bigV;
    } else {
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < J;j ++) {
                V += p[i][j];
            }
        }
    }


    auto*** x = new GRBVar **[M];
    for (int i = 0; i < M; i++) {
        x[i] = reinterpret_cast<GRBVar **>(new GRBVar **[J]);
        for(int j = 0; j < J; j++) {
            x[i][j] = reinterpret_cast<GRBVar *>(new GRBVar *[V]);
            for (int t = 0; t < V; t++) {
                x[i][j][t] = model.addVar(0, 1, 0, GRB_INTEGER);
            }
        }
    }


    // Constraint 11
    for (int i = 0; i < M; i++) {
        for (int j= 0; j< J; j++) {
            GRBLinExpr expr;
            for (int t = 0; t < V; t++) {
                expr += x[i][j][t];
            }
            model.addConstr(expr, GRB_EQUAL, 1);
            expr = 0;
        }
    }

    GRBVar C = model.addVar(0, V, 0, GRB_INTEGER);


//    model.set(GRB_DoubleAttr_ObjBound, double(lowerBound));

    GRBLinExpr objective = GRBLinExpr(C, GRB_MINIMIZE);

    // Constraint 12
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < J; j++) {
            GRBLinExpr sumT;
            for (int t = 0; t < V; t++) {
                sumT += (t + p[i][j]) * x[i][j][t];
            }
            model.addConstr(C, GRB_GREATER_EQUAL, sumT);
        }
    }

    // Constraint 13
    for (int i = 0; i < M; i++) {
        for (int t = 0; t < V; t++) {
            GRBLinExpr sumjt;
            for (int j = 0; j < J; j++) {
                for (int tAp = t - p[i][j] + 1; tAp <= t; tAp++) {
                    if (tAp < 0) {
                        continue;
                    } else {
                        sumjt += x[i][j][tAp];
                    }
                }
            }
            model.addConstr(sumjt, GRB_LESS_EQUAL, 1);
        }
    }

    for (int j = 0; j < J; j++) {
        for (int i  = 1; i < M; i++) {
            GRBLinExpr lhs;
            GRBLinExpr rhs;
            for (int t = 0; t < V; t++) {
                auto machAtHMinus1 = get<0>(jobs[j].MachineRuntimes[i-1]);
                lhs += (t + p[machAtHMinus1][j]) * x[machAtHMinus1][j][t];
                rhs += t * x[get<0>(jobs[j].MachineRuntimes[i])][j][t];
            }
            model.addConstr(lhs, GRB_LESS_EQUAL, rhs);
        }
    }

    model.setObjective(objective, GRB_MINIMIZE);
    model.set(GRB_DoubleParam_TimeLimit, timeLimit);
    model.update();
    auto startTime = std::chrono::high_resolution_clock::now();
    model.optimize();
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    bool optimal = false;
    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        optimal = true;
    }
    tuple<int, int, double, bool, double> result;
    result = make_tuple(model.get(GRB_DoubleAttr_ObjBound), model.get(GRB_DoubleAttr_ObjVal), duration.count(), optimal, model.get(GRB_DoubleAttr_MIPGap));

    return result;
}

tuple<int, int, double, bool, double> Solver::GurobiRankBasedSolver(double timeLimit) {
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);

    auto J = instance.NumberOfJobs;
    auto M = instance.NumberOfMachines;
    auto jobs = instance.Jobs;
    auto p = instance.ProcessingTime;

    int V = 0;

    // V = sum(pij)
    for (int i = 0; i < p.size(); i++) {
        for (int j = 0; j < p[i].size(); j++) {
            V += p[i][j];
        }
    }

    // Variable x[i][j][k]
    auto*** x = new GRBVar** [M];
    for (int i = 0; i < M; i++) {
        x[i] = reinterpret_cast<GRBVar **>(new GRBVar **[J]);
        for(int j = 0; j < J; j++) {
            x[i][j] = reinterpret_cast<GRBVar *>(new GRBVar *[J]);
            for (int k = 0; k < J; k++) {
                x[i][j][k] = model.addVar(0, 1, 0, GRB_INTEGER);
            }
        }
    }

    // Variable h[i][k]
    auto** h = new GRBVar* [M];
    for (int i = 0; i < M; i++) {
        h[i] = reinterpret_cast<GRBVar *>(new GRBVar* [J]);
        for (int k = 0; k < J; k++) {
            h[i][k] = model.addVar(0, V, 0, GRB_INTEGER);
        }
    }

    // Variable r[i][j][k]
    std::vector<std::vector<std::vector<int>>> r(M, std::vector<std::vector<int>>(J, std::vector<int>(M)));

    for (int j = 0; j < jobs.size(); j++) {
        auto jobIndex = j;
        auto machRuntimes = jobs[j].MachineRuntimes;
        for (int idx = 0; idx < machRuntimes.size(); idx++) {
            auto kOp = idx;
            auto machIdx = get<0>(machRuntimes[idx]);
            r[machIdx][j][kOp] = 1;
            cout << "On mach: " << machIdx << ", job: " << jobIndex << ", op: " << kOp << endl;
        }
    }

    // Constraint 17


    for (int i = 0; i < M; i++) {
        for (int k = 0; k < J; k++) {
            GRBLinExpr constr17;
            for (int j = 0; j < J; j++) {constr17 += x[i][j][k];}
            model.addConstr(constr17, GRB_EQUAL, 1);
        }
    }

    // Constraint 18
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < J; j++) {
            GRBLinExpr constr18;
            for (int k = 0; k < J; k++) {
                constr18 += x[i][j][k];
            }
            model.addConstr(constr18, GRB_EQUAL, 1);
        }
    }



    // Constraint 19
    for (int k = 0; k < J-1; k++) {
        for (int i = 0; i < M; i++) {
            GRBLinExpr sumJ;
            for (int j = 0; j < J; j++) {
                sumJ += p[i][j] * x[i][j][k];
            }
            model.addConstr(sumJ + h[i][k] , GRB_LESS_EQUAL, h[i][k+1]);
        }
    }


    // Constraint 20
//    for (int i = 0; i < M; i++) {
//        GRBLinExpr rhBase;
//        GRBLinExpr rp;
//        GRBLinExpr rxBase;
//        GRBLinExpr rxPlus;
//        GRBLinExpr rhPlus;
//        for (int j = 0; j < J; j++) {
//            for (int l = 0; l < M - 1; l++) {
//                for (int k = 0; k < J; k++) {
//                    rhBase += r[i][j][l] * h[i][k];
//                    rxBase += r[i][j][l] * x[i][j][k];
//                }
//                rp += r[i][j][l] * p[i][j];
//                // kAp: k with apostrophe k'
//                for (int kAp = 0; kAp < J; kAp++) {
//                    rxPlus += r[i][j][l+1] * x[i][j][kAp];
//                    rhPlus += r[i][j][l+1] * h[i][kAp];
//                }
//            }
//
//        }
//        model.addConstr(rhBase + rp , GRB_LESS_EQUAL, V * (1 - rxBase) + V * (1 - rxPlus) + rhPlus);
//        rhBase = 0;
//        rp = 0;
//        rxBase = 0;
//        rxPlus = 0;
//        rhPlus = 0;
//    }
    for (int j = 0; j < J; j++) {
        for (int k = 0; k < J; k++) {
            for (int kAp = 0; kAp < J; kAp++) {
                for (int l = 0; l < M-1; l++) {
                    GRBLinExpr sum1;
                    GRBLinExpr sum2;
                    GRBLinExpr sum3;
                    GRBLinExpr sum4;
                    GRBLinExpr sum5;
                    for (int i = 0; i < M; i++) {
                        sum1 += r[i][j][l] * h[i][k];
                        sum2 += r[i][j][l] * p[i][j];
                        sum3 += r[i][j][l] * x[i][j][k];
                        sum4 += r[i][j][l+1] * x[i][j][kAp];
                        sum5 += r[i][j][l+1] * h[i][kAp];
                    }
                    model.addConstr(sum1 + sum2, GRB_LESS_EQUAL, V * (1 - sum3) + V * (1 - sum4) + sum5);
                }
            }
        }
    }
    // Constraint 21
    GRBVar C = model.addVar(0, V, 0, GRB_INTEGER);
    GRBLinExpr objective = GRBLinExpr(C);
    for (int i = 0; i < M; i++) {
        GRBLinExpr sumPx;
        for (int j = 0; j < J; j++) {
            sumPx += p[i][j] * x[i][j][J-1];
        }
        model.addConstr(C, GRB_GREATER_EQUAL, h[i][J-1] + sumPx);
    }

//    for (int i = 0; i < M; i++) {
//        for (int k = 0; k < J; k++) {
//            model.addConstr(h[i][k], GRB_GREATER_EQUAL, 0);
//        }
//    }
//    for (int i = 0; i < M; i++) {
//        for (int j = 0; j < J; j++) {
//            for (int k = 0; k < J; k++) {
//                model.addConstr(0, GRB_LESS_EQUAL, x[i][j][k]);
//                model.addConstr(x[i][j][k], GRB_LESS_EQUAL, 1);
//            }
//        }
//    }

    model.setObjective(objective, GRB_MINIMIZE);
    model.set(GRB_DoubleParam_TimeLimit, timeLimit);
    model.update();
    auto startTime = std::chrono::high_resolution_clock::now();
    model.optimize();
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    bool optimal = false;
    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        optimal = true;
    }
    tuple<int, int, double, bool, double> result;
    result = make_tuple(model.get(GRB_DoubleAttr_ObjBound), model.get(GRB_DoubleAttr_ObjVal), duration.count()/1000, optimal, model.get(GRB_DoubleAttr_MIPGap));

    return result;
}



