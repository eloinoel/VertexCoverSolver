#ifndef VERTEXCOVERSOLVER_H
#define VERTEXCOVERSOLVER_H
#include <stdio.h>
#include <unordered_map>
#include <math.h>          // INFINITY
#include <chrono>

class BucketGraph;

/* Helper */

bool outOfTime(std::chrono::time_point<std::chrono::high_resolution_clock> startTime, double timeoutCap);

void resetGraphAfterBranching(BucketGraph* G, std::unordered_map<int, bool>* vc);

/* Exact Solver */

std::unordered_map<int, bool>* vcSolverRecursive(BucketGraph* G, int* numRec, bool printDebug);
std::unordered_map<int, bool>* vcSolverConstrained(BucketGraph* G, int* numRec, bool printDebug);


/* Heuristic Code */

std::unordered_map<int, bool>* maxHeuristicSolver(BucketGraph* G, int* numRec, bool applyReductions, bool randomise);

std::unordered_map<int, bool>* chooseSmallestHeuristicSolution(BucketGraph* G, int* numRec, std::unordered_map<int, bool>** currentSmallestVC, bool applyReductions = true,
 bool includeRandomsWithReductions = true, int numRandomSolverCalls = INT_MAX, int timeoutSoftCap = 20, bool printBestSolution = false);

std::unordered_map<int, bool>* fastVC(BucketGraph* G, std::unordered_map<int, bool>* vc, int* numRec, double timeout);

#endif