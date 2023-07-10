#include <stdio.h>
#include <iostream>
#include <fstream>  //ifstream file opening
#include <math.h>          // INFINITY
#include <chrono>
#include <random>

#include <signal.h>
#include <stdlib.h>

#include "utils/ColorPrint.h"
#include "utils/BucketGraph.h"
#include "VertexCoverSolver.h"
#include "utils/Benchmark.h"

using namespace std;
typedef ColorPrint cp;

/* SIG INT PRINTING */
BucketGraph* bucketGraph = nullptr;
bool printReductionDiff = false;
bool interrupted_by_sig = false;
bool printing_sol = false;
unordered_map<int, bool>* heuristicVC = nullptr;
int heuristicNumRecursions = 0;

void my_sig_handler(sig_atomic_t s)
{
    interrupted_by_sig = true;
    //print original edges if timeout in data reduction for (b) on exercise sheet
    if(printReductionDiff && bucketGraph != nullptr && !printing_sol)
    {
        printing_sol = true;
        std::vector<std::string>* str = bucketGraph->getOriginalEdgesToConsoleString();
        for(int i = 0; i < (int) str->size(); i++)
        {
            cout << str->at(i);
        }
        cout << "#difference: " << 0 << '\n';
        //free pointers
        //if(str) { delete str; }
        //if(bucketGraph) { bucketGraph->freeGraph(); }
    }

    //timeout at computing heuristic solution, output last best solution
    if(heuristicVC != nullptr && !printing_sol)
    {
        printing_sol = true;
        bucketGraph->printVertices(heuristicVC);
        cout << "#recursive steps: " << heuristicNumRecursions << endl;
        //free pointers
        //if(heuristicVC) { delete heuristicVC; }
        //if(bucketGraph) { bucketGraph->freeGraph(); }
    }

    exit(0);
}

/** Execute specific version of program with optional arguments for more prints
 * version
 * 0: Heuristic Solver
 * 1: Recursive Solver
 * 2: Constrained Solver

 * 5: (b) apply data reduction and output smaller graph and diff in vc size
 * 6: SAT Solver
 * ....
*/
void chooseImplementationAndOutput(int version = 1, bool printGraph = false, bool printMappings = false,
bool printDebug = false, bool printVCSize = false, bool printVC = true, bool printBounds = false)
{
    // Heuristic Solver
    if(version == 0)
    {
        auto startGraph = std::chrono::high_resolution_clock::now();
        bucketGraph = BucketGraph::readStandardInput();
        auto endGraph = std::chrono::high_resolution_clock::now();
        double graphConstructionDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endGraph - startGraph).count() /  1000) / (double) 1000;

        if (bucketGraph == nullptr)
            throw invalid_argument("Error constructing graph from input file.");
        if (printGraph)
            bucketGraph->print();

		if(printVC)
        {
            const int MAX_TIME_BUDGET = 60;
            const int INITIAL_SOLUTION_GENERATION_TIME_CAP = 40;//20; //in seconds
            const int HEURISTIC_SOLVER_TIME_CAP = INITIAL_SOLUTION_GENERATION_TIME_CAP - graphConstructionDuration;
            const int NUM_RANDOM_SOLUTION_GENERATIONS = 30;
            const int PRINT_TIME = 5;

            //std::cout << "before maxHeuristicSolver" << std::endl;
            heuristicNumRecursions = 0;
            auto startHeuristicWrapper = std::chrono::high_resolution_clock::now();
            //first generate a fast heuristic solution
            //cout << "before heuristic solver" << endl;
            heuristicVC = maxHeuristicSolver(bucketGraph, &heuristicNumRecursions, false, false);
            //see if we can find a better initial solution
            heuristicVC = chooseSmallestHeuristicSolution(bucketGraph, &heuristicNumRecursions, &heuristicVC, true, true, NUM_RANDOM_SOLUTION_GENERATIONS, HEURISTIC_SOLVER_TIME_CAP);
            auto endHeuristicWrapper = std::chrono::high_resolution_clock::now();
            double heuristicWrapperDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endHeuristicWrapper - startHeuristicWrapper).count() /  1000) / (double) 1000;
            //auto localSearchVC = fastVC(bucketGraph, heuristicVC, MAX_TIME_BUDGET);


            //set graph to state that it is in when vc vertices are inactive --> for fastVC() method
            //std::cout << "before setInactive" << std::endl;
            for(auto it = heuristicVC->begin(); it != heuristicVC->end(); ++it)
            {
                bucketGraph->setInactive(it->first);
            }
            double currentDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endHeuristicWrapper - startGraph).count() /  1000) / (double) 1000;
            //std::cout << "before fastVC" << std::endl;
            // TODO: find suitable timout value for localSearch
            const int LOCAL_SEARCH_TIME_CAP = MAX_TIME_BUDGET - currentDuration - 5;
            auto localSearchVC = fastVC(bucketGraph, heuristicVC, &heuristicNumRecursions, LOCAL_SEARCH_TIME_CAP);
            int localSearchVCSize = localSearchVC->size();
            int heuristicVCSize = heuristicVC->size();
            //int localSearchVCSize = heuristicVC->size();
            delete heuristicVC;
            heuristicVC = localSearchVC;

            auto startPrintSolution = std::chrono::high_resolution_clock::now();
            currentDuration = (std::chrono::duration_cast<std::chrono::microseconds>(startPrintSolution - startGraph).count() /  1000) / (double) 1000;

            //cout << "before print" << endl;
            if(!interrupted_by_sig)
            {
                //safely print solution, otherwise wait SIG
                if(currentDuration < MAX_TIME_BUDGET - PRINT_TIME)
                {
                    printing_sol = true;
                    bucketGraph->printVertices(heuristicVC);
                    cout << "#recursive steps: " << heuristicNumRecursions << endl;
                    //cout << "#recursive steps: " << currentDuration<< endl;
                    //cout << "#solution size: " << localSearchVCSize << endl;
                    //cout << "#recursive steps: " << 1000000 + localSearchVCSize - heuristicVCSize/* vcMax->size() - vcMaxRandom->size() */ << endl;
                }
                else
                {
                    while(true) {}
                }
            }

            auto endPrintSolution = std::chrono::high_resolution_clock::now();
            double printDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endPrintSolution - startPrintSolution).count() /  1000) / (double) 1000;
            //std::cout << "Total duration: " << graphConstructionDuration + heuristicWrapperDuration + printDuration << " seconds, Graph construction:" << graphConstructionDuration << " seconds, HeuristicWrapper: " << heuristicWrapperDuration << " seconds, Print solution: " << printDuration << "\n";
        }

        //free pointers
        //if(bucketGraph) { bucketGraph->freeGraph(); }
        //if(heuristicVC) { delete heuristicVC; heuristicVC = NULL; }
    }
    // Recursive Solver
    else if(version == 1)
    {
        auto startGraph = std::chrono::high_resolution_clock::now();
        BucketGraph* G = BucketGraph::readStandardInput();
        auto endGraph = std::chrono::high_resolution_clock::now();
        double graphConstructionDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endGraph - startGraph).count() /  1000) / (double) 1000;
        if(printDebug)
            std::cout << "#Constructed Graph of size n=" << G->getTotalNumVertices() << ", m=" << G->getNumEdges() << " in " << graphConstructionDuration << " seconds" << std::endl;
        if (G == nullptr)
            throw invalid_argument("Error constructing graph from input file.");
        if (printGraph)
        {
            G->print();
            //G->printActiveList();
            //G->printBucketQueue();
        }
        unordered_map<int, bool>* vc = nullptr;

        if(printVC)
        {
            int numRecursiveSteps = 0;
            vc = vcSolverRecursive(G, &numRecursiveSteps, printDebug);
            //writeSolutionToConsole(G->getStringsFromVertexIndices(vc));
            G->printVertices(vc);
            cout << "#recursive steps: " << numRecursiveSteps << endl;

            if(printVCSize)
            {
                cout << "#vc size: " << vc->size() << endl;
            }
        }

        if(printBounds)
        {
            int bound = G->getLowerBoundVC();
            cout << "#recursive steps: " << bound << endl;
        }

        //free pointers
        if(vc) { delete vc; }
        if(G) { delete G; }
    }
    // Constrained Solver
    else if(version == 2)
    {
        auto startGraph = std::chrono::high_resolution_clock::now();
        BucketGraph* G = BucketGraph::readStandardInput();
        auto endGraph = std::chrono::high_resolution_clock::now();
        double graphConstructionDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endGraph - startGraph).count() /  1000) / (double) 1000;
        if(printDebug)
            std::cout << "#Constructed Graph of size n=" << G->getTotalNumVertices() << ", m=" << G->getNumEdges() << " in " << graphConstructionDuration << " seconds" << std::endl;
        if (G == nullptr)
            throw invalid_argument("Error constructing graph from input file.");
        if (printGraph)
        {
            G->print();
            //G->printActiveList();
//            G->printBucketQueue();
//            G->printBucketSizes();
        }
        unordered_map<int, bool>* vc = nullptr;

        if(printVC)
        {
            int numRecursiveSteps = 0;
            vc = vcSolverConstrained(G, &numRecursiveSteps, printDebug);
            G->printVertices(vc);
            cout << "#recursive steps: " << numRecursiveSteps << endl;

            if(printVCSize)
            {
                cout << "#vc size: " << vc->size() << endl;
            }
        }

        if(printBounds)
        {
            int bound = G->getLowerBoundVC();
            cout << "#recursive steps: " << bound << endl;
        }

        //free pointers
        if(vc) { delete vc; }
        if(G) { delete G; }
    }
    else if(version == 5)
    {
        printReductionDiff = true;
        bucketGraph = BucketGraph::readStandardInput();
        if (bucketGraph == nullptr)
            throw invalid_argument("Error constructing graph from input file.");

        int numRecursiveSteps = 0;
        //unordered_map<int, bool>* vc = vcSolverRecursive(bucketGraph, &numRecursiveSteps);
        unordered_map<int, bool>* vc;

        //int tmpK = vc->size();
        int original_k = bucketGraph->getNumVertices();
        int new_k = original_k;
        bucketGraph->preprocess(&new_k);

        printing_sol = true;
        std::vector<std::string>* str = bucketGraph->getEdgesToConsoleString();

        if(!interrupted_by_sig)
        {
            for(int i = 0; i < (int) str->size(); i++)
            {
                cout << str->at(i);
            }
            cout << "#difference: " << to_string(original_k - new_k) << endl;
        }
        else
        {
            cout << "Select correct version please.";
        }
        bucketGraph->freeGraph();

        //free pointers
        //if(vc) { delete vc; }
        //if(bucketGraph) { bucketGraph->freeGraph(); }
    }
}

/*----------------------------------------------------------*/
/*-----------------------   Main   -------------------------*/
/*----------------------------------------------------------*/

int main(int argc, char* argv[]) {
    //By default, cin/cout waste time synchronizing themselves with the C library’s stdio buffers, so that you can freely intermix calls to scanf/printf with operations on cin/cout
    std::ios::sync_with_stdio(false);
	try
	{
        //TODO: disable printDebug for final submission
        chooseImplementationAndOutput(1, false, false, true, false, true, false);
//        chooseImplementationAndOutput(2, false, false, true, true, true, false);
        //chooseImplementationAndOutput(1, true, false, false, true, true, false); //print alot
        //chooseImplementationAndOutput(5, false, false, false, false, true, false);
    }
	catch (const exception& e)
	{
		cerr << ColorPrint::dye("Error while running vertex cover solver.\n", 'r');
        cerr << ColorPrint::dye(e.what() ,'r') << '\n';
	}
}