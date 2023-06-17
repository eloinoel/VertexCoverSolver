#include <iostream>
#include <fstream>  //ifstream file opening
#include <stack>          // std::stack
#include <math.h>          // INFINITY
#include "utils/ColorPrint.h"
#include "utils/BucketGraph.h"
#include "utils/SATSolver.h"
#include <chrono>
#include <random>


#include <signal.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;
typedef ColorPrint cp;

/*----------------------------------------------------------*/
/*---------------   Exercise 4 Solver Code   ---------------*/
/*----------------------------------------------------------*/

void resetGraphAfterBranching(BucketGraph* G, unordered_map<int, bool>* vc)
{
    for(auto it = vc->begin(); it != vc->end(); ++it)
    {
        G->setActive(it->first);
    }
}

unordered_map<int, bool>* maxHeuristicSolver(BucketGraph* G, int* numRec, bool applyReductions = false, bool randomise = false)
{
    // Apply Reduction Rules for the first time
    int numPreprocessingVCVertices = 0;
    if(applyReductions)
    {
        vector<bool> rulesToApply = vector<bool> ({true, true, true, false});
        G->preprocess(&numPreprocessingVCVertices, rulesToApply);
        numPreprocessingVCVertices = -numPreprocessingVCVertices;
    }

    unordered_map<int, bool>* vc = new unordered_map<int, bool>();
    while(true)
    {
        (*numRec)++;
        int vertex = -1;
        if(randomise)
        {
            vertex = G->getRandomMaxDegreeVertex(10);
        }
        else
        {
            vertex = G->getMaxDegreeVertex();
        }

        //no vertices left
        if (vertex == -1)
        {
            break;
        }

        int vertexDeg = G->getVertexDegree(vertex);

        //graph has no edges left
        if (vertexDeg == 0)
        {
            break;
        }

        //take vertex into vc
        G->setInactive(vertex);
        vc->insert({vertex, true});
    }
    
    resetGraphAfterBranching(G, vc);
    //undo reductions
    if(applyReductions)
    {
        int k = 0;
        G->unreduce(&k, vc->size()+numPreprocessingVCVertices, vc);
        G->resetLPBoundDataStructures();
    }

    return vc;
}

/* bool outOfTime(std::chrono::steady_clock::time_point startTime, int timeoutCap)
{
    auto currentTime = std::chrono::high_resolution_clock::now();
    double currentDuration = (std::chrono::duration_cast<std::chrono::microseconds>(currentTime - startTime).count() /  1000) / (double) 1000;
    return currentDuration > timeoutCap;
}

unordered_map<int, bool>* chooseSmallestHeuristicSolution(BucketGraph* G, int* numRec, bool applyReductions = false, int numRandomSolverCalls = 20, int timeoutSoftCap = 30)
{
    auto startTime = std::chrono::high_resolution_clock::now();

    unordered_map<int, bool>* vcMax = maxHeuristicSolver(G, numRec, false);

    if(outOfTime(startTime, timeoutSoftCap))
        return vcMax;

    //preprocessing on
    int vcMaxPreprocessingNumRec = 0;
    unordered_map<int, bool>* vcMaxPreprocessing = nullptr;
    if(applyReductions)
    {
        vcMaxPreprocessing = maxHeuristicSolver(G, &vcMaxPreprocessingNumRec, true);
        //select the best solution
        if(vcMaxPreprocessing->size() < vcMax->size())
        {
            delete vcMax;
            vcMax = vcMaxPreprocessing;
            *numRec = vcMaxPreprocessingNumRec;
        }
        else
        {
            delete vcMaxPreprocessing;
        }
    }

    if(outOfTime(startTime, timeoutSoftCap))
        return vcMax;


    //selecting random max degree vertices
    int vcMaxRandomNumRec;
    unordered_map<int, bool>* vcMaxRandom = nullptr;
    if(numRandomSolverCalls > 0)
    {
        for(int i = 0; i < numRandomSolverCalls - 1; ++i)
        {
            vcMaxRandomNumRec = 0;
            if(vcMaxRandom == nullptr)
            {
                vcMaxRandom = maxHeuristicSolver(G, &vcMaxRandomNumRec, true);
            }

            if(vcMaxRandom->size() < vcMax->size())
            {
                delete vcMax;
                vcMax = vcMaxRandom;
                *numRec = vcMaxRandomNumRec;
            }
            else
            {
                delete vcMaxRandom;
            }

            if(outOfTime(startTime, timeoutSoftCap))
                return vcMax;
        }
    }

    return vcMax;
} */

//https://ieeexplore.ieee.org/abstract/document/6486444
unordered_map<int, bool>* minHeuristicSolver(BucketGraph* G, int* numRec)
{
    unordered_map<int, bool>* vc = new unordered_map<int, bool>();
    while(true)
    {
        (*numRec)++;

        int vertex = G->getMinDegreeVertex();
        //cout << "looking at vertex " << vertex << "\n";
        //no vertices left or edges left
        if (vertex == -1)
        {
            break;
        }

        //take neighbours into vc
        vector<int>* neighbours = G->getNeighbours(vertex);
        for(int i = 0; i < (int) neighbours->size(); ++i)
        {
            G->setInactive(neighbours->at(i));
            vc->insert({neighbours->at(i), true});
        }
        delete neighbours;
    }

    int redundanceCount = 0;
    //scan the solution space for redundant vertices, this is probably dumb as fuck
    /* for(auto it = vc->begin(); it != vc->end(); ++it)
    {
        Vertex* v = G->getVertex(it->first);
        vector<int>* neighboursOfV = v->getAdj();
        bool redundantVertex = true;
        //all neighbours in vc as well
        for(int i = 0; i < (int) neighboursOfV->size(); i++)
        {
            int neighbourOfV = neighboursOfV->at(i);
            if(vc->find(neighbourOfV) == vc->end())
            {
                redundantVertex = false;
                break;
            }
        }
        if(redundantVertex)
        {
            vc->erase(it);
            ++redundanceCount;
            cout << cp::dye("minHeuristicSolver: found redundant vertex", 'r');
        }
    } */

    return vc;
}

void initGainLoss(BucketGraph* G, unordered_map<int, bool>* vc, std::vector<int>* gain, std::vector<int>* loss)
{
    for (int i=0; i<gain->size(); ++i)
    {
        (*gain)[i] = 0;
        (*loss)[i] = 0;
    }
    for (auto it = vc->begin(); it != vc->end(); ++it)
    {
        Vertex* v = G->getVertex(it->first);
        for (auto u = v->getAdj()->begin(); u != v->getAdj()->end(); ++u)
        {
            if(vc->find(*u) == vc->end())
            {
                (*loss)[it->first]++;
            }
        }
    }
}

int getMinLossIndex(unordered_map<int, bool>* vc, std::vector<int>* loss)
{
    int min = INT32_MAX;
    int minIndex = -1;
    for (int i=0; i<vc->size(); ++i)
    {
        if (loss->at(i) < min)
        {
            min = loss->at(i);
            minIndex = i;
            continue;
        }
    }
    return minIndex;
}

void removeMinLossVCVertex(BucketGraph* G, unordered_map<int, bool>* vc, std::vector<int>* gain, std::vector<int>* loss)
{
    int u_index = getMinLossIndex(vc, loss);
    auto u = vc->find(u_index);
    if (u == vc->end()) { throw invalid_argument("Iterator of u not found"); }
    vc->erase(u);
    (*gain)[u_index] = (*loss)[u_index];
    (*loss)[u_index] = 0;
    Vertex* uVertex = G->getVertex(u_index);
    for (auto neighbour = uVertex->getAdj()->begin(); neighbour != uVertex->getAdj()->end(); ++neighbour)
    {
        (*gain)[*neighbour]++;
    }
}

void addRandomUncoveredEdgeMaxGainEndpointVertex(BucketGraph* G, unordered_map<int, bool>* vc, std::vector<int>* gain, std::vector<int>* loss)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0, G->getBucketQueue()->size());
    int rand = (int) distr(gen) + (G->getNumVertices() - G->getNumConnectedVertices());

    // TODO:
}

/*
*   timeout in milliseconds
*/
unordered_map<int, bool>* fastVC(BucketGraph* G, unordered_map<int, bool>* vc, int timeout)
{
    int n = G->getTotalNumVertices();
    std::vector<int> gain = std::vector<int>(n);
    std::vector<int> loss = std::vector<int>(n);
    unordered_map<int, bool> bestVC (*vc);
    unordered_map<int, bool> currentVC (*vc);
    initGainLoss(G, &currentVC, &gain, &loss);
    auto startFastVC = std::chrono::high_resolution_clock::now();
    while (true)
    {
        auto currentFastVC = std::chrono::high_resolution_clock::now();
        int fastVC = std::chrono::duration_cast<std::chrono::microseconds>(currentFastVC - startFastVC).count() / 1000;
        if (fastVC >= timeout) { break; }
        // we found an improved vc! time to overwrite our previous best solution and search further
        if (G->getMaxDegree() <= 0)
        {
            // TODO: check that this is overridden by VALUE!
            bestVC = currentVC;
            removeMinLossVCVertex(G, &currentVC, &gain, &loss);
        }
        // TODO: replace removeMinLossVCVertex with BMS function, only checking a set number of vertices and taking that minimum
        removeMinLossVCVertex(G, &currentVC, &gain, &loss);


    }
    return &bestVC;
}

/*----------------------------------------------------------*/
/*---------------   Exercise 3 Solver Code   ---------------*/
/*----------------------------------------------------------*/

unordered_map<int, bool>* vcVertexBranchingRecursive(BucketGraph* G, int k, int depth, int* numRec)
{
    (*numRec)++;
	if (k < 0)
    {
		return nullptr;
	}
    int previousK = k;
    bool cut = false;
    //if(depth /* % 10 */ == 0) { cut = G->reduce(&k); }
    //std::cout << "> cutting through data reductions " << '\n';
    cut = G->reduce(&k);
    if(cut)
    {
        //std::cout << "> cutting through data reductions " << '\n';
        G->unreduce(&k, previousK);
        return nullptr;
    }

    //std::cout << "> calculated LPBound: " << G->getLPBound() << " with k=" << k << '\n';
    if (k < G->getLPBound()) {
        G->unreduce(&k, previousK);
        return nullptr;
    }

    //cout << "before getMaxDegreeVertex" << endl;
	int vertex = G->getMaxDegreeVertex();
    //no vertices left
    if (vertex == -1)
    {
        unordered_map<int, bool>* vc = new unordered_map<int, bool>();
        G->unreduce(&k, previousK, vc);
        return vc;
    }
    //cout << "before getVertexDegree: " << vertex << endl;
    int vertexDeg = G->getVertexDegree(vertex);
    //cout << "got maxDegree vertex: " << vertex << endl;
    //G->print();
	//graph has no edges left
	if (vertexDeg == 0)
	{
		unordered_map<int, bool>* vc = new unordered_map<int, bool>();
        G->unreduce(&k, previousK, vc);
        return vc;
	}

    //cout << cp::dye("branching: choosing vertex: " + std::to_string(vertex), 'b') << endl;
	//delete first vertex from graph and explore solution
    G->setInactive(vertex);
    //cout << "before branching" << endl;
	unordered_map<int, bool>* S = vcVertexBranchingRecursive(G, k - 1, depth+1, numRec);
	if (S != nullptr)
	{
        //revert changes for multiple executions of the algorithm
        G->setActive(vertex);
        S->insert({vertex, true}); //push results
        G->unreduce(&k, previousK, S); //unreduce needs correct vc für unmerge of deg2rule
		//S->push_back(vertex);
		return S;
	}
	else
	{
        //cout << "before setActive" << endl;
		//revert changes to graph
		G->setActive(vertex);
        G->unreduce(&k, previousK);
	}
    //cout << cp::dye("restoring vertex: ", 'g') << vertex << endl;

	//cannot fully explore neighbours
    if (vertexDeg > k)
    {
        //G->unreduce(&k, previousK);
        return nullptr;
    }

    //cout << "deleting neighbourhood of vertex " << vertex << ": ";
    vector<int>* neighbours = G->getNeighbours(vertex);
    /* cout << ColorPrint::dye("branching: choosing neighbours of vertex " + std::to_string(vertex) + ": ", 'b');
    for(int i = 0; i < (int) neighbours->size(); i++)
    {
        cout <<  ColorPrint::dye(std::to_string(neighbours->at(i)) + ", ", 'b');
    }
    cout << endl; */
    //if(vertex == 24) { G->print(); }
    G->setInactive(neighbours);
    //cout << "pinc " << '\n';
	S = vcVertexBranchingRecursive(G, k - neighbours->size(), depth+1, numRec);
    //cout << "prec " << '\n';
	if (S != nullptr)
	{
        //revert changes for multiple executions of the algorithm
        G->setActive(neighbours);
        for (int i = 0; i < (int) neighbours->size(); i++)
        {
            //S->push_back(neighbours->at(i));
            S->insert({neighbours->at(i), true}); //push results
        }
        G->unreduce(&k, previousK, S); //unreduce needs correct vc für unmerge of deg2rule
        return S;
	}
	else
	{
		//revert changes to graph
		G->setActive(neighbours);
        G->unreduce(&k, previousK);
	}
    /* cout << "restoring neighbourhood of vertex " << vertex << ": ";
    for(int i = 0; i < (int) neighbours->size(); i++)
    {
        cout << neighbours->at(i) << ", ";
    }
    cout << endl; */
    // free neighbours
    delete neighbours;

    G->unreduce(&k, previousK);
    return nullptr;
}

unordered_map<int, bool>* vcSolverRecursive(BucketGraph* G, int* numRec)
{
    int numPreprocessingVCVertices = 0;
	int k = 0;
    // Apply Reduction Rules for the first time
    G->preprocess(&numPreprocessingVCVertices);
    numPreprocessingVCVertices = -numPreprocessingVCVertices;
    k = G->getLowerBoundVC();
    //G->printBucketSizes();

//    G->printBucketQueue();
//    G->print();

	unordered_map<int, bool>* vc;

	while (true)
	{

        // Begin Branching
		vc = vcVertexBranchingRecursive(G, k, 0, numRec);
		if (vc != nullptr)
		{
            // Add Reduced Vertices to Vertex Cover
            G->unreduce(&k, vc->size()+numPreprocessingVCVertices, vc);
			return vc;
		}
        k++;
	}
}

/*----------------------------------------------------------*/
/*-----------------   Helper Functions   -------------------*/
/*----------------------------------------------------------*/

void writeSolutionToFile(string fileName, vector<string>* vc)
{
	ofstream outfile(fileName);
	for (auto it = vc->begin(); it != vc->end(); ++it)
	{
		outfile << *it << '\n';
	}
	outfile.close();
}

void writeSolutionToConsole(vector<string>* vc)
{
	for (auto it = vc->begin(); it != vc->end(); ++it)
	{
		cout << *it << '\n';
	}
}

/*----------------------------------------------------------*/
/*---------------------   Execution   ----------------------*/
/*----------------------------------------------------------*/


BucketGraph* bucketGraph = nullptr;
bool interrupted_by_sig = false;
bool printing_sol = false;

void my_sig_handler(sig_atomic_t s)
{
    //print original edges if timeout in data reduction for (b) on exercise sheet
    if(bucketGraph != nullptr && !printing_sol)
    {
        std::vector<std::string>* str = bucketGraph->getOriginalEdgesToConsoleString();
        for(int i = 0; i < (int) str->size(); i++)
        {
            cout << str->at(i);
        }
        cout << "#difference: " << 0 << '\n';
        interrupted_by_sig = true;
    }
    exit(0);
}

/** Execute specific version of program with optional arguments for more prints
 * version
 * 0: Heuristic Solver
 * 1: Bucketgraph

 * 5: (b) apply data reduction and output smaller graph and diff in vc size
 * 6: SAT Solver
 * ....
*/
void chooseImplementationAndOutput(int version = 1, bool printGraph = false, bool printMappings = false, 
bool printDebug = false, bool printVCSize = false, bool printVC = true, bool printBounds = false)
{
    if(version == 0)
    {
        auto startGraph = std::chrono::high_resolution_clock::now();
        BucketGraph* G = BucketGraph::readStandardInput();
        auto endGraph = std::chrono::high_resolution_clock::now();
        if (G == nullptr)
            throw invalid_argument("Error constructing graph from input file.");
        if (printGraph)
            G->print();

		if(printVC)
        {
            int numRecursions = 0;
            auto startHeuristicWrapper = std::chrono::high_resolution_clock::now();
            unordered_map<int, bool>* vc = maxHeuristicSolver(G, &numRecursions, false);
            //unordered_map<int, bool>* vc = chooseSmallestHeuristicSolution(G, &numRecursions, true, 10); //TODO: doesn't compile yet
            auto endHeuristicWrapper = std::chrono::high_resolution_clock::now();

            G->printVertices(vc);

            cout << "#recursive steps: " << numRecursions << endl;
            double heuristicWrapperDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endHeuristicWrapper - startHeuristicWrapper).count() /  1000) / (double) 1000;
            //std::cout << "HeuristicWrapper: " << heuristicWrapperDuration << " seconds, \n"; 

            
            /* 
            auto startHeuristic = std::chrono::high_resolution_clock::now();
            unordered_map<int, bool>* vcMax = maxHeuristicSolver(G, &numRecursions, false);
            auto endHeuristic = std::chrono::high_resolution_clock::now();

            auto startPreprocess = std::chrono::high_resolution_clock::now();
            unordered_map<int, bool>* vcMaxPreprocess = maxHeuristicSolver(G, &numRecursions, true);
            auto endPreprocess = std::chrono::high_resolution_clock::now();

            auto startRandomHeuristic = std::chrono::high_resolution_clock::now();
            unordered_map<int, bool>* vcMaxRandom = maxHeuristicSolver(G, &numRecursions, false, true);
            auto endRandomHeuristic = std::chrono::high_resolution_clock::now(); */


            //unordered_map<int, bool>* vcMin = minHeuristicSolver(G, &numRecursions);
            //vector<int>* vc = heuristicSolver(G, &numRecursions);
            
            //unordered_map<int, bool>* bestVC = vcMax;
            /* if(vcMaxNoPreprocess->size() < bestVC->size())
            {
                bestVC = vcMaxNoPreprocess;
            } */
            /* if(vcMaxRandom->size() < bestVC->size())
            {
                bestVC = vcMaxRandom;
            } */
            //cout << "before print vertices" << endl;
            /* auto startPrintSolution = std::chrono::high_resolution_clock::now();
            G->printVertices(bestVC);
            auto endPrintSolution = std::chrono::high_resolution_clock::now();

            double Graph = (std::chrono::duration_cast<std::chrono::microseconds>(endGraph - startGraph).count() /  1000) / (double) 1000;
            double Heur = (std::chrono::duration_cast<std::chrono::microseconds>(endHeuristic - startHeuristic).count() /  1000) / (double) 1000;
            double HeurRandom = (std::chrono::duration_cast<std::chrono::microseconds>(endRandomHeuristic - startRandomHeuristic).count() /  1000) / (double) 1000;
            double HeurPreprocess = (std::chrono::duration_cast<std::chrono::microseconds>(endPreprocess - startPreprocess).count() /  1000) / (double) 1000;
            double Print = (std::chrono::duration_cast<std::chrono::microseconds>(endPrintSolution - startPrintSolution).count() /  1000) / (double) 1000;
            std::cout << "Computed heuristic solution in " << Graph + Heur + HeurRandom + HeurPreprocess + Print << " seconds (Graph construction: " << Graph << " + Heuristic solving: " << Heur + HeurRandom + HeurPreprocess<< " + Printing solution: " << Print << ")" << '\n';
            std::cout << "Heuristic: " << Heur << " seconds, Preprocess + Heuristic: " << HeurPreprocess << " seconds, Random Heuristic: " << HeurRandom << " seconds" << '\n'; */
            //double HeurInt = (std::chrono::duration_cast<std::chrono::microseconds>(endHeuristic - startHeuristic).count() /  1000);
            //double GraphInt = (std::chrono::duration_cast<std::chrono::microseconds>(endGraph - startGraph).count() /  1000);
            //double PrintInt = (std::chrono::duration_cast<std::chrono::microseconds>(endPrintSolution - startPrintSolution).count() /  1000);
            //cout << "#recursive steps: " << /* GraphInt + HeurInt +  */PrintInt << endl;
            //cout << "#recursive steps: " << 1000000 + vcMax->size() - vcMaxNoPreprocess->size() << endl;
            //cout << "#recursive steps: " << 1000000 + vcMax->size() - vcMaxRandom->size() << endl;
            //cout << vcMax->size() << endl;
            //cout << vcMaxNoPreprocess->size() << endl;
            //cout << "#recursive steps: " << numRecursions << endl;
            /* if (printVCSize)
                cout << "VC size: " << vcMax->size() << endl; */
        }
    }
    else if(version == 1)
    {
        BucketGraph* G = BucketGraph::readStandardInput();
        if (G == nullptr)
            throw invalid_argument("Error constructing graph from input file.");
        if (printGraph)
        {
            G->print();
            //G->printActiveList();
            //G->printBucketQueue();
        }

        if(printVC)
        {
            int numRecursiveSteps = 0;
            unordered_map<int, bool>* vc = vcSolverRecursive(G, &numRecursiveSteps);
            //writeSolutionToConsole(G->getStringsFromVertexIndices(vc));
            G->printVertices(vc);
            cout << "#recursive steps: " << numRecursiveSteps << endl;

            if(printVCSize)
            {
                cout << "vc size: " << vc->size() << endl;
            }
        }

        if(printBounds)
        {
            int bound = G->getLowerBoundVC();
            cout << "#recursive steps: " << bound << endl;
        }

    }
    else if(version == 5)
    {
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
    }
    else if(version == 6)
    {
        // For now it all runs from its constructor function
        SATSolver SATSolver;
    }
    else
    {
        cout << "Select correct version please.";
    }
}

/*----------------------------------------------------------*/
/*-----------------------   Main   -------------------------*/
/*----------------------------------------------------------*/

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false); //By default, cin/cout waste time synchronizing themselves with the C library’s stdio buffers, so that you can freely intermix calls to scanf/printf with operations on cin/cout
	try
	{
        chooseImplementationAndOutput(1, false, false, false, false, true, false);
        //chooseImplementationAndOutput(1, true, false, false, true, true, false); //print alot
        //signal(SIGINT, my_sig_handler); //catches SIGINT to output anything
        //chooseImplementationAndOutput(5, false, false, false, false, true, false);
	}
	catch (const exception& e)
	{
		cerr << ColorPrint::dye("Error while running vertex cover solver.\n", 'r');
        cerr << ColorPrint::dye(e.what() ,'r');
	}
}