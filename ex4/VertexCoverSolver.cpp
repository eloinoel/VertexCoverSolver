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

bool outOfTime(std::chrono::time_point<std::chrono::high_resolution_clock> startTime, double timeoutCap)
{
    auto currentTime = std::chrono::high_resolution_clock::now();
    double currentDuration = std::chrono::duration_cast<std::chrono::microseconds>(currentTime - startTime).count() / (double) 1000000;
    return currentDuration > timeoutCap;
}

/*
* will try generating multiple solutions to take the best one
* will stop if current time elapsed given @timeoutSoftCap
* call with nullptr as @currentSmallestVC if want to generate simple and fast max heuristic solution at beginning
//TODO: possible optimisation to only calculate preprocessing once at beginning and copy graph and then just add reduction vertices to solutions
*/
unordered_map<int, bool>* chooseSmallestHeuristicSolution(BucketGraph* G, int* numRec, unordered_map<int, bool>** currentSmallestVC, bool applyReductions = true, bool includeRandomsWithReductions = true, int numRandomSolverCalls = 20, int timeoutSoftCap = 20)
{
    int best_solution = 0; //0: maxHeuristic, 1: with preprocessing, 2: randomised, 3: randomised with preprocessing
    auto startTime = std::chrono::high_resolution_clock::now();

    if(*currentSmallestVC == nullptr)
    {
        *currentSmallestVC = maxHeuristicSolver(G, numRec);
    }

    //cout << "maxHeuristicSolver size: " << (*currentSmallestVC)->size() << ", addr: " << *currentSmallestVC << endl;

    if(outOfTime(startTime, timeoutSoftCap)) { return *currentSmallestVC; }

    //preprocessing on
    int vcMaxPreprocessingNumRec = 0;
    unordered_map<int, bool>* vcMaxPreprocessing = nullptr;
    if(applyReductions)
    {
        vcMaxPreprocessing = maxHeuristicSolver(G, &vcMaxPreprocessingNumRec, true);
        //select the best solution
        if(vcMaxPreprocessing->size() < (*currentSmallestVC)->size())
        {
            delete *currentSmallestVC;
            *currentSmallestVC = vcMaxPreprocessing;
            *numRec = vcMaxPreprocessingNumRec;
            best_solution = 1;
        }
        else
        {
            delete vcMaxPreprocessing;
        }
    }

    //cout << "after Preprocessing size: " << (*currentSmallestVC)->size() << ", addr: " << *currentSmallestVC << endl;

    if(outOfTime(startTime, timeoutSoftCap)) { return *currentSmallestVC; }

    //selecting random max degree vertices
    int vcMaxRandomNumRec;
    unordered_map<int, bool>* vcMaxRandom = nullptr;
    if(numRandomSolverCalls > 0)
    {
        for(int i = 0; i < numRandomSolverCalls; ++i)
        {
            vcMaxRandomNumRec = 0;
            bool usedPreprocessing = false;
            if(includeRandomsWithReductions && (i % 2 == 1)) //every second random iteration also uses preprocessing
            {
                vcMaxRandom = maxHeuristicSolver(G, &vcMaxRandomNumRec, true, true);
                usedPreprocessing = true;
            }
            else
            {
                vcMaxRandom = maxHeuristicSolver(G, &vcMaxRandomNumRec, false, true);
            }

            if(vcMaxRandom->size() < (*currentSmallestVC)->size())
            {
                delete *currentSmallestVC;
                *currentSmallestVC = vcMaxRandom;
                *numRec = vcMaxRandomNumRec;
                if(usedPreprocessing)
                {
                    best_solution = 3;
                }
                else
                {
                    best_solution = 2;
                }
            }
            else
            {
                delete vcMaxRandom;
            }

            if(outOfTime(startTime, timeoutSoftCap)) { return *currentSmallestVC; }
        }
    }
    //cout << "after random size: " << (*currentSmallestVC)->size() << ", addr: " << *currentSmallestVC << endl;
    //cout << "#recursive steps: " << best_solution << endl;

    return *currentSmallestVC;
}

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

void initGainLossAdditionTime(BucketGraph* G, unordered_map<int, bool>* vc, std::vector<int>* gain, std::vector<int>* loss, std::vector<int>* additionTime)
{
    for (int i = 0; i < (int) gain->size(); ++i)
    {
        (*gain)[i] = 0;
        (*loss)[i] = 0;
        (*additionTime)[i] = 0;
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
    if (vc == nullptr) { throw invalid_argument("getMinLossIndex: passed vc is nullptr"); }
    if (vc->empty()) { throw invalid_argument("getMinLossIndex: passed vc is empty"); }

    int min = INT32_MAX;
    int minIndex = -1;

    for (auto it = vc->begin(); it != vc->end(); ++it)
    {
        if (loss->at(it->first) < min)
        {
            min = loss->at(it->first);
            minIndex = it->first;
            continue;
        }
    }
    return minIndex;
}

int getBMSMinLossIndex(unordered_map<int, bool>* vc, std::vector<int>* loss, int k)
{
    if (vc == nullptr) { throw invalid_argument("getMinLossIndex: passed vc is nullptr\n"); }
    if (vc->empty()) { throw invalid_argument("getMinLossIndex: passed vc is empty\n"); }
    if (k > (int) vc->size()) { throw invalid_argument("getMinLossIndex: passed k is larger than passed vc\n"); }
    if (k < 1) { throw invalid_argument("getMinLossIndex: passed k that smaller than 1\n"); }

    int min = INT32_MAX;
    int minIndex = -1;

    std::random_device rd;
    std::mt19937 gen(rd());
    int partitionSize = vc->size() / k;
    //std::uniform_int_distribution<> distr(0, vc->size()-1);
    std::uniform_int_distribution<> distr(0, partitionSize-1);

    for (int i = 0; i < k; ++i)
    {

    }
    int i=0;
    for (auto it = vc->begin(); it != vc->end() && i<k; ++i, ++it)
    {
        int rndIndex = (int) distr(gen);
        if((i+1)*partitionSize > (int) vc->size())
        {
            rndIndex = rndIndex % ((i+1)*partitionSize - vc->size());
        }

        for(int j=0; j<rndIndex; ++j)
        {
            ++it;
        }

        if (loss->at(it->first) < min)
        {
            min = loss->at(it->first);
            minIndex = it->first;
            continue;
        }
    }
    return minIndex;
}

void removeMinLossVCVertex(BucketGraph* G, unordered_map<int, bool>* vc, std::vector<int>* gain, std::vector<int>* loss)
{
    int u_index = getMinLossIndex(vc, loss);
    auto u = vc->find(u_index);
    /* std::cout << "uIndex: " << u_index << std::endl;
    std::cout << "vc size: " << vc->size() << std::endl;
    int i=0;
    for (auto it = vc->begin(); it != vc->end(); ++it)
    {
        std::cout << "vc[" << i << "]: " << it->first << std::endl;
        i++;
    } */
    if (u == vc->end()) { throw invalid_argument("removeMinLossVCVertex: Iterator of u not found\n"); }
    vc->erase(u);
    G->setActive(u_index);
    (*gain)[u_index] = (*loss)[u_index];
    (*loss)[u_index] = 0;
    Vertex* uVertex = G->getVertex(u_index);
    for (auto neighbour = uVertex->getAdj()->begin(); neighbour != uVertex->getAdj()->end(); ++neighbour)
    {
        if (!(G->getVertex(*neighbour)->getActive())) { continue; }
        (*gain)[*neighbour]++;
    }
    //std::cout << cp::dye("Removed " + std::to_string(u_index), 'r') << " from vc" << std::endl;
}

void removeBMSMinLossVCVertex(BucketGraph* G, unordered_map<int, bool>* vc, std::vector<int>* gain, std::vector<int>* loss, int k)
{
    int u_index = getBMSMinLossIndex(vc, loss, k);
    auto u = vc->find(u_index);

    if (u == vc->end()) { throw invalid_argument("removeMinLossVCVertex: Iterator of u not found\n"); }
    vc->erase(u);
    G->setActive(u_index);
    (*gain)[u_index] = (*loss)[u_index];
    (*loss)[u_index] = 0;
    Vertex* uVertex = G->getVertex(u_index);
    for (auto neighbour = uVertex->getAdj()->begin(); neighbour != uVertex->getAdj()->end(); ++neighbour)
    {
        if (!(G->getVertex(*neighbour)->getActive())) { continue; }
        (*gain)[*neighbour]++;
    }
}

void addRandomUncoveredEdgeMaxGainEndpointVertex(BucketGraph* G, unordered_map<int, bool>* vc, std::vector<int>* gain, std::vector<int>* loss, std::vector<int>* additionTime, int numRecursions)
{
    //std::cout << "Adding to vc..." << std::endl;
    int v = G->getRandomConnectedVertex(10);
    if (v == -1) { return; }
    Vertex* vertex = G->getVertex(v);
    if (vertex->getDegree() == 0) { throw invalid_argument("addRandomUncoveredEdgeMaxGainEndpointVertex: Got vertex of degree 0\n"); }
    //std::cout << cp::dye("Chose edge adjacent to vertex: " + std::to_string(v), 'y') << std::endl;

    auto neighbours = G->getNeighbours(v);  // TODO: slow getNeighbours
    //std::cout << "neighbours size: " << neighbours->size() << std::endl;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0, neighbours->size()-1);
    int rnd = (int) distr(gen);
    //std::cout << cp::dye("Chose vertex neighbours index: " + std::to_string(rnd), 'y') << std::endl;
    int neighbour = neighbours->at(rnd);
    Vertex* neighbourVertex = G->getVertex(neighbour);
    //std::cout << cp::dye("Chose vertex neighbour: " + std::to_string(neighbour), 'y') << std::endl;

    // add vertex with higher gain into vc
    if((*gain)[v] > (*gain)[neighbour] || ((*gain)[v] == (*gain)[neighbour] && (*additionTime)[v] <= (*additionTime)[neighbour]))
    {
        for (auto u = vertex->getAdj()->begin(); u != vertex->getAdj()->end(); ++u)
        {
            if (!(G->getVertex(*u)->getActive())) { continue; }
            (*gain)[*u]--;
        }
        G->setInactive(v);
        vc->insert({v, true});
        (*gain)[v] = 0;
        (*additionTime)[v] = numRecursions;
        //std::cout << cp::dye("Added " + std::to_string(v), 'g') << " to vc" << std::endl;
    }
    else
    {
        for (auto u = neighbourVertex->getAdj()->begin(); u != neighbourVertex->getAdj()->end(); ++u)
        {
            if (!(G->getVertex(*u)->getActive())) { continue; }
            (*gain)[*u]--;
        }
        G->setInactive(neighbour);
        vc->insert({neighbour, true});
        (*gain)[neighbour] = 0;
        (*additionTime)[neighbour] = numRecursions;
        //std::cout << cp::dye("Added " + std::to_string(neighbour), 'g') << " to vc" << std::endl;
    }
}

/*
*   vertices of vertex cover vc need to be inactive in graph G
*   timeout in seconds
*   vc is not changed, you need to update it manually after the function call
*/
unordered_map<int, bool>* fastVC(BucketGraph* G, unordered_map<int, bool>* vc, double timeout)
{
    //std::cout << "before guards" << std::endl;
    if (vc == nullptr) { throw invalid_argument("fastVC: passed vc is nullptr\n"); }
    for (auto it = vc->begin(); it != vc->end(); ++it)
    {
        if (G->getVertex(it->first)->getActive())
        {
            throw invalid_argument("fastVC: passed vc vertices must be set to inactive in graph\n");
        }
    }

    //std::cout << "before inits" << std::endl;
    int numRecursions = 0;
    int n = G->getTotalNumVertices();
    std::vector<int> gain = std::vector<int>(n);
    std::vector<int> loss = std::vector<int>(n);
    std::vector<int> additionTime = std::vector<int>(n);
    unordered_map<int, bool>* currentVC = new unordered_map<int, bool>(*vc);
    unordered_map<int, bool>* bestVC = new unordered_map<int, bool>(*vc);
    initGainLossAdditionTime(G, currentVC, &gain, &loss, &additionTime);
    auto startFastVC = std::chrono::high_resolution_clock::now();
    while (true && vc->size() > 1)
    {
        if (outOfTime(startFastVC, timeout)) { break; }
        numRecursions++;
        // we found an improved vc! time to overwrite our previous best solution and search further
        if (G->getMaxDegree() <= 0)
        {
            //std::cout << cp::dye("Found vc of size: " + std::to_string(currentVC->size()), 'g') << std::endl;
            delete bestVC;
            bestVC = new unordered_map<int, bool>(*currentVC);
            removeMinLossVCVertex(G, currentVC, &gain, &loss);
        }
        // TODO: test different k's and find suitable k value
        removeBMSMinLossVCVertex(G, currentVC, &gain, &loss, std::min(100, (int) currentVC->size())/* std::max(1, (int) currentVC->size() / 4) */);
        addRandomUncoveredEdgeMaxGainEndpointVertex(G, currentVC, &gain, &loss, &additionTime, numRecursions);
    }
    delete currentVC;
    return bestVC;
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

/* SIG INT PRINTING */
BucketGraph* bucketGraph = nullptr;
bool interrupted_by_sig = false;
bool printing_sol = false;
unordered_map<int, bool>* heuristicVC = nullptr;
int heuristicNumRecursions = 0;

void my_sig_handler(sig_atomic_t s)
{
    interrupted_by_sig = true;
    //print original edges if timeout in data reduction for (b) on exercise sheet
    if(bucketGraph != nullptr && !printing_sol)
    {
        printing_sol = true;
        std::vector<std::string>* str = bucketGraph->getOriginalEdgesToConsoleString();
        for(int i = 0; i < (int) str->size(); i++)
        {
            cout << str->at(i);
        }
        cout << "#difference: " << 0 << '\n';
    }

    //timeout at computing heuristic solution, output last best solution
    if(heuristicVC != nullptr && !printing_sol)
    {
        printing_sol = true;
        bucketGraph->printVertices(heuristicVC);
        cout << "#recursive steps: " << heuristicNumRecursions << endl;
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
            const int INITIAL_SOLUTION_GENERATION_TIME_CAP = 20; //in seconds
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
            //heuristicVC = chooseSmallestHeuristicSolution(bucketGraph, &heuristicNumRecursions, &heuristicVC, true, true, NUM_RANDOM_SOLUTION_GENERATIONS, HEURISTIC_SOLVER_TIME_CAP);
            chooseSmallestHeuristicSolution(bucketGraph, &heuristicNumRecursions, &heuristicVC, true, true, NUM_RANDOM_SOLUTION_GENERATIONS, HEURISTIC_SOLVER_TIME_CAP);
            auto endHeuristicWrapper = std::chrono::high_resolution_clock::now();
            double heuristicWrapperDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endHeuristicWrapper - startHeuristicWrapper).count() /  1000) / (double) 1000;

            //auto localSearchVC = fastVC(bucketGraph, heuristicVC, MAX_TIME_BUDGET);

            double currentDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endHeuristicWrapper - startGraph).count() /  1000) / (double) 1000;
            //set graph to state that it is in when vc vertices are inactive --> for fastVC() method
            //std::cout << "before setInactive" << std::endl;
            for(auto it = heuristicVC->begin(); it != heuristicVC->end(); ++it)
            {
                bucketGraph->setInactive(it->first);
            }
            //std::cout << "before fastVC" << std::endl;
            // TODO: find suitable timout value for localSearch
            auto localSearchVC = fastVC(bucketGraph, heuristicVC, 5/* MAX_TIME_BUDGET - currentDuration - PRINT_TIME */);
            int localSearchVCSize = localSearchVC->size();
            int heuristicVCSize = heuristicVC->size();
            delete heuristicVC;
            heuristicVC = localSearchVC;

            currentDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endHeuristicWrapper - startGraph).count() /  1000) / (double) 1000;
            auto startPrintSolution = std::chrono::high_resolution_clock::now();

            //cout << "before print" << endl;
            if(!interrupted_by_sig)
            {
                //safely print solution, otherwise wait SIG
                if(currentDuration < MAX_TIME_BUDGET - PRINT_TIME)
                {
                    printing_sol = true;
                    bucketGraph->printVertices(heuristicVC);
                    //cout << "#recursive steps: " << heuristicNumRecursions << endl;
                    cout << "#recursive steps: " << 1000000 + localSearchVCSize - heuristicVCSize/* vcMax->size() - vcMaxRandom->size() */ << endl;
                }
                else
                {
                    while(true) {}
                }
            }

            auto endPrintSolution = std::chrono::high_resolution_clock::now();
            double printDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endPrintSolution - startPrintSolution).count() /  1000) / (double) 1000;
            //std::cout << "Total duration: " << graphConstructionDuration + heuristicWrapperDuration + printDuration << " seconds, Graph construction:" << graphConstructionDuration << " seconds, HeuristicWrapper: " << heuristicWrapperDuration << " seconds, Print solution: " << printDuration << "\n";

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

        std::string solutionSAT = SATSolver.solver();

        if(solutionSAT == "-1")
            return;

        SATSolver.writeOutputSolutionToOutput(solutionSAT);

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
        chooseImplementationAndOutput(0, false, false, false, false, true, false);
        //chooseImplementationAndOutput(1, true, false, false, true, true, false); //print alot
        //signal(SIGINT, my_sig_handler); //catches SIGINT to output anything
        //chooseImplementationAndOutput(5, false, false, false, false, true, false);
        //chooseImplementationAndOutput(6, false, false, false, false, true, false);
    }
	catch (const exception& e)
	{
		cerr << ColorPrint::dye("Error while running vertex cover solver.\n", 'r');
        cerr << ColorPrint::dye(e.what() ,'r');
	}
}