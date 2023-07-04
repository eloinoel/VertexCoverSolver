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


typedef ColorPrint cp;

/*----------------------------------------------------------*/
/*--------------------   Ex 5 Solver   ---------------------*/
/*----------------------------------------------------------*/

int vcVertexBranchingConstrained(BucketGraph* G, std::unordered_map<int, bool>* vc, bool* foundVC, int k, int c, int u, int depth, int* numRec, bool printDebug = false)
{
    (*numRec)++;
    int previousK = k;
    bool cut = false;
    //std::cout << "> cutting through data reductions " << '\n';
    cut = G->dynamicReduce(&k, depth, printDebug);
    if(cut)
    {
        //std::cout << "> cutting through data reductions " << '\n';
        G->unreduce(&k, previousK);
        return u;
    }

    //std::cout << "> calculated LPBound: " << G->getLPBound() << " with k=" << k << '\n';
    if (c + G->getLPBound() >= u) {
        G->unreduce(&k, previousK);
        return u;
    }

    //cout << "before getMaxDegreeVertex" << endl;
	int vertex = G->getMaxDegreeVertex();
    //no vertices left
    if (vertex == -1)
    {
        (*foundVC) = true;
        G->unreduce(&k, previousK, vc);
        return c;
    }
    //cout << "before getVertexDegree: " << vertex << endl;
    int vertexDeg = G->getVertexDegree(vertex);
    //cout << "got maxDegree vertex: " << vertex << endl;
    //G->print();
	//graph has no edges left
	if (vertexDeg == 0)
	{
		(*foundVC) = true;
        G->unreduce(&k, previousK, vc);
        return c;
	}

    // TODO: Solve connected components independently

    //cout << cp::dye("branching: choosing vertex: " + std::to_string(vertex), 'b') << endl;
	//delete first vertex from graph and explore solution
    G->setInactive(vertex);
    bool b1vcFound = false;
    std::unordered_map<int, bool> b1vc = std::unordered_map<int, bool>();
    //cout << "before branching" << endl;
	u = vcVertexBranchingConstrained(G, &b1vc, &b1vcFound, k - 1, c + 1, u, depth+1, numRec);
	if (b1vcFound)
	{
        //revert changes for multiple executions of the algorithm
        vc->insert(std::pair(vertex, true));
        G->setActive(vertex);
        vc->merge(b1vc); //push results
        G->unreduce(&k, previousK, vc); //unreduce needs correct vc für unmerge of deg2rule
		//S->push_back(vertex);
		return u;
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
    if (vertexDeg > u)
    {
        //G->unreduce(&k, previousK);
        return u;
    }

    //cout << "deleting neighbourhood of vertex " << vertex << ": ";
    vector<int>* neighbours = G->getNeighbours(vertex);
    /* cout << ColorPrint::dye("branching: choosing neighbours of vertex " + std::to_string(vertex) + ": ", 'b');
    for(int i = 0; i < (int) neighbours->size(); i++)
    {
        cout <<  ColorPrint::dye(std::to_string(neighbours->at(i)) + ", ", 'b');
    }
    cout << endl; */
    G->setInactive(neighbours);
    bool b2vcFound = false;
    std::unordered_map<int, bool> b2vc = std::unordered_map<int, bool>();
    //cout << "pinc " << '\n';
	u = vcVertexBranchingConstrained(G, &b2vc, &b2vcFound, k - neighbours->size(), c + neighbours->size(), u, depth+1, numRec);
    //cout << "prec " << '\n';
	if (b2vcFound)
	{
        //revert changes for multiple executions of the algorithm
        G->setActive(neighbours);
        for(int i = 0; i <= (int) neighbours->size(); i++)
        {
            vc->insert(std::pair(neighbours->at(i), true));
        }
        vc->merge(b1vc); //push results
        G->unreduce(&k, previousK, vc); //unreduce needs correct vc für unmerge of deg2rule
        return u;
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
    return u;
}

unordered_map<int, bool>* vcSolverConstrained(BucketGraph* G, int* numRec, bool printDebug)
{
    int numPreprocessingVCVertices = 0;
	int k = 0;
    // Apply Reduction Rules for the first time
    auto startPreprocess = std::chrono::high_resolution_clock::now();
    G->preprocess(&numPreprocessingVCVertices, printDebug);
    numPreprocessingVCVertices = -numPreprocessingVCVertices;
    auto endPreprocess = std::chrono::high_resolution_clock::now();
    double preprocessDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endPreprocess - startPreprocess).count() /  1000) / (double) 1000;
    if(printDebug)
        std::cout << "#Preprocessed Graph to size n=" << G->getNumVertices() << ", m=" << G->getNumEdges() << " in " << preprocessDuration << " seconds" << " (reduced by " << numPreprocessingVCVertices << " vertices)" << std::endl;

    // Get upper bound
    int u = getUpperBound(G, 5);

    bool foundVC = false;
	unordered_map<int, bool>* vc = new unordered_map<int, bool>();
    auto startBranching = std::chrono::high_resolution_clock::now();
    vcVertexBranchingConstrained(G, vc, &foundVC, G->getTotalNumVertices(), 0, u, 0, numRec, printDebug);
    if (vc != nullptr)
    {
        // Add Reduced Vertices to Vertex Cover
        G->unreduce(&k, vc->size()+numPreprocessingVCVertices, vc);
        auto endBranching = std::chrono::high_resolution_clock::now();
        double branchingDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endPreprocess - startPreprocess).count() /  1000) / (double) 1000;
        if(printDebug)
            std::cout << "#Finished branching in " << branchingDuration << " seconds" << std::endl;
        return vc;
    }
}

int getUpperBound(BucketGraph* G, double timeCap)
{
    const double MAX_TIME_BUDGET = timeCap;
    const double HEURISTIC_SOLVER_TIME_CAP = (MAX_TIME_BUDGET*2)/3; //in seconds
    const int NUM_RANDOM_SOLUTION_GENERATIONS = 30;

    //std::cout << "before maxHeuristicSolver" << std::endl;
    int heuristicNumRecursions = 0;
    auto startHeuristicWrapper = std::chrono::high_resolution_clock::now();
    //first generate a fast heuristic solution
    //cout << "before heuristic solver" << endl;
    unordered_map<int, bool>* heuristicVC = maxHeuristicSolver(G, &heuristicNumRecursions, false, false);
    //see if we can find a better initial solution 
    heuristicVC = chooseSmallestHeuristicSolution(G, &heuristicNumRecursions, &heuristicVC, true, true, NUM_RANDOM_SOLUTION_GENERATIONS, HEURISTIC_SOLVER_TIME_CAP, false);
    auto endHeuristicWrapper = std::chrono::high_resolution_clock::now();
    double heuristicWrapperDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endHeuristicWrapper - startHeuristicWrapper).count() /  1000) / (double) 1000;
    //auto localSearchVC = fastVC(bucketGraph, heuristicVC, MAX_TIME_BUDGET);


    //set graph to state that it is in when vc vertices are inactive --> for fastVC() method
    //std::cout << "before setInactive" << std::endl;
    for(auto it = heuristicVC->begin(); it != heuristicVC->end(); ++it)
    {
        G->setInactive(it->first);
    }
    double currentDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endHeuristicWrapper - startHeuristicWrapper).count() /  1000) / (double) 1000;
    //std::cout << "before fastVC" << std::endl;
    // TODO: find suitable timout value for localSearch
    const int LOCAL_SEARCH_TIME_CAP = MAX_TIME_BUDGET - currentDuration;
    auto localSearchVC = fastVC(G, heuristicVC, &heuristicNumRecursions, LOCAL_SEARCH_TIME_CAP);
    int localSearchVCSize = localSearchVC->size();
    delete heuristicVC;
    heuristicVC = localSearchVC;
    delete localSearchVC;
    return localSearchVCSize;
}

/*----------------------------------------------------------*/
/*---------------   Exercise 3 Solver Code   ---------------*/
/*----------------------------------------------------------*/

unordered_map<int, bool>* vcVertexBranchingRecursive(BucketGraph* G, int k, int depth, int* numRec, bool printDebug = false)
{
    (*numRec)++;
	if (k < 0)
    {
		return nullptr;
	}
    int previousK = k;
    bool cut = false;
    //std::cout << "> cutting through data reductions " << '\n';
    cut = G->dynamicReduce(&k, depth, printDebug);
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

unordered_map<int, bool>* vcSolverRecursive(BucketGraph* G, int* numRec, bool printDebug)
{
    int numPreprocessingVCVertices = 0;
	int k = 0;
    // Apply Reduction Rules for the first time
    auto startPreprocess = std::chrono::high_resolution_clock::now();
    G->preprocess(&numPreprocessingVCVertices, printDebug);
    numPreprocessingVCVertices = -numPreprocessingVCVertices;
    auto endPreprocess = std::chrono::high_resolution_clock::now();
    double preprocessDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endPreprocess - startPreprocess).count() /  1000) / (double) 1000;
    if(printDebug)
        std::cout << "#Preprocessed Graph to size n=" << G->getNumVertices() << ", m=" << G->getNumEdges() << " in " << preprocessDuration << " seconds" << " (reduced by " << numPreprocessingVCVertices << " vertices)" << std::endl;
    auto startLowerBound = std::chrono::high_resolution_clock::now();
    k = G->getLowerBoundVC();
    auto endLowerBound = std::chrono::high_resolution_clock::now();
    double lowerBoundDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endLowerBound - startLowerBound).count() /  1000) / (double) 1000;
    if(printDebug)
    {
        std::cout << "#Calculated lower bound k=" << k << " in " << lowerBoundDuration << " seconds" << std::endl;
    }
	unordered_map<int, bool>* vc;
    auto startBranching = std::chrono::high_resolution_clock::now();
	while (true)
	{

        // Begin Branching
		vc = vcVertexBranchingRecursive(G, k, 0, numRec, printDebug);
		if (vc != nullptr)
		{
            // Add Reduced Vertices to Vertex Cover
            G->unreduce(&k, vc->size()+numPreprocessingVCVertices, vc);
            auto endBranching = std::chrono::high_resolution_clock::now();
            double branchingDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endPreprocess - startPreprocess).count() /  1000) / (double) 1000;
            if(printDebug)
                std::cout << "#Finished branching in " << branchingDuration << " seconds" << std::endl;
			return vc;
		}
        k++;
	}
}


/*----------------------------------------------------------*/
/*-------------------   HEURISTIC CODE   -------------------*/
/*----------------------------------------------------------*/

bool outOfTime(std::chrono::time_point<std::chrono::high_resolution_clock> startTime, double timeoutCap)
{
    auto currentTime = std::chrono::high_resolution_clock::now();
    double currentDuration = std::chrono::duration_cast<std::chrono::microseconds>(currentTime - startTime).count() / (double) 1000000;
    return currentDuration > timeoutCap;
}

void resetGraphAfterBranching(BucketGraph* G, std::unordered_map<int, bool>* vc)
{
    for(auto it = vc->begin(); it != vc->end(); ++it)
    {
        G->setActive(it->first);
    }
}

std::unordered_map<int, bool>* maxHeuristicSolver(BucketGraph* G, int* numRec, bool applyReductions = false, bool randomise = false)
{
    // Apply Reduction Rules for the first time
    int numPreprocessingVCVertices = 0;
    if(applyReductions)
    {
        std::vector<bool> rulesToApply = std::vector<bool> ({true, true, true, false});
        G->preprocess(&numPreprocessingVCVertices, rulesToApply);
        numPreprocessingVCVertices = -numPreprocessingVCVertices;
    }
    std::unordered_map<int, bool>* vc = new std::unordered_map<int, bool>();
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

/*
* will try generating multiple solutions to take the best one
* will stop if current time elapsed given @timeoutSoftCap
* call with nullptr as @currentSmallestVC if want to generate simple and fast max heuristic solution at beginning
//TODO: possible optimisation to only calculate preprocessing once at beginning and copy graph and then just add reduction vertices to solutions
*/
std::unordered_map<int, bool>* chooseSmallestHeuristicSolution(BucketGraph* G, int* numRec, std::unordered_map<int, bool>** currentSmallestVC, bool applyReductions, bool includeRandomsWithReductions, int numRandomSolverCalls,
 int timeoutSoftCap, bool printBestSolution)
{
    int best_solution = 0; //0: maxHeuristic, 1: with preprocessing, 2: randomised, 3: randomised with preprocessing
    auto startTime = std::chrono::high_resolution_clock::now();

    if(*currentSmallestVC == nullptr)
    {
        *currentSmallestVC = maxHeuristicSolver(G, numRec);
    }

    //cout << "maxHeuristicSolver size: " << (*currentSmallestVC)->size() << ", addr: " << *currentSmallestVC << endl;

    if(outOfTime(startTime, timeoutSoftCap)) {
        if(printBestSolution)
        {
            std::cout << "#recursive steps: " << best_solution << std::endl;
        }
        return *currentSmallestVC;
    }

    //preprocessing on
    int vcMaxPreprocessingNumRec = 0;
    std::unordered_map<int, bool>* vcMaxPreprocessing = nullptr;
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

    if(outOfTime(startTime, timeoutSoftCap)) {
        if(printBestSolution)
        {
            std::cout << "#recursive steps: " << best_solution << std::endl;
        }
        return *currentSmallestVC;
    }

    //selecting random max degree vertices
    int vcMaxRandomNumRec;
    std::unordered_map<int, bool>* vcMaxRandom = nullptr;
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

            if(outOfTime(startTime, timeoutSoftCap)) {
                if(printBestSolution)
                {
                    std::cout << "#recursive steps: " << best_solution << std::endl;
                }
                return *currentSmallestVC;
            }
        }
    }
    //cout << "after random size: " << (*currentSmallestVC)->size() << ", addr: " << *currentSmallestVC << endl;
    if(printBestSolution)
    {
        std::cout << "#recursive steps: " << best_solution << std::endl;
    }
    return *currentSmallestVC;
}

void initGainLossAge(BucketGraph* G, std::unordered_map<int, bool>* vc, std::vector<int>* gain, std::vector<int>* loss, std::vector<int>* age)
{
    for (int i = 0; i < (int) gain->size(); ++i)
    {
        (*gain)[i] = 0;
        (*loss)[i] = 0;
        (*age)[i] = 0;
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

int getMinLossIndex(std::unordered_map<int, bool>* vc, std::vector<int>* loss, std::vector<int>* age)
{
    if (vc == nullptr) { throw std::invalid_argument("getMinLossIndex: passed vc is nullptr"); }
    if (vc->empty()) { throw std::invalid_argument("getMinLossIndex: passed vc is empty"); }

    int min = INT32_MAX;
    int minIndex = -1;

    for (auto it = vc->begin(); it != vc->end(); ++it)
    {
        // TODO: use age as tiebreaker?
        if (loss->at(it->first) < min || (loss->at(it->first) == min && age->at(it->first) < age->at(minIndex)))
        {
            min = loss->at(it->first);
            minIndex = it->first;
            continue;
        }
    }
    return minIndex;
}

int getBMSMinLossIndex(std::unordered_map<int, bool>* vc, std::vector<int>* loss, std::vector<int>* age, int k)
{
    if (vc == nullptr) { throw std::invalid_argument("getMinLossIndex: passed vc is nullptr\n"); }
    if (vc->empty()) { throw std::invalid_argument("getMinLossIndex: passed vc is empty\n"); }
    if (k > (int) vc->size()) { throw std::invalid_argument("getMinLossIndex: passed k is larger than passed vc\n"); }
    if (k < 1) { throw std::invalid_argument("getMinLossIndex: passed k that smaller than 1\n"); }

    int min = INT32_MAX;
    int minIndex = -1;

    std::random_device rd;
    std::mt19937 gen(rd());
    int partitionSize = vc->size() / k;
    std::uniform_int_distribution<> distr(0, partitionSize-1);

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

        if (loss->at(it->first) < min || (loss->at(it->first) == min && age->at(it->first) < age->at(minIndex)))
        {
            min = loss->at(it->first);
            minIndex = it->first;
            continue;
        }
    }
    return minIndex;
}

void removeMinLossVCVertex(BucketGraph* G, std::unordered_map<int, bool>* vc, std::vector<int>* gain, std::vector<int>* loss, std::vector<int>* age, int numRecursions)
{
    int u_index = getMinLossIndex(vc, loss, age);
    auto u = vc->find(u_index);
    /* std::cout << "uIndex: " << u_index << std::endl;
    std::cout << "vc size: " << vc->size() << std::endl;
    int i=0;
    for (auto it = vc->begin(); it != vc->end(); ++it)
    {
        std::cout << "vc[" << i << "]: " << it->first << std::endl;
        i++;
    } */
    if (u == vc->end()) { throw std::invalid_argument("removeMinLossVCVertex: Iterator of u not found\n"); }
    vc->erase(u);
    G->setActive(u_index);
    (*gain)[u_index] = (*loss)[u_index];
    (*loss)[u_index] = 0;
    (*age)[u_index] = numRecursions;
    Vertex* uVertex = G->getVertex(u_index);
    for (auto neighbour = uVertex->getAdj()->begin(); neighbour != uVertex->getAdj()->end(); ++neighbour)
    {
        if (!(G->getVertex(*neighbour)->getActive())) { continue; }
        (*gain)[*neighbour]++;
    }
    //std::cout << cp::dye("Removed " + std::to_string(u_index), 'r') << " from vc" << std::endl;
}

void removeBMSMinLossVCVertex(BucketGraph* G, std::unordered_map<int, bool>* vc, std::vector<int>* gain, std::vector<int>* loss, std::vector<int>* age, int numRecursions, int k)
{
    int u_index = getBMSMinLossIndex(vc, loss, age, k);
    auto u = vc->find(u_index);

    if (u == vc->end()) { throw std::invalid_argument("removeMinLossVCVertex: Iterator of u not found\n"); }
    vc->erase(u);
    G->setActive(u_index);
    (*gain)[u_index] = (*loss)[u_index];
    (*loss)[u_index] = 0;
    (*age)[u_index] = numRecursions;
    Vertex* uVertex = G->getVertex(u_index);
    for (auto neighbour = uVertex->getAdj()->begin(); neighbour != uVertex->getAdj()->end(); ++neighbour)
    {
        if (!(G->getVertex(*neighbour)->getActive())) { continue; }
        (*gain)[*neighbour]++;
    }
}

void addRandomUncoveredEdgeMaxGainEndpointVertex(BucketGraph* G, std::unordered_map<int, bool>* vc, std::vector<int>* gain, std::vector<int>* loss, std::vector<int>* age, int numRecursions)
{
    //std::cout << "Adding to vc..." << std::endl;
    //int v = G->getRandomConnectedVertex(10);    // TODO: consider pseudo randomness (are all vertices selectable?)
    int v = G->getRandomConnectedVertex(INT32_MAX); // TODO: test more what limit is appropriate

    if (v == -1) { return; }
    Vertex* vertex = G->getVertex(v);
    if (vertex->getDegree() == 0) { throw std::invalid_argument("addRandomUncoveredEdgeMaxGainEndpointVertex: Got vertex of degree 0\n"); }
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
    if((*gain)[v] > (*gain)[neighbour] || ((*gain)[v] == (*gain)[neighbour] && (*age)[v] <= (*age)[neighbour]))
    {
        for (auto u = vertex->getAdj()->begin(); u != vertex->getAdj()->end(); ++u)
        {
            if (!(G->getVertex(*u)->getActive())) { continue; }
            (*gain)[*u]--;
            (*loss)[v]++;
        }
        G->setInactive(v);
        vc->insert({v, true});
        (*gain)[v] = 0;
        (*age)[v] = numRecursions;
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
        (*age)[neighbour] = numRecursions;
        //std::cout << cp::dye("Added " + std::to_string(neighbour), 'g') << " to vc" << std::endl;
    }
}

/*
*   vertices of vertex cover vc need to be inactive in graph G
*   timeout in seconds
*   vc is not changed, you need to update it manually after the function call
*/
std::unordered_map<int, bool>* fastVC(BucketGraph* G, std::unordered_map<int, bool>* vc, int* numRec, double timeout)
{
    //std::cout << "before guards" << std::endl;
    if (vc == nullptr) { throw std::invalid_argument("fastVC: passed vc is nullptr\n"); }
    for (auto it = vc->begin(); it != vc->end(); ++it)
    {
        if (G->getVertex(it->first)->getActive())
        {
            throw std::invalid_argument("fastVC: passed vc vertices must be set to inactive in graph\n");
        }
    }

    //std::cout << "before inits" << std::endl;
    int numRecursions = 0;
    int n = G->getTotalNumVertices();
    std::vector<int> gain = std::vector<int>(n);
    std::vector<int> loss = std::vector<int>(n);
    std::vector<int> age = std::vector<int>(n);
    std::unordered_map<int, bool>* currentVC = new std::unordered_map<int, bool>(*vc);
    std::unordered_map<int, bool>* bestVC = new std::unordered_map<int, bool>(*vc);
    initGainLossAge(G, currentVC, &gain, &loss, &age);
    auto startFastVC = std::chrono::high_resolution_clock::now();
    while (true && vc->size() > 1)
    {
        if (outOfTime(startFastVC, timeout)) { break; }
        auto foundVCCheckStart = std::chrono::high_resolution_clock::now();
        numRecursions++;
        // we found an improved vc! time to overwrite our previous best solution and search further
        if (G->getMaxDegree() <= 0)
        {
            //std::cout << cp::dye("Found vc of size: " + std::to_string(currentVC->size()), 'g') << std::endl;
            delete bestVC;
            bestVC = new std::unordered_map<int, bool>(*currentVC);
            removeMinLossVCVertex(G, currentVC, &gain, &loss, &age, numRecursions);
        }
        // TODO: test different k's and find suitable k value


        auto rmMinLossStart = std::chrono::high_resolution_clock::now();
        //removeBMSMinLossVCVertex(G, currentVC, &gain, &loss, &age, numRecursions, std::min(500000, (int) currentVC->size())/* std::max(1, (int) currentVC->size() / 4) */);
        removeMinLossVCVertex(G, currentVC, &gain, &loss, &age, numRecursions);
        auto addUncEdgeEndpointStart = std::chrono::high_resolution_clock::now();
        addRandomUncoveredEdgeMaxGainEndpointVertex(G, currentVC, &gain, &loss, &age, numRecursions);
        auto addUncEdgeEndpointEnd = std::chrono::high_resolution_clock::now();
        /* double foundVCCheck = std::chrono::duration_cast<std::chrono::microseconds>(rmMinLossStart - foundVCCheckStart).count() / (double) 1000000;
        double rmMinLoss = std::chrono::duration_cast<std::chrono::microseconds>(addUncEdgeEndpointStart - rmMinLossStart).count() / (double) 1000000;
        double addUncEdgeEndpoint = std::chrono::duration_cast<std::chrono::microseconds>(addUncEdgeEndpointEnd - addUncEdgeEndpointStart).count() / (double) 100000;
        std::cout << "Completed local search iteration in " << foundVCCheck + rmMinLoss + addUncEdgeEndpoint << " seconds (VCCheck: " << foundVCCheck << " + rmMinLoss: " << rmMinLoss << " + addUncEdgeEndpoint: " << addUncEdgeEndpoint << ")" << '\n'; */
    }
    (*numRec) = (*numRec) + numRecursions;
    delete currentVC;
    return bestVC;
}

