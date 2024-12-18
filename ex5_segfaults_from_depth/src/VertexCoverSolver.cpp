#include <iostream>
#include <fstream>  //ifstream file opening
#include <stack>          // std::stack
#include <math.h>          // INFINITY
#include "utils/ColorPrint.h"
#include "utils/BucketGraph.h"
#include "utils/Packing.h"
#include <chrono>
#include <random>


#include <signal.h>
#include <stdlib.h>
#include <stdio.h>


typedef ColorPrint cp;

/*----------------------------------------------------------*/
/*---------------   Exercise 3 Solver Code   ---------------*/
/*----------------------------------------------------------*/

std::unordered_map<int, bool>* vcVertexBranchingRecursive(BucketGraph* G, int k, int depth, int* numRec, bool printDebug = false)
{

    (*numRec)++;
    // Save recursion depth to unreduce degree-3 at correct depth
	int currentRecursion = *numRec;

    if (k < 0)
    {
		return nullptr;
	}
    int previousK = k;
    bool cut = false;
    if(printDebug)
        std::cout << "#--> applying data reductions " << '\n';
    cut = G->dynamicReduce(&k, depth, printDebug);
    if(cut)
    {
        if(printDebug)
            std::cout << "#--> cutting through data reductions " << '\n';
        G->unreduce(&k, previousK, depth);
        return nullptr;
    }

    int lowerBound = G->getLPBound();
    if(printDebug)
        std::cout << "#--> calculated LPBound: " << lowerBound << " with k=" << k << '\n';
    if (k < lowerBound) {
        G->unreduce(&k, previousK, depth);
        return nullptr;
    }

    //cout << "before getMaxDegreeVertex" << endl;
	int vertex = G->getMaxDegreeVertex();
    //no vertices left
    if (vertex == -1)
    {
        std::unordered_map<int, bool>* vc = new std::unordered_map<int, bool>();
        G->unreduce(&k, previousK, depth, vc);
        return vc;
    }
    //cout << "before getVertexDegree: " << vertex << endl;
    int vertexDeg = G->getVertexDegree(vertex);
    //cout << "got maxDegree vertex: " << vertex << endl;
    //G->print();
	//graph has no edges left
	if (vertexDeg == 0)
	{
		std::unordered_map<int, bool>* vc = new std::unordered_map<int, bool>();
        G->unreduce(&k, previousK, depth, vc);
        return vc;
	}

    //cout << cp::dye("branching: choosing vertex: " + std::to_string(vertex), 'b') << endl;
	//delete first vertex from graph and explore solution
    G->setInactive(vertex);

    //int pDeg3 = G->period_deg3;
    //int pUnc = G->period_unc;
    //int pLP = G->period_lp;

    //cout << "before branching" << endl;
	std::unordered_map<int, bool>* S = vcVertexBranchingRecursive(G, k - 1, depth+1, numRec);
	if (S != nullptr)
	{
        //revert changes for multiple executions of the algorithm
        G->setActive(vertex);
        S->insert({vertex, true}); //push results
        G->unreduce(&k, previousK, depth, S); //unreduce needs correct vc für unmerge of deg2rule
		//S->push_back(vertex);
		return S;
	}
	else
	{
        //cout << "before setActive" << endl;
		//revert changes to graph
		G->setActive(vertex);
	}
    //cout << cp::dye("restoring vertex: ", 'g') << vertex << endl;

    //G->period_deg3 = pDeg3;
    //G->period_unc = pUnc;
    //G->period_lp = pLP;

	//cannot fully explore neighbours
    if (vertexDeg > k)
    {
        G->unreduce(&k, previousK, depth);
        return nullptr;
    }

    //cout << "deleting neighbourhood of vertex " << vertex << ": ";
    std::vector<int>* neighbours = G->getNeighbours(vertex);
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
            S->insert({neighbours->at(i), true}); //push results
        }
        G->unreduce(&k, previousK, depth, S); //unreduce needs correct vc für unmerge of deg2rule
        return S;
	}
	else
	{
		//revert changes to graph
		G->setActive(neighbours);
	}
    /* cout << "restoring neighbourhood of vertex " << vertex << ": ";
    for(int i = 0; i < (int) neighbours->size(); i++)
    {
        cout << neighbours->at(i) << ", ";
    }
    cout << endl; */
    // free neighbours
    delete neighbours;

    //G->initRulePeriods();

    G->unreduce(&k, previousK, depth);
    return nullptr;
}

std::unordered_map<int, bool>* vcSolverRecursive(BucketGraph* G, int* numRec, bool printDebug)
{
    int numPreprocessingVCVertices = 0;
	int k = 0;

    // Apply Reduction Rules for the first time
    auto startPreprocess = std::chrono::high_resolution_clock::now();
                                                    //    0    1     2     3     4    5       6       7     8     9     10   11
    std::vector<bool> rulesToApply = std::vector<bool>{true, true, false, true, true, false, false, false, false, false, false, false};
    G->preprocess(&numPreprocessingVCVertices, rulesToApply, printDebug);
    numPreprocessingVCVertices = -numPreprocessingVCVertices;
    auto endPreprocess = std::chrono::high_resolution_clock::now();
    double preprocessDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endPreprocess - startPreprocess).count() /  1000) / (double) 1000;
    if(printDebug)
        std::cout << "#Preprocessed Graph to size n=" << G->getNumVertices() << ", m=" << G->getNumEdges() << " in " << preprocessDuration << " seconds" << " (reduced by " << numPreprocessingVCVertices << " vertices)" << std::endl;
    
    auto startLowerBound = std::chrono::high_resolution_clock::now();
    k = G->getLowerBoundVC();
    //G->upperBound = 2 * k;
    //G->initRulePeriods();
//    std::cout << "# k =" << k << '\n';


    auto endLowerBound = std::chrono::high_resolution_clock::now();
    double lowerBoundDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endLowerBound - startLowerBound).count() /  1000) / (double) 1000;
    if(printDebug)
    {
        std::cout << "#Calculated lower bound k=" << k << " in " << lowerBoundDuration << " seconds" << std::endl;
    }

	std::unordered_map<int, bool>* vc;
    auto startBranching = std::chrono::high_resolution_clock::now();
	while (true)
	{
//        std::cout << "# Trying k =" << k << '\n';
        // Begin Branching
		vc = vcVertexBranchingRecursive(G, k, 0, numRec, printDebug);
		if (vc != nullptr)
		{
            // Add Reduced Vertices to Vertex Cover
            G->unreduce(&k, vc->size()+numPreprocessingVCVertices, -1, vc);
            if(G->getReductionStackSize() > 0)
            {
                throw std::invalid_argument("vcSolverRecursive: reduction rule stack isn't fully popped after final unreduce");
            }
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
                                                            //    0    1     2     3     4    5       6       7     8     9     10   11
        //std::vector<bool> rulesToApply = std::vector<bool>{true, true, false, true, true, false, false, true, true, true, true, true};
        std::vector<bool> rulesToApply = std::vector<bool> ({true, true, true, false, false, false , false , false , false, false, true, true});
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
        G->unreduce(&k, vc->size()+numPreprocessingVCVertices, -1, vc);
        if(G->getReductionStackSize() > 0)
        {
            throw std::invalid_argument("vcSolverRecursive: reduction rule stack isn't fully popped after final unreduce");
        }
        G->resetMatching();
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
            int vertex = u->first;
            if(vc->find(vertex) == vc->end())
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
    if (u == vc->end()) { throw std::invalid_argument("removeMinLossVCVertex: Iterator of u not found\n"); }
    vc->erase(u);
    G->setActive(u_index);
    (*gain)[u_index] = (*loss)[u_index];
    (*loss)[u_index] = 0;
    (*age)[u_index] = numRecursions;
    Vertex* uVertex = G->getVertex(u_index);
    for (auto neighbourIt = uVertex->getAdj()->begin(); neighbourIt != uVertex->getAdj()->end(); ++neighbourIt)
    {
        int neighbour = neighbourIt->first;
        if (!(G->getVertex(neighbour)->getActive())) { continue; }
        (*gain)[neighbour]++;
    }
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
    for (auto neighbourIt = uVertex->getAdj()->begin(); neighbourIt != uVertex->getAdj()->end(); ++neighbourIt)
    {
        int neighbour = neighbourIt->first;
        if (!(G->getVertex(neighbour)->getActive())) { continue; }
        (*gain)[neighbour]++;
    }
}

void addRandomUncoveredEdgeMaxGainEndpointVertex(BucketGraph* G, std::unordered_map<int, bool>* vc, std::vector<int>* gain, std::vector<int>* loss, std::vector<int>* age, int numRecursions)
{
    //int v = G->getRandomConnectedVertex(10);    // TODO: consider pseudo randomness (are all vertices selectable?)
    int v = G->getRandomConnectedVertex(INT32_MAX); // TODO: test more what limit is appropriate

    if (v == -1) { return; }
    Vertex* vertex = G->getVertex(v);
    if (vertex->getDegree() == 0) { throw std::invalid_argument("addRandomUncoveredEdgeMaxGainEndpointVertex: Got vertex of degree 0\n"); }

    auto neighbours = G->getNeighbours(v);  // TODO: slow getNeighbours
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0, neighbours->size()-1);
    int rnd = (int) distr(gen);
    int neighbour = neighbours->at(rnd);
    Vertex* neighbourVertex = G->getVertex(neighbour);

    // add vertex with higher gain into vc
    if((*gain)[v] > (*gain)[neighbour] || ((*gain)[v] == (*gain)[neighbour] && (*age)[v] <= (*age)[neighbour]))
    {
        for (auto u = vertex->getAdj()->begin(); u != vertex->getAdj()->end(); ++u)
        {
            int vert = u->first;
            if (!(G->getVertex(vert)->getActive())) { continue; }
            (*gain)[vert]--;
            (*loss)[v]++;
        }
        G->setInactive(v);
        vc->insert({v, true});
        (*gain)[v] = 0;
        (*age)[v] = numRecursions;
    }
    else
    {
        for (auto u = neighbourVertex->getAdj()->begin(); u != neighbourVertex->getAdj()->end(); ++u)
        {
            int vert = u->first;
            if (!(G->getVertex(vert)->getActive())) { continue; }
            (*gain)[vert]--;
        }
        G->setInactive(neighbour);
        vc->insert({neighbour, true});
        (*gain)[neighbour] = 0;
        (*age)[neighbour] = numRecursions;
    }
}

/*
*   vertices of vertex cover vc need to be inactive in graph G
*   timeout in seconds
*   vc is not changed, you need to update it manually after the function call
*/
std::unordered_map<int, bool>* fastVC(BucketGraph* G, std::unordered_map<int, bool>* vc, int* numRec, double timeout)
{
    if (vc == nullptr) { throw std::invalid_argument("fastVC: passed vc is nullptr\n"); }
    for (auto it = vc->begin(); it != vc->end(); ++it)
    {
        if (G->getVertex(it->first)->getActive())
        {
            throw std::invalid_argument("fastVC: passed vc vertices must be set to inactive in graph\n");
        }
    }

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
    // TODO: resetting graph may be inefficient or redundant
    for(auto it = currentVC->begin(); it != currentVC->end(); ++it)
    {
        if(G->isActive(it->first)) { continue; }
        G->setActive(it->first);
    }
    for(auto it = bestVC->begin(); it != bestVC->end(); ++it)
    {
        if(!G->isActive(it->first)) { continue; }
        G->setInactive(it->first);
    }
    delete currentVC;
    return bestVC;
}

/*----------------------------------------------------------*/
/*--------------------   Ex 5 Solver   ---------------------*/
/*----------------------------------------------------------*/

/*
    G should be preprocessed
*/
int getUpperBoundFromRandomHeuristic(BucketGraph* graphCopy, BucketGraph* G, int numPreprocessingVCVertices, double timeoutSoftCap)
{
    int vcNumRec = 0;
    std::unordered_map<int, bool>* currentBestVC = nullptr;
    int currentBestVCSize = INT_MAX;

    //start generating solutions
    auto startTime = std::chrono::high_resolution_clock::now();
    currentBestVC = maxHeuristicSolver(G, &vcNumRec);
    currentBestVCSize = (int)currentBestVC->size() + numPreprocessingVCVertices;

    if(outOfTime(startTime, timeoutSoftCap)) {
        delete currentBestVC;
        return currentBestVCSize;
    }

    //selecting random max degree vertices
    std::unordered_map<int, bool>* vcRandom = nullptr;
    for(int i = 0; i < INT_MAX; ++i)
    {
        vcNumRec = 0;
        int vcSize = INT_MAX;
        if(i % 3 == 1) //every third random iteration doesn't use preprocessed graph
        {
            vcRandom = maxHeuristicSolver(graphCopy, &vcNumRec, false, true);
            vcSize = vcRandom->size();
        }
        else
        {
            //graph with preprocessing is used
            vcRandom = maxHeuristicSolver(G, &vcNumRec, false, true);
            vcSize = vcRandom->size() + numPreprocessingVCVertices;
        }

        if(vcSize < currentBestVCSize)
        {
            delete currentBestVC;
            currentBestVC = vcRandom;
            currentBestVCSize = vcSize;
        }
        else
        {
            delete vcRandom;
        }

        if(outOfTime(startTime, timeoutSoftCap)) {
            delete currentBestVC;
            return currentBestVCSize;
        }
    }

    delete currentBestVC;
    return currentBestVCSize;
}

/*
    G should be preprocessed
*/
std::unordered_map<int, bool>* getUpperBoundFromRandomHeuristicVC(BucketGraph* graphCopy, BucketGraph* G, int numPreprocessingVCVertices, double timeoutSoftCap)
{
    int vcNumRec = 0;
    std::unordered_map<int, bool>* currentBestVC = nullptr;
    int currentBestVCSize = INT_MAX;
    bool currentBestVCIsPreprocessed = true;
    bool randomVCIsPreprocessed = false;

    // preprocessing
    int k = 0;
    BucketGraph* GC = graphCopy->copy();
    GC->preprocess(&k);

    //start generating solutions
    auto startTime = std::chrono::high_resolution_clock::now();
    currentBestVC = maxHeuristicSolver(G, &vcNumRec);
    currentBestVCSize = (int)currentBestVC->size() + numPreprocessingVCVertices;

    if(outOfTime(startTime, timeoutSoftCap)) {
        if(currentBestVCIsPreprocessed)
        {
            GC->unreduce(&k, 0, -1, currentBestVC);
        }
        return currentBestVC;
    }

    //selecting random max degree vertices
    std::unordered_map<int, bool>* vcRandom = nullptr;
    for(int i = 0; i < INT_MAX; ++i)
    {
        vcNumRec = 0;
        int vcSize = INT_MAX;
        if(i % 3 == 1) //every third random iteration doesn't use preprocessed graph
        {
            vcRandom = maxHeuristicSolver(graphCopy, &vcNumRec, false, true);
            vcSize = vcRandom->size();
            randomVCIsPreprocessed = false;
        }
        else
        {
            //graph with preprocessing is used
            vcRandom = maxHeuristicSolver(GC, &vcNumRec, false, true);
            vcSize = vcRandom->size() + numPreprocessingVCVertices;
            randomVCIsPreprocessed = true;
            GC->preprocess(&k);
        }

        if(vcSize < currentBestVCSize)
        {
            delete currentBestVC;
            currentBestVC = vcRandom;
            currentBestVCSize = vcSize;
            currentBestVCIsPreprocessed = randomVCIsPreprocessed;
        }
        else
        {
            delete vcRandom;
        }

        if(outOfTime(startTime, timeoutSoftCap)) {
            if(currentBestVCIsPreprocessed)
            {
                GC->unreduce(&k, 0, -1, currentBestVC);
            }
            return currentBestVC;
        }
    }

    if(currentBestVCIsPreprocessed)
    {
        GC->unreduce(&k, 0, -1, currentBestVC);
    }
    return currentBestVC;
}

/* G should already be preprocessed */
int getUpperBound(BucketGraph* graphCopy, BucketGraph* G, int numPreprocessingVCVertices, double timeCap, bool printDebug)
{
    const double MAX_TIME_BUDGET = timeCap;
    const double HEUR_INITIAL_SOLUTION_TIME_CAP = (MAX_TIME_BUDGET*2)/3; //in seconds

    auto startHeuristicWrapper = std::chrono::high_resolution_clock::now();

    //first generate a fast heuristic solution
    int heuristicNumRecursions = 0;
    std::unordered_map<int, bool>* heuristicVC = maxHeuristicSolver(graphCopy, &heuristicNumRecursions, false, false);
    int currentBestVCSize = heuristicVC->size();
    
    //see if we can find a better heuristical solution
    int randomHeurSize = getUpperBoundFromRandomHeuristic(graphCopy, G, numPreprocessingVCVertices, HEUR_INITIAL_SOLUTION_TIME_CAP);
    if(randomHeurSize < currentBestVCSize)
        currentBestVCSize = randomHeurSize;
    auto endHeuristicWrapper = std::chrono::high_resolution_clock::now();
    double heuristicWrapperDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endHeuristicWrapper - startHeuristicWrapper).count() /  1000) / (double) 1000;

    if(printDebug) { std::cout << "#getUpperBound: initial best heuristical solution size: " << currentBestVCSize << std::endl; }

    //set graph to state that it is in when vc vertices are inactive --> for fastVC() method
    for(auto it = heuristicVC->begin(); it != heuristicVC->end(); ++it)
    {
        graphCopy->setInactive(it->first);
    }
    const int LOCAL_SEARCH_TIME_CAP = MAX_TIME_BUDGET - heuristicWrapperDuration;

    auto localSearchVC = fastVC(graphCopy, heuristicVC, &heuristicNumRecursions, LOCAL_SEARCH_TIME_CAP);
    int localSearchVCSize = localSearchVC->size();
    if(localSearchVCSize < currentBestVCSize)
        currentBestVCSize = localSearchVCSize;

    /* //no need because graphcopy
     for(auto it = localSearchVC->begin(); it != localSearchVC->end(); ++it)
    {
        if(G->isActive(it->first)) { continue; }
        G->setActive(it->first);
    } */

    if(printDebug) { std::cout << "#getUpperBound: fast vc solution size: " << localSearchVCSize << std::endl; }

    delete heuristicVC;
    delete localSearchVC;
    return currentBestVCSize;
}

/* G should already be preprocessed */
std::unordered_map<int, bool>* getUpperBoundVC(BucketGraph* graphCopy, BucketGraph* G, int numPreprocessingVCVertices, double timeCap, bool printDebug)
{
    const double MAX_TIME_BUDGET = timeCap;
    const double HEUR_INITIAL_SOLUTION_TIME_CAP = (MAX_TIME_BUDGET*2)/3; //in seconds

    auto startHeuristicWrapper = std::chrono::high_resolution_clock::now();

    //first generate a fast heuristic solution
    int heuristicNumRecursions = 0;
    std::unordered_map<int, bool>* heuristicVC = maxHeuristicSolver(graphCopy, &heuristicNumRecursions, false, false);
    int heuristicVCSize = heuristicVC->size();
    int currentBestVCSize = heuristicVC->size();
    std::unordered_map<int, bool>* currentBestVC = heuristicVC;

    if(printDebug) { std::cout << "#getUpperBound: initial heuristical solution size: " << heuristicVCSize << std::endl; }

    //see if we can find a better heuristical solution
    std::unordered_map<int, bool>* randomHeurVC = getUpperBoundFromRandomHeuristicVC(graphCopy, G, numPreprocessingVCVertices, HEUR_INITIAL_SOLUTION_TIME_CAP);
    int randomHeurVCSize = randomHeurVC->size();
    if(randomHeurVCSize < currentBestVCSize) {
        delete currentBestVC;
        currentBestVCSize = randomHeurVCSize;
        currentBestVC = randomHeurVC;
    }
    auto endHeuristicWrapper = std::chrono::high_resolution_clock::now();
    double heuristicWrapperDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endHeuristicWrapper - startHeuristicWrapper).count() /  1000) / (double) 1000;

    if(printDebug) { std::cout << "#getUpperBound: initial best heuristical solution size: " << randomHeurVCSize << std::endl; }

    //set graph to state that it is in when vc vertices are inactive --> for fastVC() method
    for(auto it = currentBestVC->begin(); it != currentBestVC->end(); ++it)
    {
        graphCopy->setInactive(it->first);
    }
    const int LOCAL_SEARCH_TIME_CAP = MAX_TIME_BUDGET - heuristicWrapperDuration;

    auto localSearchVC = fastVC(graphCopy, currentBestVC, &heuristicNumRecursions, LOCAL_SEARCH_TIME_CAP);
    int localSearchVCSize = localSearchVC->size();
    if(localSearchVCSize < currentBestVCSize) {
        delete currentBestVC;
        currentBestVCSize = localSearchVCSize;
        currentBestVC = localSearchVC;
    }

    /* //no need because graphcopy
     for(auto it = localSearchVC->begin(); it != localSearchVC->end(); ++it)
    {
        if(G->isActive(it->first)) { continue; }
        G->setActive(it->first);
    } */

    if(printDebug) { std::cout << "#getUpperBound: fast vc solution size: " << localSearchVCSize << std::endl; }

    return currentBestVC;
}

std::pair<int, std::unordered_map<int, bool>*> vcVertexBranchingConstrained(BucketGraph* G, int k, int c, int u, int depth, int* numRec,
    Packing* constraints = nullptr, bool printDebug = false)
{
    (*numRec)++;
    int previousK = k;
    bool cut = false;
    cut = G->dynamicReduce(&k, depth, printDebug);
    if(cut)
    {
        G->unreduce(&k, previousK, depth);
        return std::pair<int, std::unordered_map<int, bool>*>(u, nullptr);
    }
    c += previousK - k;

    if (c + G->getLPBound() >= u) {
        G->unreduce(&k, previousK, depth);
        return std::pair<int, std::unordered_map<int, bool>*>(u, nullptr);
    }

	int vertex = G->getMaxDegreeVertex();
    //no vertices left
    if (vertex == -1)
    {
        std::unordered_map<int, bool>* vc = new std::unordered_map<int, bool>();
        G->unreduce(&k, previousK, depth, vc);
        return std::pair<int, std::unordered_map<int, bool>*>(c, vc);
    }
    int vertexDeg = G->getVertexDegree(vertex);
    //G->print();
	//graph has no edges left
	if (vertexDeg == 0)
	{
        std::unordered_map<int, bool>* vc = new std::unordered_map<int, bool>();
        G->unreduce(&k, previousK, depth, vc);
        return std::pair<int, std::unordered_map<int, bool>*>(c, vc);
	}

    std::pair<int, std::unordered_map<int, bool>*> firstSolution = std::pair<int, std::unordered_map<int, bool>*>(u, nullptr);
    std::pair<int, std::unordered_map<int, bool>*> secondSolution = std::pair<int, std::unordered_map<int, bool>*>(u, nullptr);;

    std::pair<int, std::unordered_map<int, bool>*> results;
    std::unordered_map<int, bool>* vc = nullptr;
    if(constraints->checkNeighbourhoodConstraints(vertex, G))
    {
        //delete first vertex from graph and explore solution
        G->setInactive(vertex);
        results = vcVertexBranchingConstrained(G, k - 1, c + 1, u, depth+1, numRec);
        if(u > results.first) {
            k = k - u + results.first;
            previousK = previousK - u + results.first;
            u = results.first;
        }
        vc = results.second;
        if (vc != nullptr)
        {
            //revert changes for multiple executions of the algorithm
            vc->insert(std::pair(vertex, true));
            G->setActive(vertex);
            //G->unreduce(&k, previousK, vc); //unreduce needs correct vc für unmerge of deg2rule
            //G->printVC(vc);
            firstSolution = std::pair<int, std::unordered_map<int, bool>*>(u, vc);
        }
        else
        {
            //revert changes to graph
            G->setActive(vertex);
            //constraints->undoNeighbourhoodConstraints(vertex, G);
            //G->unreduce(&k, previousK);
            firstSolution = std::pair<int, std::unordered_map<int, bool>*>(u, nullptr);
        }
        //cout << cp::dye("restoring vertex: ", 'g') << vertex << endl;
    }
    
    //TODO: replace u by k
    //std::cout << "before cannot explore neighbours cutoff" << std::endl;
	//cannot fully explore neighbours
    if (vertexDeg > u)
    {
        G->unreduce(&k, previousK, depth, firstSolution.second);
        return firstSolution;
    }

    //std::cout << "deleting neighbourhood of vertex " << vertex << ": " << std::endl;
    std::vector<int>* neighbours = G->getNeighbours(vertex);
    /* std::cout << cp::dye("branching: choosing neighbours of vertex " + std::to_string(vertex) + ": ", 'b');
    for(int i = 0; i < (int) neighbours->size(); i++)
    {
        std::cout <<  cp::dye(std::to_string(neighbours->at(i)) + ", ", 'b');
    }
    std::cout << std::endl; */
    G->setInactive(neighbours);
    //std::cout << "prebranch k=" << k << std::endl;
	results = vcVertexBranchingConstrained(G, k - neighbours->size(), c + neighbours->size(), u, depth+1, numRec);
    //std::cout << "after neighbour branching" << std::endl;
    if(u > results.first) {
        //k = k - u + results.first;
        //previousK = previousK - u + results.first;
        u = results.first;
    }
    vc = results.second;
    //std::cout << "postbranch k=" << k << '\n';
    //cout << "prec " << '\n';
	if (vc != nullptr)
	{
        //revert changes for multiple executions of the algorithm
        G->setActive(neighbours);
        for(int i = 0; i < (int) neighbours->size(); i++)
        {
            vc->insert(std::pair(neighbours->at(i), true));
        }
        //G->printVC(vc);
        //G->unreduce(&k, previousK, vc); //unreduce needs correct vc für unmerge of deg2rule
        secondSolution = std::pair<int, std::unordered_map<int, bool>*>(u, vc);
	}
	else
	{
		//revert changes to graph
		G->setActive(neighbours);
        //G->unreduce(&k, previousK);
        secondSolution = std::pair<int, std::unordered_map<int, bool>*>(u, nullptr);
	}
    /* cout << "restoring neighbourhood of vertex " << vertex << ": ";
    for(int i = 0; i < (int) neighbours->size(); i++)
    {
        cout << neighbours->at(i) << ", ";
    }
    cout << endl; */
    // free neighbours
    delete neighbours;

    if(firstSolution.second != nullptr && firstSolution.first <= secondSolution.first)
    {
        if(secondSolution.second != nullptr) { delete secondSolution.second; }
        G->unreduce(&k, previousK, depth, firstSolution.second);
        return firstSolution;
    }
    else
    {
        if(firstSolution.second != nullptr) { delete firstSolution.second; }
        G->unreduce(&k, previousK, depth, secondSolution.second);
        return secondSolution;
    }
}


std::unordered_map<int, bool>* vcSolverConstrained(BucketGraph* G, int* numRec, bool printDebug)
{
    std::vector<bool> rulesToApply = std::vector<bool>{true, true, false, true, true, false, false, false, false, false, true, true};
    BucketGraph* graphCopy = G->copy();
    BucketGraph* graphCopyForPre = G->copy();

    int numPreprocessingVCVerticesCopy = 0;
    graphCopyForPre->preprocess(&numPreprocessingVCVerticesCopy, rulesToApply, printDebug);
    numPreprocessingVCVerticesCopy = -numPreprocessingVCVerticesCopy;

    double UPPER_BOUND_TIME_CAP = 0.1; //in seconds

    /* auto startUpper = std::chrono::high_resolution_clock::now();
    int numHeurRecursions = 0;
    std::unordered_map<int, bool>* heuristicVC = maxHeuristicSolver(G, &numHeurRecursions, false, false);
    int upperBound = heuristicVC->size();
    delete heuristicVC;
    //int u = 15000;//1730;//getUpperBound(G, 5);
    auto endUpper = std::chrono::high_resolution_clock::now();
    double upperBoundDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endUpper - startUpper).count() /  1000) / (double) 1000;
 */
    //preprocessing
    int numPreprocessingVCVertices = 0;
	int k = 0;
    auto startPreprocess = std::chrono::high_resolution_clock::now();
    
    G->preprocess(&numPreprocessingVCVertices, rulesToApply, printDebug);
    numPreprocessingVCVertices = -numPreprocessingVCVertices;
    auto endPreprocess = std::chrono::high_resolution_clock::now();
    double preprocessDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endPreprocess - startPreprocess).count() /  1000) / (double) 1000;
    if(printDebug) {
        std::cout << "#Preprocessed Graph to size n=" << G->getNumVertices() << ", m=" << G->getNumEdges() << " in " << preprocessDuration << " seconds" << " (reduced by " << numPreprocessingVCVertices << " vertices)" << std::endl;
    }

    //get upper bound
    auto startUpper = std::chrono::high_resolution_clock::now();
    int u = getUpperBound(graphCopy, graphCopyForPre, numPreprocessingVCVerticesCopy, UPPER_BOUND_TIME_CAP, printDebug);
    auto endUpper = std::chrono::high_resolution_clock::now();
    double upperBoundDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endUpper - startUpper).count() /  1000) / (double) 1000;
    if(printDebug)
    {
        std::cout << "#Calculated upper Bound u=" << u << " in " << upperBoundDuration << " seconds" << std::endl;
    }

    //init packing
    Packing* packingConstraints = new Packing(G->getVertexReferencesSize());

	std::unordered_map<int, bool>* vc;
    auto startBranching = std::chrono::high_resolution_clock::now();
    auto results = vcVertexBranchingConstrained(G, u, 0, u, 0, numRec, packingConstraints, printDebug);
    vc = results.second;

    if(vc == nullptr)
    {
        vc = new std::unordered_map<int, bool>();
    }
    // Add Reduced Vertices to Vertex Cover
    G->unreduce(&k, vc->size()+numPreprocessingVCVertices, -1, vc);
    if(G->getReductionStackSize() > 0)
    {
        throw std::invalid_argument("vcSolverRecursive: reduction rule stack isn't fully popped after final unreduce");
    }
    auto endBranching = std::chrono::high_resolution_clock::now();
    double branchingDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endBranching - startPreprocess).count() /  1000) / (double) 1000;
    if(printDebug)
        std::cout << "#Finished branching in " << branchingDuration << " seconds" << std::endl;
    return vc;
}

std::pair<int, std::unordered_map<int, bool>*> vcVertexBranchingConstrainedB(BucketGraph* G, int k, int c, int u, int depth, int* numRec,
    Packing* constraints = nullptr, bool printDebug = false)
{
    (*numRec)++;
    int previousK = k;
    bool cut = false;
    //G->printReductionStack();
    //std::cout << "preReduce" << std::endl;
    //G->print();
    cut = G->dynamicReduce(&k, depth, printDebug);
    int dataredK = k;
    //std::cout << "postReduce" << std::endl;
    //G->print();
    if(cut)
    {
        //std::cout << "Cutting through data reductions" << std::endl;
        //std::cout << cp::dye("data reduction cutoff unreduce: ", 'y') << std::endl;
        G->unreduce(&k, previousK, depth);
        //std::cout << "postUnreduce" << std::endl;
        //G->print();
        return std::pair<int, std::unordered_map<int, bool>*>(u, nullptr);
    }
    c += previousK - k;

    if (c + G->getLPBound() >= u) {
        //std::cout << "Cutting through at least as good initial solution" << std::endl;
        //std::cout << cp::dye("constrained cutoff unreduce: ", 'y') << std::endl;
        G->unreduce(&k, previousK, depth);
        return std::pair<int, std::unordered_map<int, bool>*>(u, nullptr);
    }

	int vertex = G->getMaxDegreeVertex();
    //no vertices left
    if (vertex == -1)
    {
        //std::cout << "Found through no vertex" << std::endl;
        //std::cout << cp::dye("no vertex unreduce: ", 'y') << vertex << std::endl;
        std::unordered_map<int, bool>* vc = new std::unordered_map<int, bool>();
        G->unreduce(&k, previousK, depth, vc);
        return std::pair<int, std::unordered_map<int, bool>*>(c, vc);
    }
    int vertexDeg = G->getVertexDegree(vertex);
    //G->print();
	//graph has no edges left
	if (vertexDeg == 0)
	{
        //std::cout << "Found through deg0 vertex" << std::endl;
        //std::cout << cp::dye("deg0 unreduce: ", 'y') << vertex << std::endl;
        std::unordered_map<int, bool>* vc = new std::unordered_map<int, bool>();
        G->unreduce(&k, previousK, depth, vc);
        return std::pair<int, std::unordered_map<int, bool>*>(c, vc);
	}

    std::pair<int, std::unordered_map<int, bool>*> firstSolution = std::pair<int, std::unordered_map<int, bool>*>(u, nullptr);
    std::pair<int, std::unordered_map<int, bool>*> secondSolution = std::pair<int, std::unordered_map<int, bool>*>(u, nullptr);

    std::pair<int, std::unordered_map<int, bool>*> results;
    std::unordered_map<int, bool>* vc = nullptr;
    if(constraints->checkNeighbourhoodConstraints(vertex, G))
    {
        //std::cout << "MaxDeg Branching" << std::endl;
        //if(vertex == 33) { std::cout << "Setting vertex 33 inactive during branching" << std::endl; }
        //delete first vertex from graph and explore solution
        G->setInactive(vertex);
        results = vcVertexBranchingConstrainedB(G, u - c/* k - 1 */, c + 1, u, depth+1, numRec);
        //std::cout << "After MaxDeg Branching" << std::endl;
        if(u > results.first) {
            /* k = k - u + results.first;
            previousK = previousK - u + results.first; */
            u = results.first;
            /* u = results.first;
            previousK = previousK - k + u - c;
            k = u - c; */
        }
        vc = results.second;
        //if(vertex == 33) { std::cout << "Setting vertex 33 active after branching" << std::endl; }
        if (vc != nullptr)
        {
            //revert changes for multiple executions of the algorithm
            vc->insert(std::pair(vertex, true));
            G->setActive(vertex);
            //G->unreduce(&k, previousK, vc); //unreduce needs correct vc für unmerge of deg2rule
            //G->printVC(vc);
            firstSolution = std::pair<int, std::unordered_map<int, bool>*>(u, vc);
        }
        else
        {
            //revert changes to graph
            G->setActive(vertex);
            //constraints->undoNeighbourhoodConstraints(vertex, G);
            //G->unreduce(&k, previousK);
            firstSolution = std::pair<int, std::unordered_map<int, bool>*>(u, nullptr);
        }
        //cout << cp::dye("restoring vertex: ", 'g') << vertex << endl;
    }

//TODO: replace u by k
    //std::cout << "before cannot explore neighbours cutoff" << std::endl;
	//cannot fully explore neighbours
    if (vertexDeg > u)
    {
        //std::cout << cp::dye("explore neighbours unreduce: ", 'y') << vertex << std::endl;
        //std::cout << "Cutting through Neighbourhood bigger than upper bound" << std::endl;
        if(firstSolution.second != nullptr) {
            G->unreduce(&k, previousK, depth, firstSolution.second);
        }
        else { 
            /* G->printReductionStack();
            G->print();
            std::cout << std::endl; */
            G->unreduce(&k, previousK, depth); }
        return firstSolution;
    }

    //std::cout << "Neighbourhood Branching" << std::endl;
    //std::cout << "deleting neighbourhood of vertex " << vertex << ": " << std::endl;
    std::vector<int>* neighbours = G->getNeighbours(vertex);
    /* std::cout << cp::dye("branching: choosing neighbours of vertex " + std::to_string(vertex) + ": ", 'b');
    for(int i = 0; i < (int) neighbours->size(); i++)
    {
        std::cout <<  cp::dye(std::to_string(neighbours->at(i)) + ", ", 'b');
    }
    std::cout << std::endl; */
    //if(G->vertexHasEdgeTo(vertex, 33)) { std::cout << "Setting vertex 33 inactive during neighbour branching" << std::endl; }
    G->setInactive(neighbours);
    //std::cout << "prebranch k=" << k << std::endl;
	results = vcVertexBranchingConstrainedB(G, u - c/* k - neighbours->size() */, c + neighbours->size(), u, depth+1, numRec);
    //std::cout << "after neighbour branching" << std::endl;
    if(u > results.first) {
        //k = k - u + results.first;
        //previousK = previousK - u + results.first;
        u = results.first;
    }
    vc = results.second;
    //std::cout << "postbranch k=" << k << '\n';
    //cout << "prec " << '\n';
    //if(G->vertexHasEdgeTo(vertex, 33)) { std::cout << "Setting vertex 33 active after neighbour branching" << std::endl; }
	if (vc != nullptr)
	{
        //revert changes for multiple executions of the algorithm
        G->setActive(neighbours);
        for(int i = 0; i < (int) neighbours->size(); i++)
        {
            vc->insert(std::pair(neighbours->at(i), true));
        }
        //G->printVC(vc);
        //G->unreduce(&k, previousK, vc); //unreduce needs correct vc für unmerge of deg2rule
        secondSolution = std::pair<int, std::unordered_map<int, bool>*>(u, vc);
	}
	else
	{
		//revert changes to graph
		G->setActive(neighbours);
        //G->unreduce(&k, previousK);
        secondSolution = std::pair<int, std::unordered_map<int, bool>*>(u, nullptr);
	}
    //G->print();
    /* cout << "restoring neighbourhood of vertex " << vertex << ": ";
    for(int i = 0; i < (int) neighbours->size(); i++)
    {
        cout << neighbours->at(i) << ", ";
    }
    cout << endl; */
    // free neighbours
    delete neighbours;

    if(firstSolution.second != nullptr && firstSolution.first <= secondSolution.first)
    {
        if(secondSolution.second != nullptr) { delete secondSolution.second; }
        //std::cout << cp::dye("chose first solution unreduce: ", 'y') << vertex << std::endl;
        if(firstSolution.second != nullptr) { G->unreduce(&k, previousK, depth, firstSolution.second); }
        else { G->unreduce(&k, previousK, depth); }
        return firstSolution;
    }
    else
    {
        if(firstSolution.second != nullptr) { delete firstSolution.second; }
        //std::cout << cp::dye("chose secondsolution unreduce: ", 'y') << vertex << std::endl;
        if(secondSolution.second != nullptr) { G->unreduce(&k, previousK, depth, secondSolution.second); }
        else { G->unreduce(&k, previousK, depth); }
        return secondSolution;
    }
}

std::unordered_map<int, bool>* vcSolverConstrainedB(BucketGraph* G, int* numRec, bool printDebug)
{
    std::vector<bool> rulesToApply = std::vector<bool>{true, true, false, true, true, false, false, false, false, false, false, false};
    BucketGraph* graphCopy = G->copy();
    BucketGraph* graphCopyForPre = G->copy();

    int numPreprocessingVCVerticesCopy = 0;
    graphCopyForPre->preprocess(&numPreprocessingVCVerticesCopy, rulesToApply, printDebug);
    numPreprocessingVCVerticesCopy = -numPreprocessingVCVerticesCopy;

    double UPPER_BOUND_TIME_CAP = 0.1;//1.f; //in seconds

    /* auto startUpper = std::chrono::high_resolution_clock::now();
    int numHeurRecursions = 0;
    std::unordered_map<int, bool>* heuristicVC = maxHeuristicSolver(G, &numHeurRecursions, false, false);
    int upperBound = heuristicVC->size();
    delete heuristicVC;
    //int u = 15000;//1730;//getUpperBound(G, 5);
    auto endUpper = std::chrono::high_resolution_clock::now();
    double upperBoundDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endUpper - startUpper).count() /  1000) / (double) 1000;
 */
    //preprocessing
    int numPreprocessingVCVertices = 0;
	int k = 0;
    auto startPreprocess = std::chrono::high_resolution_clock::now();

    G->preprocess(&numPreprocessingVCVertices, rulesToApply, printDebug);
    numPreprocessingVCVertices = -numPreprocessingVCVertices;
    auto endPreprocess = std::chrono::high_resolution_clock::now();
    double preprocessDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endPreprocess - startPreprocess).count() /  1000) / (double) 1000;
    if(printDebug) {
        std::cout << "#Preprocessed Graph to size n=" << G->getNumVertices() << ", m=" << G->getNumEdges() << " in " << preprocessDuration << " seconds" << " (reduced by " << numPreprocessingVCVertices << " vertices)" << std::endl;
    }

    //get upper bound
    auto startUpper = std::chrono::high_resolution_clock::now();
    std::unordered_map<int, bool>* heuristicVC = getUpperBoundVC(graphCopy, graphCopyForPre, numPreprocessingVCVerticesCopy, UPPER_BOUND_TIME_CAP, printDebug);
    int u = heuristicVC->size();
    auto endUpper = std::chrono::high_resolution_clock::now();
    double upperBoundDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endUpper - startUpper).count() /  1000) / (double) 1000;
    if(printDebug)
    {
        std::cout << "#Calculated upper Bound u=" << u << " in " << upperBoundDuration << " seconds" << std::endl;
    }

    //G->print();
    /* G->printReductionStack();
    std::cout << std::endl; */

    //init packing
    Packing* packingConstraints = new Packing(G->getVertexReferencesSize());

	std::unordered_map<int, bool>* vc;
    auto startBranching = std::chrono::high_resolution_clock::now();
    auto results = vcVertexBranchingConstrainedB(G, u, 0, u, 0, numRec, packingConstraints, true || printDebug);
    vc = results.second;

    bool wasNullptr = false;
    if(vc == nullptr)
    {
        wasNullptr = true;
        vc = new std::unordered_map<int, bool>();
    }
    // Add Reduced Vertices to Vertex Cover
    G->unreduce(&k, vc->size()+numPreprocessingVCVertices, -1, vc);
    if(G->getReductionStackSize() > 0)
    {
        throw std::invalid_argument("vcSolverRecursive: reduction rule stack isn't fully popped after final unreduce");
    }
    auto endBranching = std::chrono::high_resolution_clock::now();
    double branchingDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endBranching - startPreprocess).count() /  1000) / (double) 1000;
    if(printDebug)
        std::cout << "#Finished branching in " << branchingDuration << " seconds" << std::endl;
    if(vc->size() == 0 && wasNullptr)
    {
        // TODO: if no solution found
        return heuristicVC;
    }
    return vc;
}

std::pair<int, std::unordered_map<int, bool>*> vcVertexBranchingConstrainedEloi(BucketGraph* G, int c, int u, int depth, std::unordered_map<int, bool>* u_solution, 
 int* numRec, bool printDebug = false)
{
    if(printDebug)
        std::cout << "#-->Branching: c=" << c << ", u=" << u << ", depth=" << depth << std::endl;

    (*numRec)++;
    int k = u - c;
    int previousK = k;
    bool cut = false;
    cut = G->dynamicReduce(&k, depth, printDebug);

    if(printDebug)
        std::cout << "#After reduce: kDecrement=" << previousK - k << std::endl;

    //insufficient budget
    if(cut)
    {
        if(printDebug)
            std::cout << "#reduction cut: insufficient budget"<< std::endl;

        G->unreduce(&k, previousK, depth);
        return std::pair<int, std::unordered_map<int, bool>*>(u, u_solution);
    }
    c += previousK - k;

    int lowerBound = G->getLPBound();
    if(printDebug)
        std::cout << "#Before upper bound cut: lp=" << lowerBound << ", c=" << c << std::endl;

    //cannot find better solution in this branch ---> cut
    if (c + lowerBound >= u) {
        if(printDebug)
            std::cout << "#upper bound cut"<< std::endl;
        G->unreduce(&k, previousK, depth);
        return std::pair<int, std::unordered_map<int, bool>*>(u, u_solution);
    }

    int vertex = G->getMaxDegreeVertex();
    //no vertices left
    if (vertex == -1)
    {
        if(printDebug)
            std::cout << "#no vertices left, return"<< std::endl;

        std::unordered_map<int, bool>* vc = new std::unordered_map<int, bool>();
        G->unreduce(&k, previousK, depth, vc);
        return std::pair<int, std::unordered_map<int, bool>*>(c, vc);
    }
    int vertexDeg = G->getVertexDegree(vertex);

	//graph has no edges left
	if (vertexDeg == 0)
	{
        if(printDebug)
            std::cout << "#no edges left, return"<< std::endl;

        std::unordered_map<int, bool>* vc = new std::unordered_map<int, bool>();
        G->unreduce(&k, previousK, depth, vc);
        return std::pair<int, std::unordered_map<int, bool>*>(c, vc);
	}

    if(printDebug)
        std::cout << "#Before branching on vertex=" << vertex << std::endl;

    std::pair<int, std::unordered_map<int, bool>*> firstSolution = std::pair<int, std::unordered_map<int, bool>*>(u, nullptr);
    std::pair<int, std::unordered_map<int, bool>*> secondSolution = std::pair<int, std::unordered_map<int, bool>*>(u, nullptr);

    G->setInactive(vertex);
    firstSolution = vcVertexBranchingConstrainedEloi(G, c + 1, u, depth+1, u_solution, numRec, printDebug);
    bool firstWasImprovement = false;

    if(printDebug)
            std::cout << "#After branching on vertex=" << vertex << std::endl;
    //found better solution, construct it in back propagation
    if(firstSolution.first < u)
    {
        firstSolution.second->insert(std::pair(vertex, true));
        //u = firstSolution.first; //TODO: reenable
        firstWasImprovement = true;
        G->setActive(vertex);
    }
    else //just revert graph
    {
        G->setActive(vertex);
    }

    if(printDebug)
        std::cout << "#before cannot explore neighbours cutoff" << std::endl;
	//no budget, cannot fully explore neighbours
    /* if (vertexDeg > k)
    {
        if(firstWasImprovement)
        {
            G->unreduce(&k, previousK, depth, firstSolution.second);
            return firstSolution;
        }
        else
        {
            G->unreduce(&k, previousK, depth);
            return std::pair<int, std::unordered_map<int, bool>*>(u, u_solution);
        }
    } */

    if(printDebug)
        std::cout << "#Before branching on neighbours of v=" << vertex << std::endl;
    //std::cout << "deleting neighbourhood of vertex " << vertex << ": " << std::endl;

    std::vector<int>* neighbours = G->getNeighbours(vertex);
    G->setInactive(neighbours);
	secondSolution = vcVertexBranchingConstrainedEloi(G, c + neighbours->size(), u, depth+1, u_solution, numRec, printDebug);
    bool secondWasImprovement = false;

    //found better solution, construct it in back propagation
	if (secondSolution.first < u)
	{
        for(int i = 0; i < (int) neighbours->size(); i++)
        {
            secondSolution.second->insert(std::pair(neighbours->at(i), true));
        }
        u = secondSolution.first;
        secondWasImprovement = true;
        //revert changes for multiple executions of the algorithm
        G->setActive(neighbours);
	}
	else
	{
		//revert changes to graph
		G->setActive(neighbours);
	}

    if(printDebug)
        std::cout << "#After branching on neighbours of vertex=" << vertex << std::endl;

    // free neighbours
    delete neighbours;

    //decide which solution to pass back
    if(firstSolution.first <= secondSolution.first)
    {
        if(printDebug)
            std::cout << "#Passing back first solution with u=" << firstSolution.first << std::endl;
        if(firstSolution.second != nullptr && secondWasImprovement) { delete secondSolution.second; }
        G->unreduce(&k, previousK, depth, firstSolution.second);
        return firstSolution;
    }
    else
    {
        if(printDebug)
            std::cout << "#Passing back second solution with u=" << secondSolution.first << std::endl;
        if(firstSolution.second != nullptr && firstWasImprovement) { delete firstSolution.second; }
        G->unreduce(&k, previousK, depth, secondSolution.second);
        return secondSolution;
    }
}


std::unordered_map<int, bool>* vcSolverConstrainedEloi(BucketGraph* G, int* numRec, bool printDebug)
{
    double UPPER_BOUND_TIME_CAP = 0.1; //in seconds
        //                                                  0     1      2      3    4      5      6      7      8      9      10     11
    //std::vector<bool> rulesToApply = std::vector<bool>{true, true, false, true, true, false, false, false, false, false, false, false};
    std::vector<bool> rulesToApply = std::vector<bool>{true, false, true, true, true, false, false, true, true, true, false, false};

    //BucketGraph* graphCopy = G->copy();
    //BucketGraph* graphCopyWithPre = G->copy();

    //preprocessing for heuristic graph
    /* int numPreprocessingVCVerticesCopy = 0;
    graphCopyWithPre->preprocess(&numPreprocessingVCVerticesCopy, rulesToApply, printDebug);
    numPreprocessingVCVerticesCopy = -numPreprocessingVCVerticesCopy; */

    //preprocessing for branch graph
    int numPreprocessingVCVertices = 0;
	int k = 0;
    auto startPreprocess = std::chrono::high_resolution_clock::now();

    G->preprocess(&numPreprocessingVCVertices, rulesToApply, printDebug);
    numPreprocessingVCVertices = -numPreprocessingVCVertices;
    auto endPreprocess = std::chrono::high_resolution_clock::now();
    double preprocessDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endPreprocess - startPreprocess).count() /  1000) / (double) 1000;
    if(printDebug) {
        std::cout << "#Preprocessed Graph to size n=" << G->getNumVertices() << ", m=" << G->getNumEdges() << " in " << preprocessDuration << " seconds" << " (reduced by " << numPreprocessingVCVertices << " vertices)" << std::endl;
    }

    //get upper bound
    //first generate a fast heuristic solution
    int heuristicNumRecursions = 0;
    std::unordered_map<int, bool>* heuristicVC = maxHeuristicSolver(G, &heuristicNumRecursions, false, false);
    int u = heuristicVC->size();
    if(printDebug)
    {
        std::cout << "#Calculated upper Bound u=" << u << std::endl;
    }

/*  TODO: doesnt work yet because need initial solution, not just bound
    auto startUpper = std::chrono::high_resolution_clock::now();
    int u = getUpperBound(graphCopy, graphCopyWithPre, numPreprocessingVCVerticesCopy, UPPER_BOUND_TIME_CAP, printDebug);
    auto endUpper = std::chrono::high_resolution_clock::now();
    double upperBoundDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endUpper - startUpper).count() /  1000) / (double) 1000;
    if(printDebug)
    {
        std::cout << "#Calculated upper Bound u=" << u << " in " << upperBoundDuration << " seconds" << std::endl;
    } */


    auto startBranching = std::chrono::high_resolution_clock::now();
    auto results = vcVertexBranchingConstrainedEloi(G, 0, u, 0, heuristicVC, numRec, false);
    std::unordered_map<int, bool>* vc = results.second;

    if(printDebug)
        std::cout << "#Finished branching with u = " << u << ", vc without pre = " << results.first << std::endl;

    // Add Reduced Vertices to Vertex Cover
    G->unreduce(&k, vc->size()+numPreprocessingVCVertices, -1, vc);
    if(G->getReductionStackSize() > 0)
    {
        throw std::invalid_argument("vcSolverRecursive: reduction rule stack isn't fully popped after final unreduce");
    }
    auto endBranching = std::chrono::high_resolution_clock::now();
    double branchingDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endBranching - startPreprocess).count() /  1000) / (double) 1000;
    if(printDebug)
        std::cout << "#Finished branching and unreducing preprocessing in " << branchingDuration << " seconds, u = " << u << ", vc without pre = " << results.first << ", vc size: " << vc->size() << std::endl;
    return vc;
}

int determineOptimalSolutionSizeBranching(BucketGraph* G, int c, int u, int depth, int* numRec, bool printDebug = false)
{
    if(printDebug)
        std::cout << "#-->Branching: c=" << c << ", u=" << u << ", depth=" << depth << std::endl;

    (*numRec)++;
    int k = u - c;
    int previousK = k;
    bool cut = false;
    cut = G->dynamicReduce(&k, depth, printDebug);
    int numReduced = previousK - k;

    if(printDebug)
        std::cout << "#After reduce: kDecrement=" << previousK - k << std::endl;

    //insufficient budget
    if(cut)
    {
        if(printDebug)
            std::cout << "#reduction cut: insufficient budget"<< std::endl;

        G->unreduce(&k, previousK, depth);
        return u;
    }
    c += numReduced;

    int lowerBound = G->getLPBound();
    if(printDebug)
        std::cout << "#Before upper bound cut: lp=" << lowerBound << ", c=" << c << std::endl;

    //cannot find better solution in this branch ---> cut
    if (c + lowerBound >= u) {
        if(printDebug)
            std::cout << "#upper bound cut"<< std::endl;
        G->unreduce(&k, previousK, depth);
        return u;
    }

    int vertex = G->getMaxDegreeVertex();
    //no vertices left
    if (vertex == -1)
    {
        if(printDebug)
            std::cout << "#no vertices left, return"<< std::endl;

        G->unreduce(&k, previousK, depth);
        return c;
    }
    int vertexDeg = G->getVertexDegree(vertex);

	//graph has no edges left
	if (vertexDeg == 0)
	{
        if(printDebug)
            std::cout << "#no edges left, return"<< std::endl;

        G->unreduce(&k, previousK, depth);
        return c;
	}

    if(printDebug)
        std::cout << "#Before branching on vertex=" << vertex << std::endl;

    G->setInactive(vertex);
    int firstU = determineOptimalSolutionSizeBranching(G, c + 1, u, depth+1, numRec, printDebug);

    if(printDebug)
            std::cout << "#After branching on vertex=" << vertex << std::endl;
    //found better solution, construct it in back propagation
    if(firstU < u)
    {
        u = firstU;
        G->setActive(vertex);
    }
    else //just revert graph
    {
        G->setActive(vertex);
    }

    if(printDebug)
        std::cout << "#before cannot explore neighbours cutoff" << std::endl;
	//no budget, cannot fully explore neighbours
    /* if (vertexDeg > k)
    {
        if(firstWasImprovement)
        {
            G->unreduce(&k, previousK, depth, firstSolution.second);
            return firstSolution;
        }
        else
        {
            G->unreduce(&k, previousK, depth);
            return std::pair<int, std::unordered_map<int, bool>*>(u, u_solution);
        }
    } */

    if(printDebug)
        std::cout << "#Before branching on neighbours of v=" << vertex << std::endl;
    //std::cout << "deleting neighbourhood of vertex " << vertex << ": " << std::endl;

    std::vector<int>* neighbours = G->getNeighbours(vertex);
    G->setInactive(neighbours);
	int secondU = determineOptimalSolutionSizeBranching(G, c + neighbours->size(), u, depth+1, numRec, printDebug);

    //found better solution, construct it in back propagation
	if (secondU < u)
	{
        u = secondU;
        G->setActive(neighbours);
	}
	else
	{
		//revert changes to graph
		G->setActive(neighbours);
	}

    if(printDebug)
        std::cout << "#After branching on neighbours of vertex=" << vertex << std::endl;

    // free neighbours
    delete neighbours;

    G->unreduce(&k, previousK, depth);
    return u + numReduced;
}

int determineOptimalSolutionSize(BucketGraph* G, int* numRec, bool printDebug)
{
    double UPPER_BOUND_TIME_CAP = 0.1; //in seconds
        //                                                  0     1      2      3    4      5      6      7      8      9      10     11
    //std::vector<bool> rulesToApply = std::vector<bool>{true, true, false, true, true, false, false, false, false, false, false, false};
    std::vector<bool> rulesToApply = std::vector<bool>{true, false, false, true, true, false, false, true, true, true, false, false};

    //BucketGraph* graphCopy = G->copy();
    //BucketGraph* graphCopyWithPre = G->copy();

    //preprocessing for heuristic graph
    /* int numPreprocessingVCVerticesCopy = 0;
    graphCopyWithPre->preprocess(&numPreprocessingVCVerticesCopy, rulesToApply, printDebug);
    numPreprocessingVCVerticesCopy = -numPreprocessingVCVerticesCopy; */

    //preprocessing for branch graph
    int numPreprocessingVCVertices = 0;
	int k = 0;
    auto startPreprocess = std::chrono::high_resolution_clock::now();

    G->preprocess(&numPreprocessingVCVertices, rulesToApply, printDebug);
    numPreprocessingVCVertices = -numPreprocessingVCVertices;
    auto endPreprocess = std::chrono::high_resolution_clock::now();
    double preprocessDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endPreprocess - startPreprocess).count() /  1000) / (double) 1000;
    if(printDebug) {
        std::cout << "#Preprocessed Graph to size n=" << G->getNumVertices() << ", m=" << G->getNumEdges() << " in " << preprocessDuration << " seconds" << " (reduced by " << numPreprocessingVCVertices << " vertices)" << std::endl;
    }

    //get upper bound
    //first generate a fast heuristic solution
    int heuristicNumRecursions = 0;
    std::unordered_map<int, bool>* heuristicVC = maxHeuristicSolver(G, &heuristicNumRecursions, false, false);
    int u = heuristicVC->size();
    if(printDebug)
    {
        std::cout << "#Calculated upper Bound u=" << u << std::endl;
    }

/*  TODO: doesnt work yet because need initial solution, not just bound
    auto startUpper = std::chrono::high_resolution_clock::now();
    int u = getUpperBound(graphCopy, graphCopyWithPre, numPreprocessingVCVerticesCopy, UPPER_BOUND_TIME_CAP, printDebug);
    auto endUpper = std::chrono::high_resolution_clock::now();
    double upperBoundDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endUpper - startUpper).count() /  1000) / (double) 1000;
    if(printDebug)
    {
        std::cout << "#Calculated upper Bound u=" << u << " in " << upperBoundDuration << " seconds" << std::endl;
    } */


    auto startBranching = std::chrono::high_resolution_clock::now();
    int result = determineOptimalSolutionSizeBranching(G, 0, u, 0, numRec, false);

    if(printDebug)
        std::cout << "#Finished branching with u = " << u << ", vc without pre = " << result << std::endl;

    // Add Reduced Vertices to Vertex Cover
    G->unreduce(&k, result+numPreprocessingVCVertices, -1);
    if(G->getReductionStackSize() > 0)
    {
        throw std::invalid_argument("vcSolverRecursive: reduction rule stack isn't fully popped after final unreduce");
    }
    auto endBranching = std::chrono::high_resolution_clock::now();
    double branchingDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endBranching - startPreprocess).count() /  1000) / (double) 1000;
    if(printDebug)
        std::cout << "#Finished branching and unreducing preprocessing in " << branchingDuration << " seconds, u = " << u << ", vc without pre = " << result << ", vc size: " << result + numPreprocessingVCVertices << std::endl;
    return result + numPreprocessingVCVertices;
}