#include <iostream>
#include <unordered_map> //O(1) for insert and access instead of O(log n) for ordered maps
#include <fstream>  //ifstream file opening
#include <stack>          // std::stack
#include "Graph.h"
#include "ArrayGraph.h"

using namespace std;

/*----------------------------------------------------------*/
/*---------------   Exercise 2 Solver Code   ---------------*/
/*----------------------------------------------------------*/

/* vector<int>* searchTreeSolveBNB(ArrayGraph* G)
{
	// current best number of vertices for VC
	int k = G->getLowerBoundVC();
	std::vector<int>* vc;
	std::vector<std::pair<bool, int>>* state = G->getState();

	// stack storing the current partial solutions left for evaluation
    std::stack<std::vector<int>*> S;
	std::vector<int>* partialVC;
	int current;

	S.push({});
	while (!S.empty())
	{
		// retrieve the current partial vertex cover vertices
		partialVC = S.top();
		S.pop();
		// set Graph to current partial solution
		G->setInactiveVertices(partialVC);
		// get maxDegVert of remaining active vertices
		current = G->getMaxDegreeVertex();

		// if no active vertex left in graph or no vertex with degree >= 1: (We found a solution)
		if (current == -1 || current == 0)
		{
			// if current solution is actually better than the current best solution: update k & vc
			if (k > (int) partialVC->size()) { // TODO: remove condition later when culling earlier with BnB
				vc = partialVC;
				k = vc->size();
			}
			// traverse back up the search tree
			continue;
		}

		// solve graph with maxVertDegree <= 2 in linear time
		else if (state->at(current).second <= 2)	// FIXME: state[current][1] needs to be updated on the graph side for this to be more accurate
		{

		}

		// refined search tree branching for maxVertDegree >= 3
		else if (state->at(current).second >= 3)	// FIXME: state[current][1] needs to be updated on the graph side for this to be more accurate
		{
			// if k and current partial VC size permit adding the neighbours
			if (k - (int) partialVC->size() >= state->at(current).second)
			{
				// add neighbours of the current vertex to the partial VC and push the resulting partial VC for further evaluation
				// NOTE: partial VC is copied, since we need two separate lists, one for each of the possible branching solutions
				std::vector<int> partialVCCopy = *partialVC;
				for(int neighbour : *(G->getNeighbours(current)))
				{
					partialVCCopy.push_back(neighbour);
				}
				S.push(&partialVCCopy);
			}
			// if k and current partial VC size permit adding the current vertex
			if (k - partialVC->size() >= 1)
			{
				// add current vertex to the partial VC and push the resulting partial VC for further evaluation
				partialVC->push_back(current);
				S.push(partialVC);
			}
		}
	}


	return vc;
} */

vector<int>* searchTreeSolveDifferentialBNB(ArrayGraph* G)
{
	// current best number of vertices for VC
	int k = G->getLowerBoundVC();
	std::vector<int>* vc;
	std::vector<std::pair<bool, int>>* state = G->getState();

	// stack storing the differentials of partial solutions, currently under evaluation and whether a partial solution was already expanded
    std::stack<std::pair<std::vector<int>*, bool>*> S;
	std::pair<std::vector<int>*, bool>* current;
	std::vector<int>* VC = nullptr;
	int branchVertex;
	// number of currently active vertices
	int partialVCSize = 0; // TODO: update

	//S.push({{}, false}); // TODO: different initialization
	while (!S.empty())
	{
		// retrieve the differential of the current partial vertex cover solution to its parent solution (+ expanded tag)
		current = S.top();

		// if current solution has already been expanded, revert its deletion of vertices and traverse back up to the next partial solution
		if (current->second)
		{
			G->setActive(current->first);
			partialVCSize -= current->first->size();
			S.pop();
			continue;
		}

		// otherwise delete the vertices (defined by the differential calculated in the expansion of the parent search tree node)
		G->setInactive(current->first);
		// get maxDegVert of remaining active vertices (the way the algorithm branches is determined by the choice of this vertex)
		branchVertex = G->getMaxDegreeVertex();
		// mark this partial solution as expanded and expand subsequently
		current->second = true;

		// if no active vertex left in graph or no vertex with degree >= 0: (We found a solution)
		if (branchVertex == -1 || G->getVertexDegree(branchVertex) == 0)
		{
			// if no 
			// or if current solution is actually better than the current best solution: update k & vc
			if (vc != nullptr || k > partialVCSize) { // TODO: remove condition later when culling earlier with BnB
				vc = G->getInactiveVertices();
				k = partialVCSize;
			}
			// traverse back up the search tree
			G->setActive(current->first);
			partialVCSize -= current->first->size();
			S.pop();
			continue;
		}

		// solve graph with maxVertDegree <= 2 in linear time
		else if (G->getVertexDegree(branchVertex) <= 2)
		{

		}

		// refined search tree branching for maxVertDegree >= 3
		else if (G->getVertexDegree(branchVertex) >= 3)
		{
			// if k and current partial VC size permit adding the neighbours
			if (k - partialVCSize >= G->getVertexDegree(branchVertex))
			{
				// add neighbours of the current vertex to the child differential for evaluating the partial vertex cover where all of the branchVertex's neighbours where taken into the vertex cover
				// TODO: is there no way to create this inline?
				std::pair<std::vector<int>*, bool> childDifferential = {G->getNeighbours(branchVertex), false};
				S.push(&childDifferential);
			}
			// if k and current partial VC size permit adding the current vertex
			if (k - partialVCSize >= 1)
			{
				// TODO: is there no way to create this inline?
				std::vector<int> bv = {branchVertex};
				std::pair<std::vector<int>*, bool> childDifferential = {&bv, false};
				S.push(&childDifferential);
			}
		} 
	}


	return vc;
}

// TODO: calculate minimal number of changes in activity of vertices to switch from the origin subgraph to the destination subgraph
/* std::vector<int, bool>* calcFlagSets(std::vector<int>* origin, std::vector<int>* dest)
{
	std::vector<int, bool> flagSets;
	for(int i=0; i<dest->size(); i++)
	{
		if(dest.) {

		} else {
			flagSets.push_back({});
		}
	}
	return &flagSets;
} */

/*----------------------------------------------------------*/
/*---------------   Exercise 1 Solver Code   ---------------*/
/*----------------------------------------------------------*/

void deleteStringFromVector(vector<string>* vec, string element)
{
	auto position = std::find(vec->begin(), vec->end(), element);
	if (position != vec->end()) // == myVector.end() means the element was not found
		vec->erase(position);
	return;
}

vector<string>* vcBranch(Graph* G, Graph* graphCopy, int k)
{
	
	if (k < 0)
		return NULL;

	pair<string, string>* edge = graphCopy->getFirstValidEdge();

	//graph has no edges left
	if (edge == nullptr)
	{
		return new vector<string>();
	}

	//delete first vertex from graph and explore solution
	graphCopy->deleteAdjacencyMapEntry(edge->first);
	vector<string>* S = vcBranch(G, graphCopy, k - 1);
	if (S != NULL)
	{
		//revert changes to graph
		graphCopy->addAdjacencyMapEntry(G, edge->first); //revert changes to graph
		//return results
		S->push_back(edge->first);
		return S;
	}
	else
	{
		//revert changes to graph
		graphCopy->addAdjacencyMapEntry(G, edge->first);
	}


	//delete second vertex from graph and explore solution
	graphCopy->deleteAdjacencyMapEntry(edge->second);
	S = vcBranch(G, graphCopy, k - 1);
	if (S != NULL)
	{
		//revert changes to graph
		graphCopy->addAdjacencyMapEntry(G, edge->second);
		//return results
		S->push_back(edge->second);
		return S;
	}
	else
	{
		//revert changes to graph
		graphCopy->addAdjacencyMapEntry(G, edge->second); 
	}
	return NULL;
	
	
}

vector<string>* searchTreeSolve(Graph* G)
{
	int k = 0;
	vector<string> *vc;
	Graph* graphCopy = G->copy();

	while (true)
	{
		vc = vcBranch(G, graphCopy, k);
		if (vc != NULL)
		{
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
		outfile << *it << endl;
	}
	outfile.close();
}

void writeSolutionToConsole(vector<string>* vc)
{
	for (auto it = vc->begin(); it != vc->end(); ++it)
	{
		cout << *it << endl;
	}
}

/*----------------------------------------------------------*/
/*-----------------------   Main   -------------------------*/
/*----------------------------------------------------------*/

int main(int argc, char* argv[]) {
	string input;
    cout << "Test\n";
	try
	{
		ArrayGraph* G = ArrayGraph::readStandardInput();
		if (G == NULL)
		{
			cerr << "Error constructing graph from input file.";
		}
        G->print();
		//test vc solver
		//vector<string>* vc = searchTreeSolve(G);
		//writeSolutionToConsole(vc);
	}
	catch (const exception& e)
	{
		cerr << "Error launching vertex cover solver.\n";
        cerr << e.what();
	}

}

