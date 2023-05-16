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

vector<int>* vcVertexBranchingRecursive(ArrayGraph* G, int k)
{
	
	if (k < 0)
		return nullptr;

	int vertex = G->getFirstActiveVertex(); //TODO:get max degree vertex
    int vertexDeg = G->getVertexDegree(vertex);

	//graph has no edges left
	if (vertexDeg == 0)
	{
		return new vector<int>();
	}

	//delete first vertex from graph and explore solution
    G->setInactive(vertex);
	vector<int>* S = vcVertexBranchingRecursive(G, k - 1);
	if (S != nullptr)
	{
		//revert changes to graph
		G->setActive(vertex); //TODO: not necessary????
		//return results
		S->push_back(vertex);
		return S;
	}
	else
	{
		//revert changes to graph
		G->setActive(vertex);
	}


	//cannot fully explore neighbours
    if (vertexDeg > k) 
    {
        return nullptr;
    }
    vector<int>* neighbours = G->getNeighbours(vertex);
    G->setInactive(neighbours);
	S = vcVertexBranchingRecursive(G, k - 1);
	if (S != nullptr)
	{
		//revert changes to graph
		G->setActive(neighbours); //TODO: not necessary????
		//return results
        for (int i = 0; i < (int) neighbours->size(); i++)
        {
            S->push_back(neighbours->at(i));
        }
		return S;
	}
	else
	{
		//revert changes to graph
		G->setInactive(neighbours);
	}
	return nullptr;
}

vector<int>* vertexBranchingSolverRecursive(ArrayGraph* G)
{
	int k = G->getVCLowerBound();
	vector<int> *vc;

	while (true)
	{
		vc = vcVertexBranchingRecursive(G, k);
		if (vc != nullptr)
		{
			return vc;
		}
		k++;
	}
}

vector<int>* vcVertexBranchingIterative(ArrayGraph* G, int k, std::vector<int>* vc)
{
	std::cout << "setting up BnB data structures\n";
	// stack storing the differentials of partial solutions, currently under evaluation and whether a partial solution was already expanded
    std::stack<std::pair<std::vector<int>*, bool>*> S;
	std::pair<std::vector<int>*, bool>* current;
	int branchVertex;
	// number of currently active vertices
	int partialVCSize = 0;
	
	// TODO: different initialization
	// first branching done by hand
	// get maxDegVert of remaining active vertices (the way the algorithm branches is determined by the choice of this vertex)
	/* branchVertex = G->getMaxDegreeVertex();
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
	} */
	std::vector<int> bv = {};
	std::pair<std::vector<int>*, bool> childDifferential = {&bv, false};
	S.push(&childDifferential);
	std::cout << "started BnB with stack size: " << S.size() << "\n";

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
		std::cout << "element was: " << S.size() << "\n";

		// otherwise delete the vertices (defined by the differential calculated in the expansion of the parent search tree node)
		G->setInactive(current->first);
		partialVCSize += current->first->size();
		// get maxDegVert of remaining active vertices (the way the algorithm branches is determined by the choice of this vertex)
		branchVertex = G->getMaxDegreeVertex();
		// mark this partial solution as expanded and expand subsequently
		current->second = true;

		// if no active vertex left in graph or no vertex with degree >= 0: (We found a solution)
		if (branchVertex == -1 || G->getVertexDegree(branchVertex) == 0)
		{
			// if no 
			// or if current solution is actually better than the current best solution: update k & vc
			if (vc == nullptr || k > partialVCSize) { // TODO: remove condition later when culling earlier with BnB
				vc = G->getInactiveVertices();
				std::cout << "calculated vc of size: " << vc->size() << "\n";
				k = partialVCSize;
			}
			// traverse back up the search tree
			G->setActive(current->first);
			partialVCSize -= current->first->size();
			S.pop();
			continue;
		}

		// solve graph with maxVertDegree <= 2 in linear time
		else if (G->getVertexDegree(branchVertex) <= 2 && false) // TODO: rm false
		{

		}

		// refined search tree branching for maxVertDegree >= 3
		else if (G->getVertexDegree(branchVertex) >= 1 /*3*/) // TODO: readd 3
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

vector<int>* vertexBranchingSolverIterative(ArrayGraph* G)
{
	int k = G->getVCLowerBound();
	vector<int>* vc = new vector<int>();

	while (true)
	{
		vcVertexBranchingIterative(G, k, vc);
		if(vc->size() == 0) { std::cout << "Did not find solution for k=" << k << "\n"; }
		if (vc->size() > 0)
		{
			return vc;
		}
		k++;
	}
}

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

void writeSolutionToConsole(vector<int>* vc)
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
		if (G == nullptr)
		{
			cerr << "Error constructing graph from input file.";
		}
        G->print();
		//test vc solver
		//vector<string>* vc = searchTreeSolve(G);
		//writeSolutionToConsole(vc);
		std::vector<int>* vc = vertexBranchingSolverRecursive(G);
        cout << vc->size() << std::endl;
		//writeSolutionToConsole(vc);
	}
	catch (const exception& e)
	{
		cerr << "Error launching vertex cover solver.\n";
        cerr << e.what();
	}
}

