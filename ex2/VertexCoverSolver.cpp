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

vector<int>* searchTreeSolveBNB(ArrayGraph* G)
{
	// current best number of vertices for VC
	int k = G->getLowerBoundVC();
	std::vector<int> *vc;
	std::vector<bool, int> *state = G->getState();

	// stack storing the current partial solutions left for evaluation
    std::stack<std::vector<int>*> S;
	int current;

	// TODO: add k condition, rm multiple invocation of getMaxVertex // TODO: rm first TODO
	S.push({});
	while (!S.empty())
	{
		// set Graph to current partial solution
		G->setInactiveVertices(S.top());
		// get maxDegVert of remaining active vertices
		current = G->getMaxDegreeVertex();

		// if no active vertex left in graph or no vertex with degree >= 1: (We found a solution)
		if (current == -1 || current == 0)
		{
			// if current solution is actually better than the current best solution: update k & vc
			if (k > (S.top())->size()) { // TODO: remove condition later when culling earlier with BnB
				vc = S.top();
				k = vc->size();
			}
			// traverse back up the search tree
			S.pop();
			continue;
		}
		else
		{
			// TODO: branching
		}

        /* if (vc->size() >= k)
        {
            break;
        } */
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
	try
	{
		Graph* G = Graph::readStandardInput();
		if (G == NULL)
		{
			cerr << "Error constructing graph from input file.";
		}
		//test vc solver
		vector<string>* vc = searchTreeSolve(G);
		writeSolutionToConsole(vc);
	}
	catch (const exception&)
	{
		cerr << "Error launching vertex cover solver.";
	}
}

