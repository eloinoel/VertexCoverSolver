#include <iostream>
#include <unordered_map> //O(1) for insert and access instead of O(log n) for ordered maps
#include <fstream>  //ifstream file opening
#include <stack>          // std::stack
#include <math.h>          // INFINITY
#include "Graph.h"
#include "ArrayGraph.h"
#include "ColorPrint.h"

using namespace std;

/*----------------------------------------------------------*/
/*---------------   Exercise 2 Solver Code   ---------------*/
/*----------------------------------------------------------*/

string tileStr(string toTile, int n) {
	string tiling = "";
	for (int i=0; i<n; i++)
	{
		tiling += toTile;
	}
	return tiling;
}


vector<int>* vcVertexBranchingRecursive(ArrayGraph* G, int k)
{
	if (k < 0)
		return nullptr;

	int vertex = G->getMaxDegreeVertex();
    //no vertices left
    if (vertex == -1)
    {
        return new vector<int>();
    }

    int vertexDeg = G->getVertexDegree(vertex); //TODO: maybe duplicate degree calculation
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
	int k = G->getLowerBoundVC();
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

vector<int>* VCVertexBranchingIterative(ArrayGraph* G, int k, std::vector<int>* vc, bool debug)
{
	// stack storing the differentials of partial solutions, currently under evaluation and whether a partial solution was already expanded
    std::stack<std::pair<std::vector<int>, bool>> S = std::stack<std::pair<std::vector<int>, bool>>();
	std::pair<std::vector<int>, bool> current;
	int branchVertex;
	int branchVertexDegree;
	// number of currently active vertices
	int partialVCSize = 0;
	
	// initialize stack
	std::vector<int> bv = {};
	std::pair<std::vector<int>, bool> childDifferential = {bv, false};
	S.push(childDifferential);
	if(debug)
	{
		std::cout << dye("started", 'g') << " BnB with stack size: " << S.size() << " and k=" << k << "\n";
	}

	while (!S.empty())
	{
		// retrieve the differential of the current partial vertex cover solution to its parent solution (+ expanded tag)
		current = S.top();

		if (current.second)
		{
			if(debug)
			{
				std::cout << tileStr("--", partialVCSize) << "- " << dye("Traversing", 'y') << " back up a previously expanded search tree node\n";
				std::cout << tileStr("--", partialVCSize) << "- " << dye("Restoring", 'p') << " vertices: {";
				if (current.first.size() > 0) std::cout << current.first.at(0);
				for (int i=1; i < (int) current.first.size(); i++)
				{
					std::cout << ", " << current.first.at(i);
				}
				std::cout << "}\n";
			}
			G->setActive(&current.first);
			partialVCSize -= current.first.size();
			S.pop();
			continue;
		}

		// otherwise delete the vertices (defined by the differential calculated in the expansion of the parent search tree node)
		G->setInactive(&current.first);
		partialVCSize += current.first.size();

		if(debug)
		{
			auto SP = S;
			std::cout << tileStr("--", partialVCSize) << "- " << dye("peeking", 'y') << " stack: {";
			while (!SP.empty())
			{
				auto ccurrent = SP.top();
				SP.pop();
				std::cout << "{";
				if (ccurrent.first.size() > 0) std::cout << ccurrent.first.at(0);
				for (int i=1; i< (int) ccurrent.first.size(); i++)
				{
					std::cout << ", " << ccurrent.first.at(i);
				}
				std::cout << "}";
				if(!SP.empty()) std::cout << ", ";
			}
			std::cout << "} of size: " << S.size() << "\n";	

			std::cout << tileStr("--", partialVCSize) << "- " << dye("Deleting", 'p') << " vertices: {";
			if (current.first.size() > 0) std::cout << current.first.at(0);
			for (int i=1; i< (int) current.first.size(); i++)
			{
				std::cout << ", " << current.first.at(i);
			}
			std::cout << "}\n";
			std::cout << tileStr("--", partialVCSize) << "- " << dye("checking", 'y') << " partial solution:\n";
			std::cout << tileStr("--", partialVCSize) << "- " << "best=" << k << ", partialVCSize: " << partialVCSize << "\n";
		}

		// get maxDegVert of remaining active vertices (the way the algorithm branches is determined by the choice of this vertex)
		branchVertex = G->getMaxDegreeVertex();
		branchVertexDegree = G->getVertexDegree(branchVertex);

		if(debug)
		{
			std::cout << tileStr("--", partialVCSize) << "> " << dye("selected", 'b') << " branchVertex: " << branchVertex << " with degree " << branchVertexDegree << "\n";
		}

		// if no active vertex left in graph or no vertex with degree >= 0: (We found a solution)
		if (branchVertex == -1 || branchVertexDegree == 0)
		{
			// if current solution is actually better than the current best solution: update k & vc
			if (k > partialVCSize || (vc == nullptr && k == partialVCSize)) { // TODO: remove condition later when culling earlier with BnB
				vc = G->getInactiveVertices();
				k = partialVCSize;

				if(debug)
				{
					std::cout << tileStr("--", partialVCSize) << "> " << dye("found", 'g') << " VC: {";
					if (vc->size() > 0) std::cout << vc->at(0);
					for (int i=1; i < (int) vc->size(); i++)
					{
						std::cout << ", " << vc->at(i);
					}
					std::cout << "} of size: " << partialVCSize << "\n";
				}
			}
			if(debug)
			{
				std::cout << tileStr("--", partialVCSize) << "- " << dye("Restoring", 'p') << " vertices: {";
				if (current.first.size() > 0) std::cout << current.first.at(0);
				for (int i=1; i < (int) current.first.size(); i++)
				{
					std::cout << ", " << current.first.at(i);
				}
				std::cout << "}\n";
			}
			// traverse back up the search tree
			G->setActive(&current.first);
			partialVCSize -= current.first.size();
			S.pop();
			continue;
		}

		// if maximum search depth k is reached
		// or if the current solution has already been expanded, revert its deletion of vertices and traverse back up to the next partial solution
		if (k <= partialVCSize)
		{
			if(debug)
			{
				std::cout << tileStr("--", partialVCSize) << "> " << dye("reached", 'c') << " search tree depth k=" << k << "\n";
				std::cout << tileStr("--", partialVCSize) << "- " << dye("Restoring", 'p') << " vertices: {";
				if (current.first.size() > 0) std::cout << current.first.at(0);
				for (int i=1; i < (int) current.first.size(); i++)
				{
					std::cout << ", " << current.first.at(i);
				}
				std::cout << "}\n";
			}

			G->setActive(&current.first);
			partialVCSize -= current.first.size();
			S.pop();
			continue;
		}

		// mark this partial solution as expanded and expand subsequently
		current.second = true;
		S.pop();
		S.push(current);

		// solve graph with maxVertDegree <= 2 in linear time
		if (branchVertexDegree <= 2 && false) // TODO: rm false
		{

		}

		// refined search tree branching for maxVertDegree >= 3
		else if (branchVertexDegree >= 1 /*3*/) // TODO: readd 3
		{
			// if k and current partial VC size permit adding the neighbours
			if (k - partialVCSize >= branchVertexDegree)
			{
				if(debug)
				{
					std::cout << tileStr("--", partialVCSize) << "> " << dye("pushing", 'p') << " maxDegreeVertex neighbourhood deletion: {";
					if (G->getNeighbours(branchVertex)->size() > 0) std::cout << G->getNeighbours(branchVertex)->at(0);
					for (int i=1; i<(int) G->getNeighbours(branchVertex)->size(); i++)
					{
						std::cout << ", " << G->getNeighbours(branchVertex)->at(i);
					}
					std::cout << "}\n";
				}

				// add neighbours of the current vertex to the child differential for evaluating the partial vertex cover where all of the branchVertex's neighbours where taken into the vertex cover
				// TODO: is there no way to create this inline?
				std::pair<std::vector<int>, bool> childDifferentialN({*(G->getNeighbours(branchVertex)), false});
				S.push(childDifferentialN);
			}
			// if k and current partial VC size permit adding the current vertex
			if (k - partialVCSize >= 1)
			{
				if(debug)
				{
					std::cout << tileStr("--", partialVCSize) << "> " << dye("pushing", 'p') << " maxDegreeVertex deletion: {" << branchVertex << "}\n";
				}

				// TODO: is there no way to create this inline?
				std::vector<int> bv = std::vector<int>();
				bv.push_back(branchVertex);
				std::pair<std::vector<int>, bool> childDifferentialV({bv, false});
				S.push(childDifferentialV);
			}

		}
	}

	return vc;
}


vector<int>* vertexBranchingSolverIterative(ArrayGraph* G)
{
	int k = G->getLowerBoundVC();
	int u = (int) INFINITY;
	vector<int>* vc = nullptr;

	while (true)
	{
		if(k > u) {
			//std::cout << "Did not find solution within upper bound u=" << u << "\n";
			return vc;
		}
		vc = VCVertexBranchingIterative(G, k, vc, true);
		//if(vc == nullptr) { std::cout << "Did not find solution for k=" << k << "\n\n"; }
		if (vc != nullptr)
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

void writeSolutionToConsole(vector<string>* vc)
{
	for (auto it = vc->begin(); it != vc->end(); ++it)
	{
		cout << *it << endl;
	}
}

/** Execute specific version of program with optional arguments for more prints
 * version
 * 0: ArrayGraph, iterative
 * 1: ArrayGraph, recursive
 * 2: Exercise 1 Solution TODO:
 * ....
 *  with Cycle Bound
 *  with Clique Bound
*/
void chooseImplementationAndOutput(int version = 0, bool printGraph = false, bool printMappings = false, bool printDebug = false, bool showVCSize = false)
{
    std::vector<int>* vc;
    switch(version)
    {
        case 0:
            ArrayGraph* G = ArrayGraph::readStandardInput();
            if (G == nullptr)
                throw invalid_argument("Error constructing graph from input file.");
            if (printGraph)
                G->print();

            vc = vertexBranchingSolverIterative(G, printDebug);
            writeSolutionToConsole(G->getStringsFromVertexIndices(vc));

            if (printMappings)
                G->printMappings(vc);
            if (showVCSize)
                cout << vc->size() << endl;
            break;
        case 1:
            ArrayGraph* G = ArrayGraph::readStandardInput();
            if (G == nullptr)
                throw invalid_argument("Error constructing graph from input file.");
            if (printGraph)
                G->print();

            vc = vertexBranchingSolverRecursive(G);
            writeSolutionToConsole(G->getStringsFromVertexIndices(vc));

            if (printMappings)
                G->printMappings(vc);
            if (showVCSize)
                cout << vc->size() << endl;
            break;
        case 2:
            break;
        case 3:
            break;
        case 4:
            break;
        default:
    }

}

/*----------------------------------------------------------*/
/*-----------------------   Main   -------------------------*/
/*----------------------------------------------------------*/

int main(int argc, char* argv[]) {

	try
	{
        chooseImplementationAndOutput(0);
	}
	catch (const exception& e)
	{
		cerr << "Error while running vertex cover solver.\n";
        cerr << e.what();
	}
}