#include <iostream>
#include <unordered_map> //O(1) for insert and access instead of O(log n) for ordered maps
#include <fstream>  //ifstream file opening
#include <stack>          // std::stack
#include <math.h>          // INFINITY
#include "Graph.h"
#include "ArrayGraph.h"
#include "ColorPrint.h"

using namespace std;

enum VCDebugMode {
	ExecutionTranscript,
	StackAccessTranscript,
	IterationsTaken,
	NoDebug
};

/*----------------------------------------------------------*/
/*---------------   Exercise 2 Solver Code   ---------------*/
/*----------------------------------------------------------*/

vector<int>* vcVertexBranchingRecursive(ArrayGraph* G, int k, int* numRec)
{
    (*numRec)++;
	if (k < 0)
    {
		return nullptr;
	}

	int vertex = G->getMaxDegreeVertex();
    
    //no vertices left
    if (vertex == -1)
    {
        return new vector<int>();
    }

    int vertexDeg = G->getVertexDegree(vertex); 
	//graph has no edges left
	if (vertexDeg == 0)
	{
		return new vector<int>();
	}

	//delete first vertex from graph and explore solution
    G->setInactive(vertex);
	vector<int>* S = vcVertexBranchingRecursive(G, k - 1, numRec);
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
	S = vcVertexBranchingRecursive(G, k - neighbours->size(), numRec);
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
		G->setActive(neighbours);
	}
	return nullptr;
}

vector<int>* vertexBranchingSolverRecursive(ArrayGraph* G, int* numRec)
{
	int k = G->getLowerBoundVC();
	vector<int> *vc;

	while (true)
	{
		vc = vcVertexBranchingRecursive(G, k, numRec);
		if (vc != nullptr)
		{
			return vc;
		}
		k++;
	}
}

/**
 * Iterative Branching with maxDegree heuristic
*/
vector<int>* VCVertexBranchingIterative(ArrayGraph* G, int k, std::vector<int>* vc, bool useDegLEQ2Alg, int* numRec)
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

	while (!S.empty())
	{
		(*numRec)++;
		// retrieve the differential of the current partial vertex cover solution to its parent solution (+ expanded tag)
		current = S.top();

		if (current.second)
		{
			G->setActive(&current.first);
			partialVCSize -= current.first.size();
			S.pop();
			continue;
		}

		// otherwise delete the vertices (defined by the differential calculated in the expansion of the parent search tree node)
		G->setInactive(&current.first);
		partialVCSize += current.first.size();

		// get maxDegVert of remaining active vertices (the way the algorithm branches is determined by the choice of this vertex)
		branchVertex = G->getMaxDegreeVertex();
		branchVertexDegree = G->getVertexDegree(branchVertex);

		// if no active vertex left in graph or no vertex with degree >= 0: (We found a solution)
		if (branchVertex == -1 || branchVertexDegree == 0)
		{
			// if current solution is actually better than the current best solution: update k & vc
			vc = G->getInactiveVertices();
			k = partialVCSize;
			// return vertex cover
			return vc;
		}

		// if maximum search depth k is reached
		// or if the current solution has already been expanded, revert its deletion of vertices and traverse back up to the next partial solution
		if (k <= partialVCSize)
		{
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
		if (useDegLEQ2Alg && branchVertexDegree <= 2)
		{
			// determine partial 2 VC for each connected component
			std::vector<int> deleted;
			std::vector<int> A;
			std::vector<int> CoA;
			bool maxDepthReached = false;
			while(true)
			{
				int origin = G->getConnectedVertex();
				
				// if vertex cover is found, return it
				if(origin == -1 && k >= partialVCSize) {
					vc = G->getInactiveVertices();
					k = partialVCSize;
					return vc;
				}

				if(k <= partialVCSize)
				{
					G->setActive(&deleted);
					partialVCSize -= deleted.size();
					G->setActive(&(current.first));
					partialVCSize -= current.first.size();
					S.pop();
					maxDepthReached = true;
					break;
				}

				auto neighbours = G->getNeighbours(origin);
				int current1 = -1;
				int current2 = -1;

				if(G->getVertexDegree(origin) < 2)
				{
					current1 = neighbours->at(0);
				}
				else
				{
					current1 = neighbours->at(0);
					current2 = neighbours->at(1);
				}
				int lastCurrent1 = origin;
				int lastCurrent2 = origin;
				bool addToA = false;
				A.push_back(origin);

				while(current1 != -1 || current2 != -1)
				{
					(*numRec)++;
					if(current1 == current2)
					{
						if(addToA)
						{
							A.push_back(current1);
						}
						else
						{
							CoA.push_back(current1);
						}
						break;
					}

					if(current1 != -1)
					{
						if(addToA)
						{
							A.push_back(current1);
						}
						else
						{
							CoA.push_back(current1);
						}
						
						int nextCurrent1 = -1;
						for(auto neighbour : *(G->getNeighbours(current1)))
						{
							if(lastCurrent1 != neighbour) { nextCurrent1 = neighbour; }
						}
						lastCurrent1 = current1;
						current1 = nextCurrent1;
					}

					if(current2 != -1)
					{
						if(addToA)
						{
							A.push_back(current2);
						}
						else
						{
							CoA.push_back(current2);
						}
						
						int nextCurrent2 = -1;
						for(auto neighbour : *(G->getNeighbours(current2)))
						{
							if(lastCurrent2 != neighbour) { nextCurrent2 = neighbour; }
						}
						lastCurrent2 = current2;
						current2 = nextCurrent2;
					}
					addToA = !addToA;
				}
				if(maxDepthReached) break;

				// choose smaller partial vc, delete the selected vertices from the graph and save their indices in deleted
				if(A.size() <= CoA.size())
				{
					deleted.insert(
						deleted.end(),
						std::make_move_iterator(A.begin()),
						std::make_move_iterator(A.end())
					);
					G->setInactive(&A);
					partialVCSize += A.size();
				}
				else
				{
					deleted.insert(
						deleted.end(),
						std::make_move_iterator(CoA.begin()),
						std::make_move_iterator(CoA.end())
					);
					G->setInactive(&CoA);
					partialVCSize += CoA.size();
				}
				A.clear();
				CoA.clear();
			}
		}

		// refined search tree branching for maxVertDegree >= 3
		else if ((!useDegLEQ2Alg && branchVertexDegree >= 1) || (useDegLEQ2Alg && branchVertexDegree >= 3))
		{
			// if k and current partial VC size permit adding the neighbours
			if (k - partialVCSize >= branchVertexDegree)
			{
				// add neighbours of the current vertex to the child differential for evaluating the partial vertex cover where all of the branchVertex's neighbours where taken into the vertex cover
				// TODO: is there no way to create this inline?
				std::pair<std::vector<int>, bool> childDifferentialN({*(G->getNeighbours(branchVertex)), false});
				S.push(childDifferentialN);
			}
			// if k and current partial VC size permit adding the current vertex
			if (k - partialVCSize >= 1)
			{
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

string tileStr(string toTile, int n) {
	string tiling = "";
	for (int i=0; i<n; i++)
	{
		tiling += toTile;
	}
	return tiling;
}
/**
 * Iterative Branching with lots of Debug prints
 * Don't try to understand the code from this method
*/
vector<int>* VCVertexBranchingIterativeDebug(ArrayGraph* G, int k, std::vector<int>* vc, bool useDegLEQ2Alg, int* numRec, VCDebugMode debug)
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
	if(debug == ExecutionTranscript)
	{
		std::cout << dye("started", 'g') << " BnB with stack size: " << S.size() << " and k=" << k << "\n";
	}

	while (!S.empty())
	{
		(*numRec)++;
		// retrieve the differential of the current partial vertex cover solution to its parent solution (+ expanded tag)
		current = S.top();

		if (current.second)
		{
			if(debug == ExecutionTranscript)
			{
				std::cout << tileStr("-", partialVCSize) << "- " << dye("Traversing", 'y') << " back up a previously expanded search tree node\n";
				std::cout << tileStr("-", partialVCSize) << "- " << dye("Restoring", 'p') << " vertices: {";
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

		if(debug == ExecutionTranscript)
		{
			auto SP = S;
			std::cout << tileStr("-", partialVCSize) << "- " << dye("peeking", 'y') << " stack: {";
			/* while (!SP.empty())
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
			} */
			std::cout << "} of size: " << S.size() << "\n";

			std::cout << tileStr("-", partialVCSize) << "- " << dye("Deleting", 'p') << " vertices: {";
			if (current.first.size() > 0) std::cout << current.first.at(0);
			for (int i=1; i< (int) current.first.size(); i++)
			{
				std::cout << ", " << current.first.at(i);
			}
			std::cout << "}\n";
			std::cout << tileStr("-", partialVCSize) << "- " << dye("checking", 'y') << " partial solution:\n";
			std::cout << tileStr("-", partialVCSize) << "- " << "k=" << k << ", partialVCSize: " << partialVCSize << "\n";
		}

		if(debug == StackAccessTranscript)
		{
			std::cout << tileStr("-", partialVCSize) << dye("{", 'r');
			if (current.first.size() > 0) std::cout << current.first.at(0);
			for (int i=1; i< (int) current.first.size(); i++)
			{
				std::cout << dye(", ", 'r') << current.first.at(i);
			}
			std::cout << dye("}", 'r') << "\n";
		}

		// get maxDegVert of remaining active vertices (the way the algorithm branches is determined by the choice of this vertex)
		branchVertex = G->getMaxDegreeVertex();
		branchVertexDegree = G->getVertexDegree(branchVertex);

		if(debug == ExecutionTranscript)
		{
			std::cout << tileStr("-", partialVCSize) << "> " << dye("selected", 'b') << " branchVertex: " << branchVertex << " with degree " << branchVertexDegree << "\n";
		}

		// if no active vertex left in graph or no vertex with degree >= 0: (We found a solution)
		if (branchVertex == -1 || branchVertexDegree == 0)
		{
			// if current solution is actually better than the current best solution: update k & vc
			
			vc = G->getInactiveVertices();
			k = partialVCSize;

			if(debug == ExecutionTranscript)
			{
				std::cout << tileStr("-", partialVCSize) << "> " << dye("found", 'g') << " VC: {";
				if (vc->size() > 0) std::cout << vc->at(0);
				for (int i=1; i < (int) vc->size(); i++)
				{
					std::cout << ", " << vc->at(i);
				}
				std::cout << "} of size: " << partialVCSize << "\n";
			}
			if(debug == ExecutionTranscript)
			{
				std::cout << tileStr("-", partialVCSize) << "- " << dye("Restoring", 'p') << " vertices: {";
				if (current.first.size() > 0) std::cout << current.first.at(0);
				for (int i=1; i < (int) current.first.size(); i++)
				{
					std::cout << ", " << current.first.at(i);
				}
				std::cout << "}\n";
			}
			
			if(debug == IterationsTaken) std::cout << "(k=" << k << ") solved in " << (*numRec) << " iterations.\n";
			// return vertex cover
			return vc;
		}

		// if maximum search depth k is reached
		// or if the current solution has already been expanded, revert its deletion of vertices and traverse back up to the next partial solution
		if (k <= partialVCSize)
		{
			if(debug == ExecutionTranscript)
			{
				std::cout << tileStr("-", partialVCSize) << "> " << dye("reached", 'c') << " search tree depth k=" << k << "\n";
				std::cout << tileStr("-", partialVCSize) << "- " << dye("Restoring", 'p') << " vertices: {";
				if (current.first.size() > 0) std::cout << current.first.at(0);
				for (int i=1; i < (int) current.first.size(); i++)
				{
					std::cout << ", " << current.first.at(i);
				}
				std::cout << "}\n";
			}
			if(debug == StackAccessTranscript) {
				std::cout << tileStr("-", partialVCSize) << dye("{", 'g');
				if (current.first.size() > 0) std::cout << current.first.at(0);
				for (int i=1; i< (int) current.first.size(); i++)
				{
					std::cout << dye(", ", 'g') << current.first.at(i);
				}
				std::cout << dye("}", 'g') << "\n";
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
		if (useDegLEQ2Alg && branchVertexDegree <= 2)
		{
			if(debug == ExecutionTranscript)
			{
				std::cout << tileStr("-", partialVCSize) << "> " << dye("started", 'c') << " linear time algorithm on branchVertex: " << branchVertex << "\n";
			}
			// determine partial 2 VC for each connected component
			std::vector<int> deleted;
			std::vector<int> A;
			std::vector<int> CoA;
			bool maxDepthReached = false;
			while(true)
			{
				int origin = G->getConnectedVertex();
				if(debug == ExecutionTranscript)
				{
					std::cout << tileStr("-", partialVCSize) << "> " << dye("selected", 'b') << " propagation origin: " << origin << "\n";
				}
				// if vertex cover is found, return it
				if(origin == -1 && k >= partialVCSize) {
					vc = G->getInactiveVertices();
					k = partialVCSize;
					if(debug == ExecutionTranscript)
					{
						std::cout << tileStr("-", partialVCSize) << "> " << dye("found", 'g') << " VC: {";
						if (vc->size() > 0) std::cout << vc->at(0);
						for (int i=1; i < (int) vc->size(); i++)
						{
							std::cout << ", " << vc->at(i);
						}
						std::cout << "} of size: " << partialVCSize << "\n";
					}
					return vc;
				}

				if(k <= partialVCSize)
				{
					/* if(debug == NoDebug)
					{
						auto SP = S;
						std::cout << tileStr("-", partialVCSize) << "- " << dye("peeking", 'y') << " stack: {";
						for (int j=0; !SP.empty() && j<7; j++)
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
	
						std::cout << tileStr("-", partialVCSize) << "> " << dye("reached", 'c') << " search tree depth k=" << k << "\n";
						std::cout << tileStr("-", partialVCSize) << "- " << dye("Restoring", 'p') << " vertices: {";
						if (current.first.size() > 0) std::cout << current.first.at(0);
						for (int i=1; i < (int) current.first.size(); i++)
						{
							std::cout << ", " << current.first.at(i);
						}
						std::cout << "}\n";
					} */
					//std::cout << tileStr("-", partialVCSize) << "> " << dye("reached", 'c') << " search tree depth k=" << k << " with partial VC size: " << partialVCSize << " and linVC: " << std::min(A.size(), CoA.size()) << "\n";
					G->setActive(&deleted);
					partialVCSize -= deleted.size();
					G->setActive(&(current.first));
					partialVCSize -= current.first.size();
					S.pop();
					maxDepthReached = true;
					break;
				}

				auto neighbours = G->getNeighbours(origin);
				int current1 = -1;
				int current2 = -1;

				if(G->getVertexDegree(origin) < 2)
				{
					current1 = neighbours->at(0);
				}
				else
				{
					current1 = neighbours->at(0);
					current2 = neighbours->at(1);
				}
				int lastCurrent1 = origin;
				int lastCurrent2 = origin;
				bool addToA = false;
				A.push_back(origin);

				while(current1 != -1 || current2 != -1)
				{
					(*numRec)++;
					if(current1 == current2)
					{
						if(addToA)
						{
							A.push_back(current1);
						}
						else
						{
							CoA.push_back(current1);
						}
						break;
					}

					if(current1 != -1)
					{
						if(addToA)
						{
							A.push_back(current1);
						}
						else
						{
							CoA.push_back(current1);
						}
						
						int nextCurrent1 = -1;
						for(auto neighbour : *(G->getNeighbours(current1)))
						{
							if(lastCurrent1 != neighbour) { nextCurrent1 = neighbour; }
						}
						lastCurrent1 = current1;
						current1 = nextCurrent1;
					}

					if(current2 != -1)
					{
						if(addToA)
						{
							A.push_back(current2);
						}
						else
						{
							CoA.push_back(current2);
						}
						
						int nextCurrent2 = -1;
						for(auto neighbour : *(G->getNeighbours(current2)))
						{
							if(lastCurrent2 != neighbour) { nextCurrent2 = neighbour; }
						}
						lastCurrent2 = current2;
						current2 = nextCurrent2;
					}
					addToA = !addToA;
				}
				if(maxDepthReached) break;

				// choose smaller partial vc, delete the selected vertices from the graph and save their indices in deleted
				if(A.size() <= CoA.size())
				{
					deleted.insert(
						deleted.end(),
						std::make_move_iterator(A.begin()),
						std::make_move_iterator(A.end())
					);
					G->setInactive(&A);
					partialVCSize += A.size();
					//std::cout << tileStr("-", partialVCSize) << "> " << dye("Deleting", 'y') << " A with size " << A.size() << " ===> (" << partialVCSize << "/" << k << ")\n";
				}
				else
				{
					deleted.insert(
						deleted.end(),
						std::make_move_iterator(CoA.begin()),
						std::make_move_iterator(CoA.end())
					);
					G->setInactive(&CoA);
					partialVCSize += CoA.size();
					//std::cout << tileStr("-", partialVCSize) << "> " << dye("Deleting", 'y') << " CoA with size " << CoA.size() << " ===> (" << partialVCSize << "/" << k << ")\n";
				}
				//std::cout << tileStr("-", partialVCSize) << "> " << dye("A + CoA = ", 'y') << A.size() + CoA.size() << " ===> (" << partialVCSize << "/" << k << ")" << "\n";
				if(debug == ExecutionTranscript)
				{
					std::cout << tileStr("-", partialVCSize) << "- " << dye("A =", 'y') << " {";
					if (A.size() > 0) std::cout << A.at(0);
					for (int i=1; i < (int) A.size(); i++)
					{
						std::cout << ", " << A.at(i);
					}
					std::cout << "}\n";

					std::cout << tileStr("-", partialVCSize) << "- " << dye("CoA =", 'y') << " {";
					if (CoA.size() > 0) std::cout << CoA.at(0);
					for (int i=1; i < (int) CoA.size(); i++)
					{
						std::cout << ", " << CoA.at(i);
					}
					std::cout << "}\n";
				}
				A.clear();
				CoA.clear();
			}
		}

		// refined search tree branching for maxVertDegree >= 3
		else if ((!useDegLEQ2Alg && branchVertexDegree >= 1) || (useDegLEQ2Alg && branchVertexDegree >= 3))
		{
			// if k and current partial VC size permit adding the neighbours
			if (k - partialVCSize >= branchVertexDegree)
			{
				if(debug == ExecutionTranscript)
				{
					std::cout << tileStr("-", partialVCSize) << "> " << dye("pushing", 'p') << " maxDegreeVertex neighbourhood deletion: {";
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
				if(debug == ExecutionTranscript)
				{
					std::cout << tileStr("-", partialVCSize) << "> " << dye("pushing", 'p') << " maxDegreeVertex deletion: {" << branchVertex << "}\n";
				}

				// TODO: is there no way to create this inline?
				std::vector<int> bv = std::vector<int>();
				bv.push_back(branchVertex);
				std::pair<std::vector<int>, bool> childDifferentialV({bv, false});
				S.push(childDifferentialV);
			}

		}
	}

	if(debug == IterationsTaken && vc == nullptr) { std::cout << "(k=" << k << ") determined to have no solution in " << (*numRec) << " iterations.\n"; }
	return vc;
}

vector<int>* vertexBranchingSolverIterative(ArrayGraph* G, bool useDegLEQ2Alg, int* numRec,  VCDebugMode debug)
{
	int k = G->getLowerBoundVC();
	vector<int>* vc = nullptr;

	if(G->getVertexCount() == 0) return new vector<int>();
	while (true)
	{
        if(debug != NoDebug)
        {
            vc = VCVertexBranchingIterative(G, k, vc, useDegLEQ2Alg, numRec);
        }
        else
        {
            vc = VCVertexBranchingIterativeDebug(G, k, vc, useDegLEQ2Alg, numRec, debug);
        }
		
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

vector<int>* vcBranchArrayGraph(ArrayGraph* G, int k)
{
    
    if (k < 0)
        return nullptr;

    pair<int, int>* edge = G->getFirstValidEdge();

    //graph has no edges left
    if (edge == nullptr)
    {
        return new vector<int>();
    }

    //delete first vertex from graph and explore solution
    G->setInactive(edge->first);
    vector<int>* S = vcBranchArrayGraph(G, k - 1);
    if (S != nullptr)
    {
        //revert changes to graph
        G->setActive(edge->first);
        //return results
        S->push_back(edge->first);
        return S;
    }
    else
    {
        //revert changes to graph
        G->setActive(edge->first);
    }


    //delete second vertex from graph and explore solution
    G->setInactive(edge->second);
    S = vcBranchArrayGraph(G, k - 1);
    if (S != nullptr)
    {
        //revert changes to graph
        G->setActive(edge->second);
        //return results
        S->push_back(edge->second);
        return S;
    }
    else
    {
        //revert changes to graph
        G->setActive(edge->second);
    }
    return nullptr;
}

vector<int>* searchTreeSolveArrayGraph(ArrayGraph* G)
{
    int k = 0;
    vector<int> *vc;

    while (true)
    {
        vc = vcBranchArrayGraph(G, k);
        if (vc != nullptr)
        {
            return vc;
        }
        k++;
    }
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
 * 2: Arraygraph, Exercise 1 Branching
 * 3: Graph, Exercise 1 Branching with Hashmap graph
 * ....
*/
void chooseImplementationAndOutput(int version = 0, bool printGraph = false, bool printMappings = false, bool printDebug = false, bool showVCSize = false, bool printVC = true, bool printBounds = false)
{

    if(version == 0)
    {
        std::vector<int>* vc;
        ArrayGraph* G = ArrayGraph::readStandardInput();

		if (G == nullptr)
			throw invalid_argument("Error constructing graph from input file.");
		if (printGraph)
			G->print();

        pair<int, int> lowerBounds;
        if(printBounds)
        {
            lowerBounds = G->getAllLowerBounds();
        }

		int numRecursiveSteps = 0;
		vc = vertexBranchingSolverIterative(G, true, &numRecursiveSteps,/* ExecutionTranscript */NoDebug);
		
		if(printVC)
		{
			writeSolutionToConsole(G->getStringsFromVertexIndices(vc));
			cout << "#recursive steps: " << numRecursiveSteps << endl;
		}
		if (printMappings)
			G->printMappings(vc);
		if (showVCSize)
			cout << "VC size: " << vc->size() << endl;

        if(printBounds)
        {
            //cout << "Clique bound: " << lowerBounds.first << ", Cycle bound: " << lowerBounds.second;
            //cout << "last-k: " << lowerBounds.first << "-" << lowerBounds.second;
            if(lowerBounds.first < lowerBounds.second) 
            {
                cout << "#recursive steps: " << lowerBounds.second;
            }
            else
            {
                cout << "#recursive steps: " << lowerBounds.first;
            }

        }
    }
    else if(version == 1)
    {
        std::vector<int>* vc;
        ArrayGraph* G = ArrayGraph::readStandardInput();
        if (G == nullptr)
            throw invalid_argument("Error constructing graph from input file.");
        if (printGraph)
            G->print();

        int numRecursiveSteps = 0;
        vc = vertexBranchingSolverRecursive(G, &numRecursiveSteps);
		if(printVC)
        {
            writeSolutionToConsole(G->getStringsFromVertexIndices(vc));
            cout << "#recursive steps: " << numRecursiveSteps << endl;
        }


        if (printMappings)
            G->printMappings(vc);
        if (showVCSize)
            cout << "VC size: " << vc->size() << endl;
    }
    else if(version == 2)
    {
        std::vector<int>* vc;
        ArrayGraph* G = ArrayGraph::readStandardInput();
        if (G == nullptr)
            throw invalid_argument("Error constructing graph from input file.");
        if (printGraph)
            G->print();

        vc = searchTreeSolveArrayGraph(G);
		if(printVC)
        	writeSolutionToConsole(G->getStringsFromVertexIndices(vc));
        if (printMappings)
			G->printMappings(vc);
		if (showVCSize)
			cout << "VC size: " << vc->size() << endl;
    }
    else if(version == 3)
    {
        std::vector<string>* vc;
        Graph* G = Graph::readStandardInput();
        if (G == nullptr)
            throw invalid_argument("Error constructing graph from input file.");
        if (printGraph)
            G->print();

        vc = searchTreeSolve(G);
		if(printVC)
        	writeSolutionToConsole(vc);
		if (showVCSize)
			cout << "VC size: " << vc->size() << endl;
    }
}

/*----------------------------------------------------------*/
/*-----------------------   Main   -------------------------*/
/*----------------------------------------------------------*/

//TODO: recursive steps print für Iterative anpassen
int main(int argc, char* argv[]) {

	try
	{
        //chooseImplementationAndOutput(0, false, false, false, false, true, true);
        chooseImplementationAndOutput(0);
	}
	catch (const exception& e)
	{
		cerr << "Error while running vertex cover solver.\n";
        cerr << e.what();
	}
}