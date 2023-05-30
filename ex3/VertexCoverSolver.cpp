#include <iostream>
#include <unordered_map> //O(1) for insert and access instead of O(log n) for ordered maps
#include <fstream>  //ifstream file opening
#include <stack>          // std::stack
#include <math.h>          // INFINITY
#include "utils/ArrayGraph.h"
#include "utils/ColorPrint.h"

using namespace std;

//enum VCDebugMode {
//	ExecutionTranscript,
//	StackAccessTranscript,
//	IterationsTaken,
//	NoDebug
//};

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

    //================================================================
    // Reduction Rules
    vector<ReductionVertices>* reductionVertices = new vector<ReductionVertices>;

    if(!G->applyReductionRules(&k, reductionVertices)) {
        delete reductionVertices;
        return nullptr;
    }
    //================================================================

	//delete first vertex from graph and explore solution
    G->setInactive(vertex);
	vector<int>* S = vcVertexBranchingRecursive(G, k - 1, numRec);
	if (S != nullptr)
	{
		//return results

        // Add Reduced Vertices to Vertex Cover
        G->addReducedVertices(S, reductionVertices);
        delete reductionVertices;

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
		//return results
        for (int i = 0; i < (int) neighbours->size(); i++)
        {
            S->push_back(neighbours->at(i));
        }

        // Add Reduced Vertices to Vertex Cover
        G->addReducedVertices(S, reductionVertices);
        delete reductionVertices;

        return S;
	}
	else
	{
		//revert changes to graph
		G->setActive(neighbours);
	}

    //================================================================
    // Reverse Data Reduction
    G->addBackReducedVertices(&k, reductionVertices);
    delete reductionVertices;
    //================================================================

	return nullptr;
}

vector<int>* vertexBranchingSolverRecursive(ArrayGraph* G, int* numRec)
{
//	int k = G->getLowerBoundVC();
	int k = 1;
//	int k = 2;
//	int k = 3;
//	int k = 4;

	vector<int> *vc;



	while (true)
	{
        // Reduction Rules
        vector<ReductionVertices>* reductionVertices = new vector<ReductionVertices>;

        // Apply Reduction Rules for the first time
        if(G->applyReductionRules(&k, reductionVertices))
            return nullptr;

		vc = vcVertexBranchingRecursive(G, k, numRec);
		if (vc != nullptr)
		{
            // Add Reduced Vertices to Vertex Cover
            G->addReducedVertices(vc, reductionVertices);
            delete reductionVertices;

			return vc;
		}

        G->addBackReducedVertices(&k, reductionVertices);
        delete reductionVertices;

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
void chooseImplementationAndOutput(int version = 0, bool printGraph = false, bool printMappings = false, bool printDebug = false, bool printVCSize = false, bool printVC = true, bool printBounds = false)
{
    if(version == 0)
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
        if (printVCSize)
            cout << "VC size: " << vc->size() << endl;
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

	try
	{
        chooseImplementationAndOutput(0, false, false, false, false, true, false);
	}
	catch (const exception& e)
	{
		cerr << "Error while running vertex cover solver.\n";
        cerr << e.what();
	}
}