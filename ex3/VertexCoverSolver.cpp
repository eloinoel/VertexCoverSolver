#include <iostream>
#include <fstream>  //ifstream file opening
#include <stack>          // std::stack
#include <math.h>          // INFINITY
#include "utils/ArrayGraph.h"
#include "utils/ColorPrint.h"
#include "utils/BucketGraph.h"


#include <signal.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;
typedef ColorPrint cp;

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
    cut = G->reduce(&k);
    //cut = G->reduce(&k);
    if(cut)
    {
        //std::cout << "> cutting through data reductions " << std::endl;
        G->unreduce(&k, previousK);
        return nullptr;
    }

    //std::cout << "> calculated LPBound: " << G->getLPBound() << " with k=" << k << std::endl;
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
    G->setInactive(neighbours);
	S = vcVertexBranchingRecursive(G, k - neighbours->size(), depth+1, numRec);
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
	int k = G->getLowerBoundVC();

//    G->printBucketQueue();
//    G->print();

	unordered_map<int, bool>* vc;

	while (true)
	{
        // Reduction Rules
        //vector<ReductionVertices>* reductionVertices = new vector<ReductionVertices>;

        // Apply Reduction Rules for the first time
       /*  if(G->applyReductionRules(&k, reductionVertices))
            return nullptr; */
		vc = vcVertexBranchingRecursive(G, k, 0, numRec);
		if (vc != nullptr)
		{
            // Add Reduced Vertices to Vertex Cover
            /* G->addReducedVertices(vc, reductionVertices);
            delete reductionVertices; */

			return vc;
		}
        /* G->addBackReducedVertices(&k, reductionVertices);
        delete reductionVertices; */
        k++;
	}
}


/*----------------------------------------------------------*/
/*---------------   Exercise 2 Solver Code   ---------------*/
/*----------------------------------------------------------*/

vector<int>* vcVertexBranchingRecursiveEx2(ArrayGraph* G, int k, int* numRec)
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
    /* vector<ReductionVertices>* reductionVertices = new vector<ReductionVertices>;

    if(!G->applyReductionRules(&k, reductionVertices)) {
        delete reductionVertices;
        return nullptr;
    } */
    //================================================================

	//delete first vertex from graph and explore solution
    G->setInactive(vertex);
	vector<int>* S = vcVertexBranchingRecursiveEx2(G, k - 1, numRec);
	if (S != nullptr)
	{
		//return results

        // Add Reduced Vertices to Vertex Cover
        /* G->addReducedVertices(S, merged, deletedReduced);
        delete reductionVertices; */

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
	S = vcVertexBranchingRecursiveEx2(G, k - neighbours->size(), numRec);
	if (S != nullptr)
	{
		//return results
        for (int i = 0; i < (int) neighbours->size(); i++)
        {
            S->push_back(neighbours->at(i));
        }

        // Add Reduced Vertices to Vertex Cover
        /* G->addReducedVertices(S, merged, deletedReduced);
        delete reductionVertices; */

        return S;
	}
	else
	{
		//revert changes to graph
		G->setActive(neighbours);
	}

    //================================================================
    // Reverse Data Reduction
    /* G->addBackReducedVertices(k, reductionVertices);
    delete reductionVertices; */
    //================================================================

	return nullptr;
}

vector<int>* vertexBranchingSolverRecursiveEx2(ArrayGraph* G, int* numRec)
{
	int k = G->getLowerBoundVC();
    //int k = 0;

	vector<int> *vc;

	while (true)
	{
        // Reduction Rules
        //vector<ReductionVertices>* reductionVertices = new vector<ReductionVertices>;

        // Apply Reduction Rules for the first time
       /*  if(G->applyReductionRules(&k, reductionVertices))
            return nullptr; */

		vc = vcVertexBranchingRecursiveEx2(G, k, numRec);
		if (vc != nullptr)
		{
            // Add Reduced Vertices to Vertex Cover
            /* G->addReducedVertices(vc, reductionVertices);
            delete reductionVertices; */

			return vc;
		}

        /* G->addBackReducedVertices(&k, reductionVertices);
        delete reductionVertices; */

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
/*-------------------   Data Reduction   -------------------*/
/*----------------------------------------------------------*/


BucketGraph* bucketGraph = nullptr;

void my_sig_handler(sig_atomic_t s)
{
    //print original edges if timeout in data reduction for (b) on exercise sheet
    if(bucketGraph != nullptr)
    {
        std::vector<std::string>* str = bucketGraph->getOriginalEdgesToConsoleString();
        for(int i = 0; i < (int) str->size(); i++)
        {
            cout << str->at(i);
        }
        cout << "#difference: " << 0 << endl;
    }
}

/** Execute specific version of program with optional arguments for more prints
 * version
 * 0: ArrayGraph, recursive
 * 1: Bucketgraph
 * 5: (b) apply data reduction and output smaller graph and diff in vc size
 * ....
*/
void chooseImplementationAndOutput(int version = 1, bool printGraph = false, bool printMappings = false, 
bool printDebug = false, bool printVCSize = false, bool printVC = true, bool printBounds = false)
{

    if(version == 0)
    {
        std::vector<int>* vc;
        ArrayGraph* G = ArrayGraph::readStandardInput();
        if (G == nullptr)
            throw invalid_argument("Error constructing graph from input file.");
        if (printGraph)
            G->print();


		if(printVC)
        {
            int numRecursiveSteps = 0;
            vc = vertexBranchingSolverRecursiveEx2(G, &numRecursiveSteps);
            writeSolutionToConsole(G->getStringsFromVertexIndices(vc));
            cout << "#recursive steps: " << numRecursiveSteps << endl;
        }

        if(printBounds)
        {
            int bound = G->getLowerBoundVC();
            cout << "#recursive steps: " << bound << endl;
        }


        if (printMappings)
            G->printMappings(vc);
        if (printVCSize)
            cout << "VC size: " << vc->size() << endl;
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
            writeSolutionToConsole(G->getStringsFromVertexIndices(vc));
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
        unordered_map<int, bool>* vc = vcSolverRecursive(bucketGraph, &numRecursiveSteps);

        int tmpK = vc->size();
        bucketGraph->reduce(&tmpK);
        
        std::vector<std::string>* str = bucketGraph->getEdgesToConsoleString();

        //bucketGraph->resetLPBoundDataStructures();
        //unordered_map<int, bool>* reducedVc = vcSolverRecursive(bucketGraph, &numRecursiveSteps);
        for(int i = 0; i < (int) str->size(); i++)
        {
            cout << str->at(i);
        }
        cout << "#difference: " << to_string((bucketGraph->getOriginalEdgeCount() - str->size())) << endl;
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
        chooseImplementationAndOutput(1, false, false, false, false, true, false);
        //chooseImplementationAndOutput(1, true, false, false, true, true, false); //print alot
        signal(SIGINT, my_sig_handler); //catches SIGINT to output anything
        //chooseImplementationAndOutput(5, false, false, false, false, true, false);
	}
	catch (const exception& e)
	{
		cerr << ColorPrint::dye("Error while running vertex cover solver.\n", 'r');
        cerr << ColorPrint::dye(e.what() ,'r');
	}
}