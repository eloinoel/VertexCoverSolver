// VertexCoverSolver.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <unordered_map> //O(1) for insert and access instead of O(log n) for ordered maps
#include <fstream>  //ifstream file opening
#include "Graph.h"

using namespace std;


void deleteStringFromVector(vector<string>* vec, string element)
{
	auto position = std::find(vec->begin(), vec->end(), element);
	if (position != vec->end()) // == myVector.end() means the element was not found
		vec->erase(position);
	return;
}

vector<string>* vcBranch(Graph* G, Graph* graphCopy, int k, vector<string>* deletedVertices)
{
	
	//cout << "deleted vertices: " << Graph::getVectorContentsString(deletedVertices) << "\n"; //TODO: delete debug
	if (k < 0)
		return NULL;

	pair<string, string>* edge = graphCopy->getFirstValidEdge();

	//graph has no edges left
	if (edge == nullptr)
	{
		//cout << "found vc " << Graph::getVectorContentsString(deletedVertices) << "\n"; //TODO: delete debug
		return new vector<string>();
	}
			

	//TODO: delete debug
	/*cout << "Before stuffing around: size: " << graphCopy->getSize();
	graphCopy->print();*/
	//cout << "next edge: " << edge->first << ", " << edge->second << "\n";

	//delete first vertex from graph and explore solution
	graphCopy->deleteAdjacencyMapEntry(edge->first);
	deletedVertices->push_back(edge->first);
	vector<string>* S = vcBranch(G, graphCopy, k - 1, deletedVertices);
	if (S != NULL)
	{
		//revert changes to graph
		graphCopy->addAdjacencyMapEntry(G, edge->first); //revert changes to graph
		deleteStringFromVector(deletedVertices, edge->first);
		//return results
		S->push_back(edge->first);
		return S;
	}
	else
	{
		//revert changes to graph
		graphCopy->addAdjacencyMapEntry(G, edge->first);
		deleteStringFromVector(deletedVertices, edge->first);
	}


	//delete second vertex from graph and explore solution
	graphCopy->deleteAdjacencyMapEntry(edge->second);
	deletedVertices->push_back(edge->second);
	S = vcBranch(G, graphCopy, k - 1, deletedVertices);
	if (S != NULL)
	{
		//revert changes to graph
		graphCopy->addAdjacencyMapEntry(G, edge->second);
		deleteStringFromVector(deletedVertices, edge->second);
		//return results
		S->push_back(edge->second);
		return S;
	}
	else
	{
		//revert changes to graph
		graphCopy->addAdjacencyMapEntry(G, edge->second); 
		deleteStringFromVector(deletedVertices, edge->second);
	}


	//TODO: delete debug
	//cout << "After stuffing around: size: " << graphCopy->getSize();
	//graphCopy->print();
	return NULL;
	
	
}


vector<string>* searchTreeSolve(Graph* G)
{
	int k = 0;
	vector<string> *vc;
	Graph* graphCopy = G->copy();
	vector<string> *deletedVertices = new vector<string>();

	while (true)
	{
		//cout << "-----------------------------------\nk = " << k << "\n"; //TODO: delete debug
		vc = vcBranch(G, graphCopy, k, deletedVertices);
		if (vc != NULL)
		{
			return vc;
		}
		k++;
		if (k == 4) {
			return vc;
		}
	}
}

void writeSolutionToFile(string fileName, vector<string>* vc)
{
	ofstream outfile(fileName);
	for (auto it = vc->begin(); it != vc->end(); ++it)
	{
		outfile << *it << endl;
	}
	outfile.close();
}

void printUsage()
{
	cerr << "Usage: ./VertexCoverSolver.exe <file path>";
}

int main(int argc, char* argv[]) {
	if (argc > 3)
	{
		printUsage();
	}

	string fileName = argv[1];
	string outputFileName;
	if (argc == 3)
		outputFileName = argv[2];
	else
		outputFileName = "prog_out.txt";


	try
	{
		Graph* G = Graph::readInput(fileName);
		if (G == NULL)
		{
			cerr << "Error constructing graph from input file.";
		}

		//test vc solver
		vector<string>* vc = searchTreeSolve(G);
		//cout << "Vertex Cover with k = " << vc->size() << ": " << Graph::getVectorContentsString(vc);
		writeSolutionToFile(outputFileName, vc);
	}
	catch (const exception&)
	{
		cerr << "Check if provided path to file is valid.";
	}
	
}

