#ifndef GRAPH_H
#define GRAPH_H

#include <string>
#include <unordered_map> //O(1) for insert and access instead of O(log n) for ordered maps
#include <vector>
#include <algorithm>

using namespace std;

class Graph
{
public:
	
private:
	unordered_map<string, vector<string>> adjacencyMap;
	//unordered_map<string, unordered_set<pair<string, string>>> adjacencyMap; // Map: Vertex --> { Edge(v1, v2) }

public:
	void addVertex(string s);

	void addEdge(string first, string second);

	vector<string>* deleteVertex(string v);

	void deleteAdjacencyMapEntry(string v);

	void addAdjacencyMapEntry(Graph* graphToCopyFrom, string vertex);

	//int getNumEdges();

	pair<string, string>* getFirstValidEdge();

	int getSize();

	void print();

	static string getVectorContentsString(vector<string> *vec);

	Graph* copy();

	/* generate graph from a text file */
	static Graph* readInputFromFile(string fileName);

	static Graph* readStandardInput();

	static Graph* readStandardInputNoGetLine();

	static vector<string> extractVertices(string line);

private:
	static string eraseLeadingTrailingWhitespacesFromString(string str);

	/* tests whether a char fulfills vertex naming format*/
	static bool isVertexCharacter(char c);
};

#endif
