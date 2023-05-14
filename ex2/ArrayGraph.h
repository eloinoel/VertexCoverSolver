#ifndef ARRAYGRAPH_H
#define ARRAYGRAPH_H

#include <string>
#include <unordered_map> //O(1) for insert and access instead of O(log n) for ordered maps
#include <vector>
#include <algorithm>

using namespace std;

class ArrayGraph
{
//variables
private:
    /* each index represents a vertex, the content of each index represents the adjacent vertices' indices */
    vector<vector<int>*> adjacencyList;

    /* contains meta information about the graph's vertices: 
    * 0: True = vertex in graph
    * 1: degree of vertex
    */
    vector<bool, int> graphState;

    /* maps from original vertex name from input data to index and degree */
    unordered_map<string, pair<int, int>> originalVertexNames;


//functions
private:
    static ArrayGraph* readStandardInput();

	static string eraseLeadingTrailingWhitespacesFromString(string str);

    /* tests whether a char fulfills vertex naming format*/
	static bool isVertexCharacter(char c);

    void addEdge(string first, string second);

public:

    /* calculate lower bound for VC */
    int getLowerBoundVC();

    /* get flag array for vertices being active: */
    vector<bool>* getActive();

    /* set flags in flags array exactly for the passed vertex indices and update degrees for those vertices and their neighbours: */
    void setActive(vector<int, bool>* vertexFlags);

    /* set all inactive vertices except the passed ones to be active and set all active vertices, that are among the passed ones to be inactive: */
    void setInactiveVertices(std::vector<int>* vertices);

    /* get graph state array: */
    vector<bool, int>* getState();

    /* get vertex with highest degree among active vertices:
    * returns -1 if no vertex in graph
    */
    int getMaxDegreeVertex();
    
    /* get indices of a vertices neighbours: */
    vector<int>* getNeighbours(int vertexIndex);

};

#endif