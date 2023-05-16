#ifndef ARRAYGRAPH_H
#define ARRAYGRAPH_H

#include <string>
#include <unordered_map> //O(1) for insert and access instead of O(log n) for ordered maps
#include <vector>
#include <algorithm>

class ArrayGraph
{
//variables
private:
    /* each index represents a vertex, the content of each index represents the adjacent vertices' indices */
    std::vector<std::vector<int>*> adjacencyList;

    /* contains meta information about the graph's vertices: 
    * 0: True = vertex in graph
    * 1: degree of vertex
    */
    std::vector<std::pair<bool, int>> graphState;

    /* maps from original vertex name from input data to index and degree */
    std::unordered_map<std::string, std::pair<int, int>> originalVertexNames;


//functions
private:
    inline ArrayGraph() {};

	static std::string eraseLeadingTrailingWhitespacesFromString(std::string str);

    /* tests whether a char fulfills vertex naming format*/
	static bool isVertexCharacter(char c);

public:

    static ArrayGraph* readStandardInput();

    void print();

    void printOriginalVertexNames();

    /* calculate lower bound for VC */
    int getLowerBoundVC();

    /* get flag array for vertices being active: */
    std::vector<bool>* getActive();

    /* set flags in flags array exactly for the passed vertex indices and update degrees for those vertices and their neighbours: */
    void setActive(std::vector<int, bool>* vertexFlags);

    /* set all inactive vertices except the passed ones to be active and set all active vertices, that are among the passed ones to be inactive: */
    void setInactiveVertices(std::vector<int>* vertices);

    /* get graph state array: */
    std::vector<std::pair<bool, int>>* getState();

    /* get vertex with highest degree among active vertices:
    * returns -1 if no vertex in graph
    */
    int getMaxDegreeVertex();
    
    /* get indices of a vertices neighbours: */
    std::vector<int>* getNeighbours(int vertexIndex);

};

#endif