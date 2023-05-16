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
    std::vector<std::pair<bool, int>>* graphState;

    /* maps from original vertex name from input data to index and degree */
    std::unordered_map<std::string, std::pair<int, int>> originalVertexNames;

    /* TODO: maximum vertex degree, when all vertices are active */
    unsigned int maxDegree;

//functions
private:
    inline ArrayGraph() {};

	static std::string eraseLeadingTrailingWhitespacesFromString(std::string str);

    /* tests whether a char fulfills vertex naming format*/
	static bool isVertexCharacter(char c);

    /* initially sets up activity flags & node degree data structure */
    void initGraphState();

public:

    // TODO: make more efficent, employ less copies
    static ArrayGraph* readStandardInput();

    void print();

    void printOriginalVertexNames();

    /* calculate lower bound for VC */
    int getLowerBoundVC();

    /* get flag array for vertices being active: */
    //std::vector<bool>* getActive();

    /* set vertices with the passed vertex indices to active and update degrees for those vertices and their neighbours: */
    void setActive(int vertexIndex);
    void setActive(std::vector<int>* vertexIndices);

    /* set vertices with the passed vertex indices to inactive and update degrees for their neighbours: */
    void setInactive(int vertexIndex);
    void setInactive(std::vector<int>* vertexIndices);

    /* get the indices of all vertices that are inactive */
    std::vector<int>* getInactiveVertices();

    /* get graph state array: */
    inline std::vector<std::pair<bool, int>>* getState();

    /* get vertex with highest degree among active vertices:
    * returns -1 if no vertex in graph
    */
    int getMaxDegreeVertex();

    /* get degree of the vertex with passed vertex index: */
    inline int getVertexDegree(int vertexIndex);
    
    /* get indices of a vertices neighbours: */
    std::vector<int>* getNeighbours(int vertexIndex);

};

#endif