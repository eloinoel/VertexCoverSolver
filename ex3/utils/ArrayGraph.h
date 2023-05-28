#ifndef ARRAYGRAPH_H
#define ARRAYGRAPH_H

#include <string>
#include <unordered_map> //O(1) for insert and access instead of O(log n) for ordered maps
#include <vector>
#include <algorithm>
#include <list>

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

    int numberOfVertices;
    int numberOfEdges;

public:

//functions
private:
protected:
    inline ArrayGraph() {};

	static std::string eraseLeadingTrailingWhitespacesFromString(std::string str);

    /* tests whether a char fulfills vertex naming format*/
	static bool isVertexCharacter(char c);

    /* initially sets up activity flags & node degree data structure */
    void initGraphState(int vertexCount, int edgeCount);

    void printVector(std::list<int>* vec, std::string name);

    bool vertexCanBeAddedToClique(int vertex, std::vector<int>* clique);
    int partition(std::vector<int>* toSort, int low, int high);
    void quickSort(std::vector<int>* toSort, int low, int high);
    bool contains(std::vector<int>* vertexIndices, int vertexIndex);

public:
    bool isVertexCoverFound();

    // TODO: make more efficent, employ less copies
    static ArrayGraph* readStandardInput();

    std::vector<std::vector<int>*>* getAdjacencyList() { return &adjacencyList; };
    void setAdjacencyList(std::vector<std::vector<int>*> _adjacencyList) { adjacencyList = _adjacencyList; };

    inline int getVertexCount() { return adjacencyList.size(); }

    std::vector<int>* getVerticesSortedByDegree();

    void print();

    void printOriginalVertexNames();

    inline int getNumberOfVertices() { return numberOfVertices; }
    inline int getNumberOfEdges() { return numberOfEdges; }

    /* calculate lower bound for VC */
    int getCliqueBound();
    int getLPBound();
    int getLPCycleBound();
    int getLowerBoundVC();

    std::vector<int> getAllLowerBounds();

    /* get flag array for vertices being active: */
    //std::vector<bool>* getActive();

    /* set vertices with the passed vertex indices to active and update degrees for those vertices and their neighbours: */
    void setActive(int vertexIndex);
    void setActive(std::vector<int>* vertexIndices);

    /* set vertices with the passed vertex indices to inactive and update degrees for their neighbours: */
    void setInactive(int vertexIndex);
    void setInactive(std::vector<int>* vertexIndices);

    inline bool isActive(int vertexIndex) { return graphState->at(vertexIndex).first; };

    /* get the indices of all vertices that are inactive */
    std::vector<int>* getInactiveVertices();

    /* get graph state array: */
    inline std::vector<std::pair<bool, int>>* getState() { return graphState; };

    /* get active vertex with degree at least one:
    * returns -1 if no vertex in graph
    */
    int getConnectedVertex();

    /* get active vertex with highest degree among active vertices:
    * returns -1 if no vertex in graph
    */
    int getMaxDegreeVertex();
    /* get max degree vertex out of candidates */
    int getMaxDegreeVertex(std::vector<int>* candidates);

    /* get degree of the vertex with passed vertex index: */
    inline int getVertexDegree(int vertexIndex) { return graphState->at(vertexIndex).second; };
    
    /* get indices of a vertices neighbours: */
    std::vector<int>* getNeighbours(int vertexIndex);
    std::vector<int>* getNeighbours(std::vector<int>* origins);

    std::pair<int, int>* getFirstValidEdge();

    std::vector<std::string>* getStringsFromVertexIndices(std::vector<int>* vertices);

    void printMappings(std::vector<int>* vertices);
    void printMappings();

};

#endif