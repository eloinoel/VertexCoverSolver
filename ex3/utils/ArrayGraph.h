#ifndef ARRAYGRAPH_H
#define ARRAYGRAPH_H

#include <string>
#include <unordered_map> //O(1) for insert and access instead of O(log n) for ordered maps
#include <vector>
#include <algorithm>
#include <list>

enum RULES{
    DEGREE_ZERO,        // = 0
    DEGREE_ONE,         // = 1
    DEGREE_TWO,         // = 2
    HIGH_DEGREE,        // = 3
    DOMINATION          // = 4
};

class ReductionVertices{
public:
    int rule;
    int kDecrement;
    std::vector<int> deletedVertices; // First idx is always to add in VC if(rule!=0)
    std::vector<int>* savedAdjacency;
};

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

    //================================================================
    // TODO: AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAH
    // Reduction Rules
    bool rule_HighDegree(int *k, std::vector<ReductionVertices>* reductionVertices);
    bool rule_DegreeZero(std::vector<ReductionVertices>* reductionArray);

    // TODO: Check that numberOfVertices & numberOfEdges are up to date
    bool rule_Buss(int* k);

    void rule_DegreeOne(int* k, std::vector<ReductionVertices>* reductionArray);
    void rule_DegreeTwo(int* k, std::vector<ReductionVertices>* reductionArray);
//    void rule_Domination(int* k);
    //================================================================

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

    // return bool indicating if no vertex cover possible
    bool applyReductionRules(int* k, std::vector<ReductionVertices>* reductionArray);

    // Adds the deleted vertices from the reduction rules to the vertex cover
    void addReducedVertices(std::vector<int>* S, std::vector<ReductionVertices>* reductionArray);

    // Restores the initial kernel problem
    void addBackReducedVertices(int *k, std::vector<ReductionVertices>* reductionArray);

    //TODO: ------------Functions that need implementation-------------

    /* Get list of vertices of degree d*/
    std::vector<int>* getVerticesDegree(int d);

    /* Add 2 Adjacency Lists*/
    std::vector<int>* putAdjacencyTogether(std::vector<int>* neigh1, std::vector<int>* neigh2);

    /* Changes Adjacency of vertex
     * update own & neighbours degree here!
     * => add a degree to new neighbours they will be decreased after again by setInactive in Degree Two Rule
     * */
    void setVertexAdjacency(int vertexIndex, std::vector<int>* savedAdjacency);

    /* Restore saved Adjacency of vertex
     * update neighbours degree will be done in setActive()
     * */
    void setVertexAdjacencyBack(int vertexIndex, std::vector<int>* savedAdjacency);



};

#endif