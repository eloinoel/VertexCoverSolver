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

    //================================================================
    // For Cycle Bound
    int cycleNumber = 0;
    int minMax = 1;

    std::vector<std::vector<int>>* cycles;

    std::vector<int>* color;
    std::vector<int>* par;
    //================================================================

public:

//functions
private:
    inline ArrayGraph() {};

	static std::string eraseLeadingTrailingWhitespacesFromString(std::string str);

    /* tests whether a char fulfills vertex naming format*/
	static bool isVertexCharacter(char c);

    /* initially sets up activity flags & node degree data structure */
    void initGraphState(int vertexCount, int edgeCount);

    //================================================================
    // For Cycle Bound
    int getCycleBound();
    int getNaiveCycleBound();

    void dfs_cycle(int u, int p);

    void printCycles();

    std::list<int> getDisjointCycles();

    std::list<int> getDisjointSet(std::list<int> validCycles);

    std::list<int> getAllDisjointSet(std::list<int> validCycles);

    void printVector(std::list<int>* vec, std::string name);

    //================================================================

    bool vertexCanBeAddedToClique(int vertex, std::vector<int>* clique);
    int getCliqueBound();
    int partition(std::vector<int>* toSort, int low, int high);
    void quickSort(std::vector<int>* toSort, int low, int high);
    bool contains(std::vector<int>* vertexIndices, int vertexIndex);


public:
    bool isVertexCoverFound();

    // TODO: make more efficent, employ less copies
    static ArrayGraph* readStandardInput();

    inline int getVertexCount() { return adjacencyList.size(); }

    std::vector<int>* getVerticesSortedByDegree();

    void print();

    void printOriginalVertexNames();

    inline int getNumberOfVertices() { return numberOfVertices; }
    inline int getNumberOfEdges() { return numberOfEdges; }

    /* calculate lower bound for VC */
    int getLowerBoundVC();

    std::pair<int, int> getAllLowerBounds();

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
    /* int getVertexDegree(int vertexIndex) { // TODO: this is a preliminary implementation until degrees are updated on the fly
        int degree = 0;
        for (int i = 0; i < (int) adjacencyList[vertexIndex]->size(); i++)
        {
            if(graphState->at(adjacencyList[vertexIndex]->at(i)).first)
            {
                degree++;
            }
        }
        return degree;
    }; */
    
    /* get indices of a vertices neighbours: */
    std::vector<int>* getNeighbours(int vertexIndex);
    std::vector<int>* getNeighbours(std::vector<int>* origins);

    // Find the first component with size > 1 that contains any the origin points
    std::vector<int> getFirstComponent(std::vector<int>* origins);
    std::vector<std::vector<int>>* getComponents(std::vector<int>* origins);

    std::pair<int, int>* getFirstValidEdge();

    std::vector<std::string>* getStringsFromVertexIndices(std::vector<int>* vertices);

    void printMappings(std::vector<int>* vertices);
    void printMappings();

};

#endif