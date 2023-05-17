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

    //================================================================
    // For Cycle Bound
    int cycleNumber = 0;

    int numberOfVertices;

    std::vector<std::vector<int>>* cycles;

    std::vector<int>* color;
    std::vector<int>* par;
    //================================================================

//functions
private:
    inline ArrayGraph() {};

	static std::string eraseLeadingTrailingWhitespacesFromString(std::string str);

    /* tests whether a char fulfills vertex naming format*/
	static bool isVertexCharacter(char c);

    /* initially sets up activity flags & node degree data structure */
    void initGraphState();

    //================================================================
    // For Cycle Bound
    int getCycleBound();

    void dfs_cycle(int u, int p);

    void printCycles();
    //================================================================


public:

    // TODO: make more efficent, employ less copies
    static ArrayGraph* readStandardInput();

    inline int getVertexCount() { return adjacencyList.size(); }

    void print();

    void printOriginalVertexNames();

    /* calculate lower bound for VC */
    int getLowerBoundVC();

    /* get flag array for vertices being active: */
    //std::vector<bool>* getActive();

    /* set vertices with the passed vertex indices to active and update degrees for those vertices and their neighbours: */
    inline void setActive(int vertexIndex) { graphState->at(vertexIndex).first = true; };
    void setActive(std::vector<int>* vertexIndices);

    /* set vertices with the passed vertex indices to inactive and update degrees for their neighbours: */
    inline void setInactive(int vertexIndex) { graphState->at(vertexIndex).first = false; };
    void setInactive(std::vector<int>* vertexIndices);

    inline bool isActive(int vertexIndex) { return graphState->at(vertexIndex).first; };

    /* get the indices of all vertices that are inactive */
    std::vector<int>* getInactiveVertices();

    /* get graph state array: */
    inline std::vector<std::pair<bool, int>>* getState() { return graphState; };

    /* get vertex with highest degree among active vertices:
    * returns -1 if no vertex in graph
    */
    int getMaxDegreeVertex();

    /* get degree of the vertex with passed vertex index: */
    //inline int getVertexDegree(int vertexIndex) { return graphState->at(vertexIndex).second; };
    int getVertexDegree(int vertexIndex) { // TODO: this is a preliminary implementation until degrees are updated on the fly
        int degree = 0;
        for (int i = 0; i < (int) adjacencyList[vertexIndex]->size(); i++)
        {
            if(graphState->at(adjacencyList[vertexIndex]->at(i)).first)
            {
                degree++;
            }
        }
        return degree;
    };
    
    /* get indices of a vertices neighbours: */
    std::vector<int>* getNeighbours(int vertexIndex);

    int getFirstActiveVertex();

    std::vector<std::string>* getStringsFromVertexIndices(std::vector<int>* vertices);

    void printMappings(std::vector<int>* vertices);

};

#endif