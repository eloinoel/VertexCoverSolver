#ifndef BUCKETGRAPH_H
#define BUCKETRAPH_H

#include <string>
#include <vector>
#include <queue>
#include <functional>
#include <limits.h>


#include <unordered_map> //O(1) for insert and access instead of O(log n) for ordered maps

#include <algorithm>

#include "boost/intrusive/list.hpp"

/* Rule Codierung:
 *
 * */
class ReductionVertices{
public:
    // Indicating the rule
    int rule;
    int kDecrement;
    std::vector<int> deletedVertices;
    std::vector<int>* savedAdjacency;
};


using namespace boost::intrusive;

class BucketVertex : public list_base_hook<>
{
public:
    int index;
public:
    BucketVertex(int _index)
     : index(_index)
    { }
};

class Bucket : public list_base_hook<>
{
public:
    int degree;
    list<BucketVertex> vertices;

public:
    Bucket(int _degree, std::vector<BucketVertex*> _vertices)
     : degree(_degree)
    {
        vertices = list<BucketVertex>();
        
        for (BucketVertex* vertex : _vertices)
        {
            vertices.push_back(*vertex);
        }
    }

    inline void insert(BucketVertex* vertex)
    {
        vertices.push_back(*vertex);
    }

    inline void remove(BucketVertex* vertex)
    {
        vertices.erase(vertices.iterator_to(*vertex));
    }

    inline void insert(std::vector<BucketVertex*> _vertices)
    {
        for(BucketVertex* vertex : _vertices)
        {
            vertices.push_back(*vertex);
        }
    }

    inline void remove(std::vector<BucketVertex*> _vertices)
    {
        for(BucketVertex* vertex : _vertices)
        {
            auto it = vertices.iterator_to(*vertex);
            vertices.erase(it);
        }
    }
};

class Vertex : public list_base_hook<>
{
    std::vector<int>* adj;
    //list adj_refs; TODO:
    std::string strName;
    int index;

    //temporary fields
    bool isActive;
    int degree;
    //maybe also add pointer to bucket queue ref for fast deletion
    BucketVertex* bucketVertex;

    friend class BucketGraph;
public:
    Vertex(std::vector<int>* adjVertices, std::string orginalName, int _index, int startingDegree)
     : adj(adjVertices), strName(orginalName), index(_index), degree(startingDegree)
    {
        isActive = true;
        bucketVertex = new BucketVertex(index);
    }
};


class BucketGraph
{
//variables
private:
    /* each index represents a vertex, that maps to a node object that may be contained in the activeList */
    std::vector<Vertex*> vertexReferences;

    /* doubly linked list, that acts as a list of active vertices, O(1) access, deletion and insertion */
    list<Vertex> activeList;

    /* each index represents a degree, that maps to a Bucket object that may be contained in the bucketQueue */
    std::vector<Bucket*> bucketReferences;
    /* priority queue of buckets that contain vertices of a certain degree (buckets are ordered after their degree ascendingly from front() to back()) */
    list<Bucket> bucketQueue;
    std::pair<int, int> degreeLimits; //bucketQueue empty: minDeg = INT32_MAX, maxDeg = -1

    /* used for reading in data, maps from original vertex name from input data to index and degree */
    std::unordered_map<std::string, std::pair<int, int>> originalVertexNames;

    /* Initialized to 0 and of size Vertex*/
    std::vector<int> neighbourArray;

//functions
public:
    inline BucketGraph() {  }

    /* creates and initialises a graph from standard input */
    static BucketGraph* readStandardInput();
    std::vector<std::string>* getStringsFromVertexIndices(std::vector<int>* vertices);

    void setActive(int vertexIndex);
    void setActive(std::vector<int>* vertexIndices);
    void setInactive(int vertexIndex);
    void setInactive(std::vector<int>* vertexIndices);

    std::vector<int>* getNeighbours(int vertexIndex);

    int getMaxDegree();
    int getMaxDegreeVertex();
    int getVertexDegree(int vertexIndex);
    int getVerticesOfDegree(int degree); //TODO:

    void print();
    void printActiveList();
    void printBucketQueue();

    int getLowerBoundVC();

    // return bool indicating if no vertex cover possible
    bool applyReductionRules(int* k, std::vector<ReductionVertices>* reductionArray);

private:

    //------------------------ Graph Construction ------------------------

    static std::string eraseLeadingTrailingWhitespacesFromString(std::string str);
    /* tests whether a char fulfills vertex naming format*/
	static bool isVertexCharacter(char c);

    void initActiveList(std::vector<std::pair<std::string, std::string>> edges);
    void initBucketQueue();
    
    //------------------------ Graph Utility ------------------------

    void addToBucketQueue(int degree, std::vector<BucketVertex*> vertices);
    void removeFromBucketQueue(int degree, std::vector<BucketVertex*> vertices);

    //------------------------ Bounds ------------------------

    int getCliqueBound(int k = INT_MAX);
    bool vertexCanBeAddedToClique(int vertex, std::vector<int>* clique);

    int getLPBound();

    //------------------------ Data Reduction ------------------------
    //TODO: apply data reduction to input graph and return output graph
    //================================================================
    // Reduction Rules
    bool rule_HighDegree(int *k, std::vector<ReductionVertices>* reductionVertices);
    bool rule_DegreeZero(std::vector<ReductionVertices>* reductionArray);

    // TODO: Check that numberOfVertices & numberOfEdges are up to date
    bool rule_Buss(int* k);

    void rule_DegreeOne(int* k, std::vector<ReductionVertices>* reductionArray);
    void rule_DegreeTwo(int* k, std::vector<ReductionVertices>* reductionArray);
    void rule_Domination(int* k, std::vector<ReductionVertices>* reductionArray);
    //================================================================

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