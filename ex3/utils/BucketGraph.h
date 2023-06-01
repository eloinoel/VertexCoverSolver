#ifndef BUCKETGRAPH_H
#define BUCKETRAPH_H

#include <string>
#include <vector>
#include <queue>
#include <utility>
#include <functional>
#include <limits.h>


#include <unordered_map> //O(1) for insert and access instead of O(log n) for ordered maps

#include <algorithm>

#include "boost/intrusive/list.hpp"

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
public:
    list_member_hook<> member_hook_; //for adj_refs

private:
    std::vector<int>* adj;
    //list<Vertex, member_hook<Vertex, list_member_hook<>, &Vertex::member_hook_>> adj_refs; //for O(1) check if has edge to specific vertex
    std::unordered_map<int, bool>* adj_map;

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

    /* used for reading in data, maps from original vertex name from input data to index and degree */
    std::unordered_map<std::string, std::pair<int, int>> originalVertexNames;

//functions
public:
    inline BucketGraph() {  }

    /* creates and initialises a graph from standard input */
    static BucketGraph* readStandardInput();
    std::vector<std::string>* getStringsFromVertexIndices(std::vector<int>* vertices);

    bool vertexHasEdgeTo(int vertex, int secondVertex); //O(1)
    bool isAdjMapConsistent();

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

private:

    //------------------------ Graph Construction ------------------------

    static std::string eraseLeadingTrailingWhitespacesFromString(std::string str);
    /* tests whether a char fulfills vertex naming format*/
	static bool isVertexCharacter(char c);

    void initActiveList(std::vector<std::pair<std::string, std::string>> edges); //--|
    void initAdjMap();                                                          //----> should be called in this order
    void initBucketQueue();                                                      //--|
    
    //-------------------------- Graph Utility --------------------------

    void addToBucketQueue(int degree, std::vector<BucketVertex*> vertices);
    void removeFromBucketQueue(int degree, std::vector<BucketVertex*> vertices);

    //------------------------ Virtual Flow Graph ------------------------
    /* create a mapping between the virtual flow graph and the bucket graph
    * in order to simulate the flow graph without instantiating it.
    * The mapping is in ascending order -> left vertices | right vertices | source | target
    */
    // Note: if the number of vertices changes all previously calculated mappings are invalid
    // Thus storing mapped indices for a period of time where this may occur is not very wise

    int virtualToRealIndex(int virtualIndex) {
        if (virtualIndex > (int) vertexReferences.size()*2-1) {
            throw std::invalid_argument("virtualToRealIndex: no real counterpart vertex exists!");
        }
        return virtualIndex % vertexReferences.size();
    }

    int realToVirtualIndex(int realIndex, bool toLeftSide) {
        if (realIndex < 0 || realIndex > (int) vertexReferences.size()-1) {
            throw std::invalid_argument("realToVirtualIndex: no real vertex with this index exists!");
        }
        if(!toLeftSide) {
            return realIndex+vertexReferences.size();
        }
        return realIndex;
    }

    std::vector<int>* getVirtualBipartitionEdges(int virtualIndex)
    {
        return {};
    }

    std::vector<int>* getVirtualFlowEdges(int virtualIndex)
    {
        return {};
    }

    bool matchingBFS(std::vector<int>* pairU, std::vector<int>* pairV, std::vector<int>* dist, int NIL);
    bool matchingDFS(std::vector<int>* pairU, std::vector<int>* pairV, std::vector<int>* dist, int u, int NIL);
    int hopcroftKarpMatchingSize();
    std::vector<int>* hopcroftKarpMatching();

    int edmondsKarp()
    {
        return -1;
    }

    //------------------------------ Bounds ------------------------------

    int getCliqueBound(int k = INT_MAX);
    bool vertexCanBeAddedToClique(int vertex, std::vector<int>* clique);

    int getLPBound();
    int getLPCycleBound();

    //------------------------ Data Reduction ------------------------
    //TODO: apply data reduction to input graph and return output graph

};

#endif