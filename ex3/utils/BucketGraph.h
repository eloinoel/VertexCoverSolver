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
            vertices.erase(vertices.iterator_to(*vertex));
        }
    }
};

class Vertex : public list_base_hook<>
{
    std::vector<int>* adj;
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

    /* priority queue of buckets that contain vertices of a certain degree (buckets are ordered after their degree ascendingly from front() to back()) */
    list<Bucket> bucketQueue;

    /* used for reading in data, maps from original vertex name from input data to index and degree */
    std::unordered_map<std::string, std::pair<int, int>> originalVertexNames;

//functions
public:
    inline BucketGraph() {  }

    /* creates and initialises a graph from standard input */
    static BucketGraph* readStandardInput();

    void setActive(int vertexIndex); //TODO: buckets
    void setInactive(int vertexIndex); //TODO: buckets

    std::vector<int>* getNeighbours(int vertexIndex);

    int getMaxDegreeVertex(); //TODO:
    int getVertexDegree(int vertexIndex);
    int getVerticesOfDegree(int degree); //TODO:

    void print();
    void printBucketQueue(); //TODO:

    int getLowerBoundVC();

private:

    //------------------------ Graph Construction ------------------------

    static std::string eraseLeadingTrailingWhitespacesFromString(std::string str);
    /* tests whether a char fulfills vertex naming format*/
	static bool isVertexCharacter(char c);

    void initActiveList(std::vector<std::pair<std::string, std::string>> edges);
    void initBucketQueue();

    void addToBucketQueue(int degree, std::vector<BucketVertex*> vertices);

    //------------------------ Graph Utility ------------------------

    void addBucket(int degree, std::vector<BucketVertex*> vertices); // TODO: Use binary search
    void removeFromBucketQueue(int degree, std::vector<BucketVertex*> vertices);

    //------------------------ Bounds ------------------------

    int getCliqueBound(int k = INT_MAX);
    bool vertexCanBeAddedToClique(int vertex, std::vector<int>* clique);

    //------------------------ Data Reduction ------------------------
    //TODO: apply data reduction to input graph and return output graph

};

#endif