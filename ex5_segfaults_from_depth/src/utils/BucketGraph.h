#ifndef BUCKETGRAPH_H
#define BUCKETRAPH_H

#include <string>
#include <vector>
#include <queue>
#include <stack>
#include <utility>
#include <functional>
#include <limits.h>


#include <unordered_map> //O(1) for insert and access instead of O(log n) for ordered maps
#include <algorithm>

#include <boost/intrusive/list.hpp>
#include <boost/functional/hash.hpp>

class Reductions;
class Reduction;

using namespace boost::intrusive;
using namespace boost;

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
    intrusive::list<BucketVertex> vertices;
private:
    intrusive::list<BucketVertex>::const_iterator stable_iterator;

public:
    Bucket(int _degree, std::vector<BucketVertex*> _vertices)
     : degree(_degree)
    {
        vertices = intrusive::list<BucketVertex>();
        for (BucketVertex* vertex : _vertices)
        {
            vertices.push_back(*vertex);
        }
        stable_iterator = vertices.iterator_to(*vertices.begin());
    }

    inline void insert(BucketVertex* vertex)
    {
        vertices.push_back(*vertex);
    }

    inline void remove(BucketVertex* vertex)
    {
        if(vertices.iterator_to(*vertex)->index == stable_iterator->index) { ++stable_iterator; }
        vertices.erase(vertices.iterator_to(*vertex));
    }

    inline void insert(std::vector<BucketVertex*> _vertices)
    {
        for(BucketVertex* vertex : _vertices)
        {
            vertices.push_back(*vertex);
        }
    }

    inline void remove(std::vector<BucketVertex*>* _vertices)
    {
        for(BucketVertex* vertex : *_vertices)
        {
            if(vertices.iterator_to(*vertex)->index == stable_iterator->index) { ++stable_iterator; }
            auto it = vertices.iterator_to(*vertex);
            vertices.erase(it);
        }
    }

    /*  The stable iterator allows for deletion from-, and insertion into a bucket, while iterating through it
    *   Whenever the element, the iterator points to is deleted, the iterator is incremented
    *   When a new element is inserted into the bucket during iteration, the iterator will iterate over it later
    */
    inline intrusive::list<BucketVertex>::const_iterator* getStableIterator()
    {
        stable_iterator = vertices.iterator_to(*vertices.begin());
        return &stable_iterator;
    }
};

class Vertex : public list_base_hook<>
{
public:
    list_member_hook<> member_hook_; //for adj_refs

private:
    std::unordered_map<int, bool>* adj_map;

    std::string strName;
    int index;

    //temporary fields
    bool isActive;
    int degree;
    //pointer to bucketQueue entry for fast deletion
    BucketVertex* bucketVertex;

    friend class BucketGraph;
    friend class Reductions;
public:

    inline std::unordered_map<int, bool>* getAdj() { return adj_map; }
    inline bool getActive() { return isActive; }
    inline int getDegree() { return degree; }
    inline int getIndex() { return index; }


    Vertex(std::unordered_map<int, bool>* adjVertices, std::string orginalName, int _index, int startingDegree)
     : adj_map(adjVertices), strName(orginalName), index(_index), degree(startingDegree)
    {
        isActive = true;
        bucketVertex = new BucketVertex(index);
    }
};


class BucketGraph
{
//variables
public:
    // Size n
    std::vector<int>* dominationHelper;

    bool deg3ind = false;
    bool deg3clique = false;
    bool deg3dom = false;
    bool deg3 = false;
    bool deg4clique = false;
    bool deg4dom = false;

    bool deg3domDeg2 = false;

    int cntDeg3Ind = 0;
    int cntDeg3Clique = 0;
    int cntDeg3Dom1 = 0;
    int cntDeg3Dom2 = 0;
    int cntDeg4Clique = 0;

    int period_deg3 = 1;
    int period_unc = 3;
    int period_lp = 2;
    bool changePeriodDeg3 = false;

    int min_period = 1;

    int upperBound = 1;
    int upperBound1_4;
    int upperBound1_3;
    int upperBound1_2;

    bool LP_INITIALISED = false;
    bool UNCONFINED_INITIALISED = false;

private:
    /* each index represents a vertex, that maps to a node object that may be contained in the activeList */
    std::vector<Vertex*> vertexReferences;

    /* doubly linked list, that acts as a list of active vertices, O(1) access, deletion and insertion */
    //intrusive::list<Vertex> activeList;

    int numEdges;
    int numVertices;

    Reductions* reductions;

    /* each index represents a degree, that maps to a Bucket object that may be contained in the bucketQueue */
    std::vector<Bucket*> bucketReferences;
    /* priority queue of buckets that contain vertices of a certain degree (buckets are ordered after their degree ascendingly from front() to back()) */
    intrusive::list<Bucket> bucketQueue;
    /* used to realize iteration over bucketQueue while concurrently modifying it */
    intrusive::list<Bucket>::const_iterator stable_bucketQueue_inc_iterator;
    intrusive::list<Bucket>::const_iterator stable_bucketQueue_dec_iterator;

    /* used for reading in data, maps from original vertex name from input data to index and degree */
    std::unordered_map<std::string, std::pair<int, int>> originalVertexNames;
    std::vector<std::pair<std::string, std::string>> edges;

    /* Current Matching */
    std::vector<int>* pairU;
    std::vector<int>* pairV;
    std::vector<int>* dist;
    std::vector<int>* unmatched = nullptr;
    std::vector<int>* next_unmatched = nullptr;
    int NIL;
    int currentLPBound;
    bool didInitialMatchingCalculation = false;
    int nv;

    /* Unconfined */
    std::vector<bool>* mayBeUnconfined;

//functions
public:
    inline BucketGraph() {  }
    ~BucketGraph() { freeGraph(); }

    /* should only be called in unreduced graph */
    BucketGraph* copy();

    inline int getVertexReferencesSize() { return vertexReferences.size(); };
    void printReductionStack();

    void initRulePeriods(){
        period_deg3 = 16;
        period_unc = 32;
        period_lp = 25;
    }
//    void initRulePeriods(){
//        upperBound1_4 = upperBound/6;
//        if(upperBound1_4 == 0) upperBound1_4++;
//        period_deg3 = upperBound1_4;
//        upperBound1_2 = upperBound/3;
//        if(upperBound1_2 == 0) upperBound1_2++;
//        period_unc = upperBound1_2;
//        upperBound1_3 = upperBound/4;
//        if(upperBound1_3 == 0) upperBound1_3++;
//        period_lp = upperBound1_3;
////        period_rule = upperBound1_4;
//    }

    /* creates and initialises a graph from standard input */
    static BucketGraph* readStandardInput(bool initReductionDataStructures = true);
    std::vector<std::string>* getStringsFromVertexIndices(std::vector<int>* vertices);
    std::vector<std::string>* getStringsFromVertexIndices(std::unordered_map<int, bool>* vertices);
    /* creates a graph from the current graph and resets its data structures */
    BucketGraph* resetGraph();
    void freeGraph();

    bool vertexHasEdgeTo(int vertex, int secondVertex); //O(1)
    int getNumConnectedVertices();
    int getTotalNumVertices();
    int getNumVertices();
    int getNumEdges();

    void setActive(int vertexIndex);
    void setActive(std::vector<int>* vertexIndices);
    void setInactive(int vertexIndex);
    void setInactive(std::vector<int>* vertexIndices);
    bool isActive(int vertexIndex){ return vertexReferences[vertexIndex]->isActive;};

    std::vector<int>* getNeighbours(int vertexIndex);
    /* specify n starting with 0 */
    int getNthActiveNeighbour(int vertex, int n);
    std::pair<int, int> getFirstTwoActiveNeighbours(int vertex);

    bool containsConnectedVertex();
    int getMaxDegree();
    inline intrusive::list<Bucket>* getBucketQueue() { return &bucketQueue; };
    inline Bucket* getBucket(int degree) { return bucketReferences[degree]; };
    int getMaxDegreeVertex();
    int getRandomMaxDegreeVertex(int randomRangeCap = -1);
    int getRandomConnectedVertex(int randomRangeCap = -1);
    /* returns min degree vertex of degree > 0 and -1 if doesn't exist */
    int getMinDegreeVertex();
    int getVertexDegree(int vertexIndex);
    intrusive::list<BucketVertex>* getVerticesOfDegree(int degree);
    //inline intrusive::list<Vertex>* getActiveList() { return &activeList; }
    /* returns -1 if no vertex of degree */
    int getFirstVertexOfDegree(int degree);
    inline Vertex* getVertex(int index) { if(index < (int) vertexReferences.size()) return vertexReferences[index]; else return nullptr; }

    int getBucketSize(int degree){return (int)getBucket(3)->vertices.size();};
    /*  The stable iterator allows for deletion from-, and insertion into the bucketQueue, while iterating through it
    *   Whenever the element, the iterator points to is deleted, the iterator is incremented/decremented
    *   When a new element is inserted into the bucket during iteration, the iterator will iterate over it later
    *   USAGE: Whenever you want to use the stable iterator, call getStableBucketQueue[Inc/Dec]Iterator,
    *   depending on which way you want to iterate through the bucketQueue Inc -> incrementing, Dec -> decrementing
    *   NOTE: please only iterate the stable iterator in one loop at a time. i.e. nesting loops that each use the same type of stable iterator will not function properly
    */
    inline intrusive::list<Bucket>::const_iterator* getStableBucketQueueIncIterator()
    {
        stable_bucketQueue_inc_iterator = bucketQueue.iterator_to(*bucketQueue.begin());
        return &stable_bucketQueue_inc_iterator;
    }
    inline intrusive::list<Bucket>::const_iterator* getStableBucketQueueDecIterator()
    {
        stable_bucketQueue_dec_iterator = --bucketQueue.iterator_to(*bucketQueue.end());
        return &stable_bucketQueue_dec_iterator;
    }

    void print();
    void printActiveList();
    void printBucketQueue();
    void printBucketSizes();
    void printVertices(std::vector<int>* vertices);
    void printVertices(std::unordered_map<int, bool>* vertices);
    std::vector<std::string>* getEdgesToConsoleString();
    std::vector<std::string>* getOriginalEdgesToConsoleString();
    int getOriginalEdgeCount();
    void printMatching();
    void printVC(std::unordered_map<int, bool>* vc);

    int getLowerBoundVC();
    int getCliqueBound(int k = INT_MAX);
    int getLPBound();

    /* apply initial data reduction rules to graph */
    void preprocess(int* k, bool printDebug = false);
    /* apply initial data reduction rules to graph and possibly omit certain rules, 0: deg1, 1: deg2, 2: domination, 3: LP, 4: unconfined etc. */
    void preprocess(int* k, std::vector<bool>& rulesToApply, bool printDebug = false);

    /* apply data reduction rules to graph depending on search depth, returns true if no vertex cover can be found for this k */
    bool dynamicReduce(int* k, int depth, bool printDebug = false);
    /* apply data reduction rules to graph, returns true if no vertex cover can be found for this k */
    bool reduce(int* k, int depth, std::vector<bool>* rulesToApply = nullptr, bool printDebug = false);
    /* vc is not nullptr, if deleted vertices should be appended to vc*/
    void unreduce(int* k, int previousK, int depth, std::unordered_map<int, bool>* vc = nullptr);
    /* merge three vertices into one for degree 2 rule, returns vertex that was merged into and its previous adjacency list */
    std::tuple<int, std::unordered_map<int, bool>*, std::vector<int>*>* merge(int v0, int v1, int v2);
    /* restores previous previously merged vertices into 3 seperate vertices */
    void unmerge(Reduction* mergeRule);

    int getReductionStackSize();
    void verifyAdjacency();
    
    /* returns true if an edge (other vertex) could be added and false if it exists already
     * throws invalid_argument exception if faulty args were provided
     */
    bool addEdgeToVertex(int vertex, int edge);
    /* delete an edge from a vertex' adj list and map if exists */
    void removeEdgeFromVertex(int vertex, int edge);//TODO: optimise graph data structure to only contain adj maps instead of lists + maps

    void removeFromMatching(int vertexIndex);
    void queueForMatching(int vertexIndex);
    void setBipartMatchingFlowComponentsInactive(std::vector<int>* L, std::vector<int>* R, int k, double maxExecTime);
    int hopcroftKarpMatchingSize();
    inline void resetMatching() { void freeMatching(); void initMatching(); };
    inline void throwMatchingInconsistency() {
        for(int i=0; i<(int)pairU->size(); i++)
        {
            if(pairU->at(i) == NIL) continue;
            if(vertexReferences[i]->degree == 0) {
                throw std::invalid_argument("Vertex of degree 0 is matched");
            }
            if(i != pairV->at(pairU->at(i))) {
                throw std::invalid_argument("Vertex of matching is not symmetric");
            }
        }
    }

    inline bool isVertexScheduledForUnconfined(int vertexIndex) { return (*mayBeUnconfined)[vertexIndex]; }
    inline void scheduleForUnconfined(int vertexIndex) { (*mayBeUnconfined)[vertexIndex] = true; }
    inline void unscheduleForUnconfined(int vertexIndex) { (*mayBeUnconfined)[vertexIndex] = false; }
    void scheduleComponentForUnconfined(int vertexIndex);

    std::vector<std::pair<std::string,std::string>> getPreprocessedEdges();

    void unreduceDeg3Ind(Reduction* rule, std::unordered_map<int, bool>* vc = nullptr);
    void unreduceDeg3Clique(Reduction* rule, std::unordered_map<int, bool>* vc = nullptr);
    void unreduceDeg3Dom(Reduction* rule, std::unordered_map<int, bool>* vc = nullptr);

    void unreduceDeg4Clique(Reduction* rule, int* k, std::unordered_map<int, bool>* vc = nullptr);
    void unreduceDeg4Dom(Reduction* rule, int* k, std::unordered_map<int, bool>* vc = nullptr);

private:

    //------------------------ Graph Construction ------------------------

    static std::string eraseLeadingTrailingWhitespacesFromString(std::string str);
    /* tests whether a char fulfills vertex naming format*/
	static bool isVertexCharacter(char c);

    void initVertexReferences();  //----> should be called in this order
    void initBucketQueue(); //--|
    void initMatching();
    void initUnconfined();

    void freeMatching();
    void freeUnconfined();

    void initDominationHelper(){ dominationHelper = new std::vector<int> (getNumVertices(), 0); };
    void clearDominationHelper(){ delete dominationHelper; }
    //-------------------------- Graph Utility --------------------------

    int bruteForceCalculateNumEdges();
    /* add vertex to bucket of the degree of the vertex. If no bucket of degree in queue, insert new bucket before given biggerBucketDegree*/
    void addToBucketQueue(int vertex);
    /*  remove vertex from bucket of the degree of the vertex, deletes empty buckets from the queue */
    void removeFromBucketQueue(int vertex);
    /* move vertex of degree to a bucket of degree + 1*/
    void moveToBiggerBucket(int degree, int vertex);
    /* move vertex of degree to a bucket of degree - 1*/
    void moveToSmallerBucket(int degree, int vertex);

    //------------------------ Virtual Flow Graph ------------------------
    /* create a mapping between the virtual flow graph and the bucket graph
    * in order to simulate the flow graph without instantiating it.
    * The mapping is in ascending order -> left vertices | right vertices | source | target
    */
    // Note: if the number of vertices changes all previously calculated mappings are invalid
    // Thus storing mapped indices for a period of time where this may occur is not very wise

    bool matchingBFS();
    bool matchingDFS(int u);
    std::pair<int, std::pair<std::vector<int>*, std::vector<int>*>> hopcroftKarpMatching();

    //------------------------------ Bounds ------------------------------

    bool vertexCanBeAddedToClique(int vertex, std::vector<int>* clique);

    //int getLPBound();
};

#endif
