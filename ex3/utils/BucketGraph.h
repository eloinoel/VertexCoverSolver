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

#include "Reductions.h"

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
    Bucket(int _degree, std::vector<BucketVertex*>& _vertices)
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
    friend class Reductions;
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

    int numEdges;
    int numVertices;

    Reductions* reductions;

    /* each index represents a degree, that maps to a Bucket object that may be contained in the bucketQueue */
    std::vector<Bucket*> bucketReferences;
    /* priority queue of buckets that contain vertices of a certain degree (buckets are ordered after their degree ascendingly from front() to back()) */
    list<Bucket> bucketQueue;

    /* used for reading in data, maps from original vertex name from input data to index and degree */
    std::unordered_map<std::string, std::pair<int, int>> originalVertexNames;
    std::vector<std::pair<std::string, std::string>> edges;

    /* Current Matching */
    std::vector<int> pairU;
    std::vector<int> pairV;
    std::vector<int> dist;
    std::vector<int>* unmatched;
    std::vector<int>* next_unmatched;
    int NIL;
    int currentLPBound;
    bool didInitialMatchingCalculation = false;

    std::vector<std::vector<int>> flow;

//functions
public:
    inline BucketGraph() {  }

    /* creates and initialises a graph from standard input */
    static BucketGraph* readStandardInput();
    std::vector<std::string>* getStringsFromVertexIndices(std::vector<int>* vertices);
    /* creates a graph from the current graph and resets its data structures */
    BucketGraph* resetGraph();

    bool vertexHasEdgeTo(int vertex, int secondVertex); //O(1)
    int getNumVertices();
    int getNumEdges();

    void setActive(int vertexIndex);
    void setActive(std::vector<int>* vertexIndices);
    void setInactive(int vertexIndex);
    void setInactive(std::vector<int>* vertexIndices);

    std::vector<int>* getNeighbours(int vertexIndex);

    int getMaxDegree();
    int getMaxDegreeVertex();
    /* heuristic from paper which generally worsens performance a bit but reduces number of recursive steps */
    int getMaxDegreeVertexMinimisingNeighbourEdges();
    int getVertexDegree(int vertexIndex);
    list<BucketVertex>* getVerticesOfDegree(int degree);
    inline list<Vertex>* getActiveList() { return &activeList; }
    /* returns -1 if no vertex of degree */
    int getFirstVertexOfDegree(int degree);
    inline Vertex* getVertex(int index) { if(index < vertexReferences.size()) return vertexReferences[index]; else return nullptr; }

    void print();
    void printActiveList();
    void printBucketQueue();
    void printEdgesToConsole();
    void printMatching();

    int getLowerBoundVC();
    int getCliqueBound(int k = INT_MAX);
    int getLPBound();
    int getFlow();

    void resetLPBoundDataStructures();

    /* apply data reduction rules to graph */
    void reduce();

private:

    //------------------------ Graph Construction ------------------------

    static std::string eraseLeadingTrailingWhitespacesFromString(std::string str);
    /* tests whether a char fulfills vertex naming format*/
	static bool isVertexCharacter(char c);

    void initActiveList();  //--|
    void initAdjMap();      //----> should be called in this order
    void initBucketQueue(); //--|
    void initMatching();
    bool isAdjMapConsistent();

    //-------------------------- Graph Utility --------------------------

    int bruteForceCalculateNumEdges();
    void addToBucketQueue(int degree, std::vector<BucketVertex*> vertices);
    /* add vertices of specific degree to bucket queue. If no bucket of degree in queue, insert new bucket before given biggerBucketDegree*/
    //void addToBucketQueueBeforeBucket(int degree, std::vector<BucketVertex*> vertices, int biggerBucketDegree);
    //void addToBucketQueueBeforeBucket(int degree, std::vector<BucketVertex*> vertices, list_iterator<bhtraits<Bucket, list_node_traits<hook_defaults::void_pointer>, safe_link, hook_defaults::tag, 1U>, false> it);
    void removeFromBucketQueue(int degree, std::vector<BucketVertex*> vertices);
    void moveToBiggerBucket(int degree, int vertex);
    void moveToSmallerBucket(int degree, int vertex);

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

    bool matchingBFS();
    bool matchingDFS(int u);
    int hopcroftKarpMatchingSize();
    std::pair<int, std::pair<std::vector<int>*, std::vector<int>*>> hopcroftKarpMatching();

    int edmondsKarpFlow()
    {
        std::cout << ".";
        //std::cout << "Beginning flow" << std::endl;
        // number of vertices helper variable (only for clarity)
        const int nv = vertexReferences.size();

        int flow_amt = 0;
        int s = nv*2;
        int t = nv*2+1;
        std::vector<int> pred = std::vector<int>(nv*2+2);
        std::queue<int> Q = std::queue<int>();
        int current;

        // init flow
        // TODO: see how much of this can be saved across iterations
        flow = std::vector<std::vector<int>>(nv*2+2);
        for(int i=0; i<nv*2+2; i++)
        {
            flow[i] = std::vector<int>(nv*2+2);
            pred[i] = -1;
            for(int j=0; j<nv*2+2; j++)
            {
                flow[i][j] = 0;
            }
        }
        // set left to right flow according to matching
        for (int i=0; i<(int) pairU.size(); i++)
        {
            if(pairU[i] != NIL)
            {
                flow[s][i] = 1;
                flow[i][pairU[i]] = 1;
                flow[pairU[i]][t] = 1;
            }
        }
        //std::cout << "Initialized flow" << std::endl;
        pred[t] = 1;
        // Until no augmenting path was found last iteration
        while (pred[t] != -1)
        {
            // clear queue & predecessors
            while (!Q.empty()) { Q.pop(); }
            for(int i=0; i<nv*2+2; i++) { pred[i] = -1; }
            // init queue with source vertex
            Q.push(s);
            //std::cout << "Pushed source vertex" << std::endl;
            // BFS to find the shortest s-t path
            while (!Q.empty())
            {
                // retrieve next queued vertex
                current = Q.front();
                Q.pop();

                // if current is s (capacity == 1)
                if(current == s)
                {
                    //std::cout << "Evaluating source vertex" << std::endl;
                    // s exclusively has edges to left vertices (indices 0 to vertexReferences.size()-1)
                    for (int i=0; i<nv; i++)
                    {
                        if(!vertexReferences[i]->isActive) { continue; }
                        // if left vertex i wasn't expanded & edge to it isn't already in flow
                        if (pred[i] == -1 && 1 > flow[current][i])
                        {
                            pred[i] = current;
                            Q.push(i);
                        }
                    }
                    //std::cout << "Evaluated source vertex" << std::endl;
                }
                // if current is a left vertex (capacity == INF for edges to right vertices)
                else if (current < nv)
                {
                    //std::cout << "Evaluating left side vertex" << std::endl;
                    // current has edges to right vertices (indices vertexReferences.size() to vertexReferences.size()*2-1)
                    for (auto v=vertexReferences[current]->adj->begin(); v != vertexReferences[current]->adj->end(); ++v)
                    {
                        if(!vertexReferences[*v]->isActive) { continue; }
                        if (pred[nv+*v] == -1)
                        {
                            pred[nv+*v] = current;
                            Q.push(nv+*v);
                        }
                    }
                    //std::cout << "Evaluated left side vertex" << std::endl;
                }
                // if current is a right vertex (capacity == 1 for edge to t and capacity == INF for the reverse edges to left vertices)
                else if (nv <= current && current < nv*2)
                {
                    //std::cout << "Evaluating right side vertex: " << current << std::endl;
                    // current has reverse edges to left vertices (indices vertexReferences.size() to vertexReferences.size()*2-1)
                    for (auto v=vertexReferences[current-nv]->adj->begin(); v != vertexReferences[current-nv]->adj->end(); ++v)
                    {
                        if(!vertexReferences[*v]->isActive) { continue; }
                        //std::cout << "Evaluating right side vertices left side neighbour" << std::endl;
                        if (pred[*v] == -1 && pairU[*v] == current)
                        {
                            pred[*v] = current;
                            Q.push(*v);
                        }
                    }
                    //std::cout << "Evaluated right side vertices reverse edges" << std::endl;
                    if (pred[t] == -1 && 1 > flow[current][t])
                    {
                        pred[t] = current;
                        break;
                        //std::cout << "Found s-t path" << std::endl;
                    }
                }
            }
            //std::cout << "Finished DFS" << std::endl;

            if (pred[t] != -1)
            {
                std::cout << "Found flow-improving path" << std::endl;
                printMatching();
                /* int df = INT32_MAX;
                for (int p = t; p != -1 && p != s; p = pred[p])
                {
                    //std::cout << "pathstep" << std::endl;
                    int cap = INT32_MAX;
                    if(p == t || pred[p] == t) { cap = 1; }
                    df = std::min(df, cap - flow[pred[p]][p]);
                }
                //std::cout << "Determined df" << std::endl;
                for (int p = t; p != -1 && p != s; p = pred[p])
                {
                    flow[pred[p]][p] += df;
                    flow[p][pred[p]] -= df;
                }
                flow_amt += df;
                */
                // always: df = 1;
                for (int p = t; p != -1 && p != s; p = pred[p])
                {
                    flow[pred[p]][p]++;
                    flow[p][pred[p]]--;
                    if(pred[p] == s || p == t) { continue; }
                    // using reverse edge
                    if (p < nv)
                    {
                        pairV[pairU[p]] = NIL;
                        pairU[pairV[pred[p]-nv]] = NIL;

                        pairU[p] = NIL;
                        pairV[pred[p]-nv] = NIL;
                        std::cout << "-" << "(" << p << ", " << pred[p]-nv << ") ";
                    }
                    // using forward edge
                    else
                    {
                        pairV[pairU[pred[p]]] = NIL;
                        pairU[pairV[p-nv]] = NIL;

                        pairU[pred[p]] = p-nv;
                        pairV[p-nv] = pred[p];
                        std::cout << "+" << "(" << pred[p] << ", " << p-nv << ") ";
                    }
                }
                std::cout << std::endl;
                printMatching();
                flow_amt++;
                //std::cout << "Processed flow-improving path" << std::endl;
            }
        }
        //std::cout << "Ending flow" << std::endl;
        return flow_amt;
    }

    void LPReduce()
    {
        std::stack<int> S = std::stack<int>();
        std::vector<bool> visited = std::vector<bool>(pairU.size());
        bool isNotComp = false;
        // for each possible connected component (We exclude degree 0 connected components here)
        for (int i=0; i<(int) pairU.size(); i++)
        {
            if (!vertexReferences[i]->isActive || pairU[i] == NIL/* this culls deg=0 rule */) { continue; }
            for(int j=0; j<(int) visited.size(); j++) { visited[j] = false; }
            while(!S.empty()) { S.pop(); }
            isNotComp = false;
            S.push(i);
            std::vector<int> componentL = std::vector<int>();
            // until component is closed
            while (!S.empty())
            {
                int current = S.top();
                // add vertex to component, if visited
                if(visited[current])
                {
                    componentL.push_back(current);
                    S.pop();
                    continue;
                }

                // expand node
                visited[current] = true;
                for (auto v=vertexReferences[current]->adj->begin(); v != vertexReferences[current]->adj->end(); ++v)
                {
                    // if we find an unmatched right vertex, abort immediately (This cannot be a component)
                    if(pairV[*v] == NIL) { isNotComp = true; break; }
                    // We skip inactive vertices and the matched right vertex and vertices of which we already visited their left matched vertex
                    if(!vertexReferences[*v]->isActive || current == pairV[*v] || visited[pairV[*v]]) { continue; }
                    // push next left vertex
                    S.push(pairV[*v]);
                }
                if(isNotComp) { break; }
            }
            if(isNotComp) { continue; }
            // found component
            // TODO: iterate through componentL and calc component right vertices
            std::cout << "Components left vertices: {";
            for (int j=0; j<(int) componentL.size(); j++)
            {
                std::cout << componentL[j] << ", ";
            }
            std::cout << "}" << std::endl;
            std::cout << "Components right vertices: {";
            for (int j=0; j<(int) componentL.size(); j++)
            {
                for (auto v=vertexReferences[componentL[j]]->adj->begin(); v != vertexReferences[componentL[j]]->adj->end(); ++v)
                {
                    std::cout << *v << ", ";
                }
            }
            std::cout << "}" << std::endl;
            // TODO: delete vertices from graph via Reductions
        }
    }

    //------------------------------ Bounds ------------------------------

    bool vertexCanBeAddedToClique(int vertex, std::vector<int>* clique);

    //int getLPBound();
    int getLPCycleBound();  // TODO: this is still trash

    //------------------------ Data Reduction ------------------------
    //TODO: apply data reduction to input graph and return output graph



};

#endif