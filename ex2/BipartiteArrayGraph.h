#ifndef BIPARTITEARRAYGRAPH_H
#define BIPARTITEARRAYGRAPH_H

#include "ArrayGraph.h"

class BipartiteArrayGraph : public ArrayGraph
{
//variables
public:

private:

    std::vector<int> pairU, pairV, dist;

//functions
public:
    static BipartiteArrayGraph* createBipartiteGraphByVertexSplit(ArrayGraph* G);

    inline BipartiteArrayGraph() {};

    inline int getLeftSize() { return getNumberOfVertices()/2; }
    inline int getRightSize() { return getNumberOfVertices()/2; }

    inline int getMaximumMatching() { return hopcroftKarp(); };
    inline int getMaximumMatchingCycleBound() { return hopcroftKarpCycleBound(); };

    bool areADJ(int u, int v);

private:

    int hopcroftKarp();
    int hopcroftKarpCycleBound();
    bool hasAugmentingPath_BFS();
    bool hasAugmentingPath_DFS(int u);
};

#endif