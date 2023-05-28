#include <queue>
#include <math.h>
#include <climits>
#include <iostream>
#include "BipartiteArrayGraph.h"
#define NIL 0
#define INF INT_MAX

BipartiteArrayGraph* BipartiteArrayGraph::createBipartiteGraphByVertexSplit(ArrayGraph* G)
{
    std::vector<std::vector<int>*>* GADJ = G->getAdjacencyList();
    if(GADJ == nullptr)
    {
        return nullptr;
    }
    int rOffset = G->getNumberOfVertices();
    BipartiteArrayGraph* GP = new BipartiteArrayGraph();
    std::vector<std::vector<int>*> GPADJ = std::vector<std::vector<int>*>((rOffset)*2+1);

    // init dummy elements
    GPADJ[0] = new std::vector<int>();
    //GPADJ[rOffset+1] = new std::vector<int>();
    for(int i=0; i < (int) GADJ->size(); i++)
    {
        // push back l vertex and r vertex
        GPADJ[i+1] = new std::vector<int>();
        GPADJ[i+rOffset+1] = new std::vector<int>();
        for(int j=0; j < (int) GADJ->at(i)->size(); j++)
        {
            GPADJ[i + 1]->push_back(GADJ->at(i)->at(j) + rOffset + 1);
            GPADJ[i + rOffset + 1]->push_back(GADJ->at(i)->at(j) + 1);
        }
    }
    GP->setAdjacencyList(GPADJ);
    GP->initGraphState(G->getNumberOfVertices()*2, G->getNumberOfEdges()*2);
    return GP;
}


/* returns true if there is an augmenting path */
bool BipartiteArrayGraph::hasAugmentingPath_BFS()
{
    std::queue<int> Q; //an integer queue
 
    // First layer of vertices (set distance as 0)
    for (int u = 1; u < (int) pairU.size(); u++)
    {
        // If this is a free vertex, add it to queue
        if (pairU[u] == NIL)
        {
            // u is not matched
            dist[u] = 0;
            Q.push(u);
        }
 
        // Else set distance as infinite so that this vertex
        // is considered next time
        else dist[u] = INF;
    }
 
    // Initialize distance to NIL as infinite
    dist[NIL] = INF;
 
    // Q is going to contain vertices of left side only.
    while (!Q.empty())
    {
        // Dequeue a vertex
        int u = Q.front();
        Q.pop();
 
        // If this node is not NIL and can provide a shorter path to NIL
        if (dist[u] < dist[NIL])
        {
            // Get all adjacent vertices of the dequeued vertex u
            //std::vector<int>::iterator i;
            for (int i=0; i<(int) getAdjacencyList()->at(u)->size(); i++)
            {
                int v = getAdjacencyList()->at(u)->at(i) - getRightSize();
 
                // If pair of v is not considered so far
                // (v, pairV[V]) is not yet explored edge.
                if (dist[pairV[v]] == INF)
                {
                    // Consider the pair and add it to queue
                    dist[pairV[v]] = dist[u] + 1;
                    Q.push(pairV[v]);
                }
            }
        }
    }
 
    // If we could come back to NIL using alternating path of distinct
    // vertices then there is an augmenting path
    return (dist[NIL] != INF);
}

/* returns true if there is an augmenting path beginning with free vertex u */
bool BipartiteArrayGraph::hasAugmentingPath_DFS(int u)
{
    if (u != NIL)
    {
        //std::vector<int>::iterator i;
        for (int i=0; i<(int) getAdjacencyList()->at(u)->size(); i++)
        {
            // Adjacent to u
            int v = getAdjacencyList()->at(u)->at(i) - getRightSize();
 
            // Follow the distances set by BFS
            if (dist[pairV[v]] == dist[u]+1)
            {
                // If dfs for pair of v also returns
                // true
                if (hasAugmentingPath_DFS(pairV[v]) == true)
                {
                    pairV[v] = u;
                    pairU[u] = v;
                    return true;
                }
            }
        }
 
        // If there is no augmenting path beginning with u.
        dist[u] = INF;
        return false;
    }
    return true;
}

int BipartiteArrayGraph::hopcroftKarp()
{
    //stores pair of u in matching where u is a vertex on left side of pipartite graph, otherwise -1
    pairU = std::vector<int>(getLeftSize() + 1);
    //stores pair of v in matching where v is a vertex on rig side of pipartite graph, otherwise -1
    pairV = std::vector<int>(getRightSize() + 1);

    dist = std::vector<int>(getLeftSize() + 1);

    //init with no matching
    for(int i = 0; i < (int) pairU.size(); i++)
    {
        pairU[i] = NIL;
    }
    for(int i = 0; i < (int) pairV.size(); i++)
    {
        pairV[i] = NIL;
    }
    int result = 0;

    //keep updating result while there is an augmenting path
    while(hasAugmentingPath_BFS())
    {
        //find free vertex
        for(int i = 1; i < (int) pairU.size(); i++)
        {
            if(pairU[i] == NIL && hasAugmentingPath_DFS(i))
            {
                result++;
            }
        }
    }
    return result;
}

bool BipartiteArrayGraph::areADJ(int u, int v)
{
    for(int i=0; i<(int) getAdjacencyList()->at(u)->size(); i++)
    {
        if(getAdjacencyList()->at(u)->at(i) == v) return true;
    }
    return false;
}

int BipartiteArrayGraph::hopcroftKarpCycleBound()
{
    //stores pair of u in matching where u is a vertex on left side of pipartite graph, otherwise -1
    pairU = std::vector<int>(getLeftSize() + 1);
    //stores pair of v in matching where v is a vertex on rig side of pipartite graph, otherwise -1
    pairV = std::vector<int>(getRightSize() + 1);

    dist = std::vector<int>(getLeftSize() + 1);

    //init with no matching
    for(int i = 0; i < (int) pairU.size(); i++)
    {
        pairU[i] = NIL;
    }
    for(int i = 0; i < (int) pairV.size(); i++)
    {
        pairV[i] = NIL;
    }

    //keep updating result while there is an augmenting path
    while(hasAugmentingPath_BFS())
    {
        //find free vertex
        for(int i = 1; i < (int) pairU.size(); i++)
        {
            if(pairU[i] == NIL && hasAugmentingPath_DFS(i))
            {

            }
        }
    }

    /* std::cout << "pairU: {";
    for(auto elem : pairU)
        std::cout << elem << ", ";
    std::cout << "}" << std::endl; */

    int LPCyclebound = 0;
    std::vector<int> leftMatches = pairU;
    std::vector<int> currentCycle = std::vector<int>();
    bool foundCycle;
    int current;
    while(true)
    {
        // find uncovered vertex index
        current = -1;
        foundCycle = false;
        for(int i=1; i<(int)leftMatches.size(); i++)
        {
            if(leftMatches[i] != -1)
            {
                current = i;
                break;
            }
        }
        if(current == -1) break;
        //std::cout << "found uncovered vertex: " << current << "\n";

        // determine cycle
        while(current != -1 && current != 0/* leftMatches[current] != -1 && leftMatches[current] != 0 */)
        {
            if(currentCycle.size() >= 2 && contains(getAdjacencyList()->at(currentCycle.back()), currentCycle.front()+pairU.size()-1)/* currentCycle.front() == current */) // TODO:
            {
               /*  std::cout << "found cycle: {";
                for(auto elem : currentCycle)
                    std::cout << elem << ", ";
                std::cout << "}" << std::endl; */
                foundCycle = true;
                break;
            }
            //std::cout << "adding " << current << " to cycle\n";
            currentCycle.push_back(current);
            leftMatches[current] = -1;
            current = pairU[current];
        }

        if(foundCycle)
        {
            // attempt to subdivide the cycle
            /* if(currentCycle.size() >= 6)
            {
                for(int c=1; (int) currentCycle.size() >= c + 4; c += 2) {
                    for(int i=0; i<(int) currentCycle.size(); i++)
                    {
                        if(areADJ(currentCycle[i], currentCycle[i+3+c])
                        && areADJ(currentCycle[i+1], currentCycle[i+2+c]))
                        {
                            LPCyclebound += 1 + c;
                            // TODO: cull i+1 -> i+2+c from currentCycle
                            currentCycle.erase(currentCycle.begin()+i+1, currentCycle.begin()+i+2+c+1);
                            i -= (2 + c); 
                            if((int) currentCycle.size() < c + 4) break;
                        }
                    }
                }
            } */
            //std::cout << "found cycle of size: " << currentCycle.size() << "\n";
            // calculate 
            LPCyclebound += std::ceil((double) currentCycle.size() / (double) 2);
        }
        /* else
        {
            LPCyclebound += std::floor((double) currentCycle.size() / (double) 2);
        } */
        /* else if(currentCycle.size() == 2)
        {
            LPCyclebound++;
        } */
        currentCycle.clear();
    }
    return LPCyclebound;
}