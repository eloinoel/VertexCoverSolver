#ifndef PACKING_H
#define PACKING_H

#include <stdio.h>
#include <vector>
#include <stdexcept>

class BucketGraph;

class Packing
{

private:
    //Constraint: sum(x_u | u is neighbour of v) <= |N(v)| - 1; x_u is 1 if u in vertex cover, 0 otherwise
    //pair < sum(x_u), |N(v)| >
    std::vector<std::pair<int, int>> neighbourhoodConstraint;
public:

    Packing(int numGraphVertices) 
    {
        neighbourhoodConstraint = std::vector<std::pair<int, int>>(numGraphVertices);
        for(int i = 0; i < (int) neighbourhoodConstraint.size(); ++i)
        {
            neighbourhoodConstraint[i] = std::pair<int, int>({-1, -1});
        }
    }

    //**********************************************************
    //*************** Neightbourhood Constraint ****************
    //**********************************************************

    void addNeighbourhoodConstraintForVertex(int vertex, int numActiveNeighbours);

    /* returns true if packing constraint still holds, returns false if the branching should be cut */
    bool incrementNeighbourhoodConstraintForVertex(int vertex);

    void decrementNeighbourhoodConstraintForVertex(int vertex);

    void deleteNeighbourhoodConstraintForVertex(int vertex);

    bool existsNeighbourhoodConstraintForVertex(int vertex);

    /* returns true if constraints hold, returns false if solver should not branch on given vertex */
    bool checkNeighbourhoodConstraints(int branchingVertex, BucketGraph* G);

    void undoNeighbourhoodConstraints(int branchingVertex, BucketGraph* G);

};

#endif