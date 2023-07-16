#ifndef PACKING_H
#define PACKING_H

#include <stdio.h>
#include <vector>
#include <stdexcept>

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

    void addNeighbourhoodConstraintForVertex(int vertex, int numActiveNeighbours)
    {
        if(vertex > (int) neighbourhoodConstraint.size() - 1 || numActiveNeighbours <= 0) { throw std::invalid_argument("Packing: addNeighbourhoodConstraintForVertex: constraint values for vertex are illegal"); }
        
        neighbourhoodConstraint[vertex].first = 0;
        neighbourhoodConstraint[vertex].second = numActiveNeighbours;
    }

    /* returns true if packing constraint still holds, returns false if the branching should be cut */
    bool incrementNeighbourhoodConstraintForVertex(int vertex)
    {
        if(vertex > (int) neighbourhoodConstraint.size() - 1) { throw std::invalid_argument("Packing: incrementNeighbourhoodConstraintForVertex: given vertex is illegal."); }
        if(neighbourhoodConstraint[vertex].first < 0 || neighbourhoodConstraint[vertex].second < 0) { throw std::invalid_argument("Packing: incrementNeighbourhoodConstraintForVertex: constraint values for vertex are illegal"); }
        
        neighbourhoodConstraint[vertex].first++;
        if(neighbourhoodConstraint[vertex].first > neighbourhoodConstraint[vertex].second - 1)
            return false;
        else
            return true;
    }

    void decrementNeighbourhoodConstraintForVertex(int vertex)
    {
        if(vertex > (int) neighbourhoodConstraint.size() - 1) { throw std::invalid_argument("Packing: incrementNeighbourhoodConstraintForVertex: given vertex is illegal."); }
        if(neighbourhoodConstraint[vertex].first < 0 || neighbourhoodConstraint[vertex].second < 0) { throw std::invalid_argument("Packing: incrementNeighbourhoodConstraintForVertex: constraint values for vertex are illegal"); }
        
        neighbourhoodConstraint[vertex].first--;
    }

    void deleteNeighbourhoodConstraintForVertex(int vertex)
    {
        if(vertex > (int) neighbourhoodConstraint.size() - 1) { throw std::invalid_argument("Packing: deleteNeighbourhoodConstraintForVertex: given vertex is illegal."); }
        
        neighbourhoodConstraint[vertex].first = -1;
        neighbourhoodConstraint[vertex].second = -1;
    }

    bool existsNeighbourhoodConstraintForVertex(int vertex)
    {
        if(vertex > (int) neighbourhoodConstraint.size() - 1) { throw std::invalid_argument("Packing: deleteNeighbourhoodConstraintForVertex: given vertex is illegal."); }
        
        if(neighbourhoodConstraint[vertex].first < 0)
            return false;
        else
            return true;
    }

};

#endif