#include "Packing.h"
#include "BucketGraph.h"
#include <chrono>

void Packing::addNeighbourhoodConstraintForVertex(int vertex, int numActiveNeighbours)
{
    if(vertex > (int) neighbourhoodConstraint.size() - 1 || numActiveNeighbours <= 0) { throw std::invalid_argument("Packing: addNeighbourhoodConstraintForVertex: constraint values for vertex are illegal"); }

    neighbourhoodConstraint[vertex].first = 0;
    neighbourhoodConstraint[vertex].second = numActiveNeighbours;
}

bool Packing::incrementNeighbourhoodConstraintForVertex(int vertex)
{
    if(vertex > (int) neighbourhoodConstraint.size() - 1) { throw std::invalid_argument("Packing: incrementNeighbourhoodConstraintForVertex: given vertex is illegal."); }
    if(neighbourhoodConstraint[vertex].first < 0 || neighbourhoodConstraint[vertex].second < 0) { throw std::invalid_argument("Packing: incrementNeighbourhoodConstraintForVertex: constraint values for vertex are illegal"); }
    
    neighbourhoodConstraint[vertex].first++;
    if(neighbourhoodConstraint[vertex].first > neighbourhoodConstraint[vertex].second - 1)
        return false;
    else
        return true;
}

void Packing::decrementNeighbourhoodConstraintForVertex(int vertex)
{
    if(vertex > (int) neighbourhoodConstraint.size() - 1) { throw std::invalid_argument("Packing: decrementNeighbourhoodConstraintForVertex: given vertex is illegal."); }
    if(neighbourhoodConstraint[vertex].first < 0 || neighbourhoodConstraint[vertex].second < 0) { throw std::invalid_argument("Packing: incrementNeighbourhoodConstraintForVertex: constraint values for vertex are illegal"); }
    
    neighbourhoodConstraint[vertex].first--;
    if(neighbourhoodConstraint[vertex].first < 0)
    {
        throw std::invalid_argument("Packing: decrementNeighbourhoodConstraintForVertex: constraint value turned negative");
    }
}

void Packing::deleteNeighbourhoodConstraintForVertex(int vertex)
{
    if(vertex > (int) neighbourhoodConstraint.size() - 1) { throw std::invalid_argument("Packing: deleteNeighbourhoodConstraintForVertex: given vertex is illegal."); }
    
    neighbourhoodConstraint[vertex].first = -1;
    neighbourhoodConstraint[vertex].second = -1;
}

bool Packing::existsNeighbourhoodConstraintForVertex(int vertex)
{
    if(vertex > (int) neighbourhoodConstraint.size() - 1) { throw std::invalid_argument("Packing: deleteNeighbourhoodConstraintForVertex: given vertex is illegal."); }
    
    if(neighbourhoodConstraint[vertex].first < 0)
        return false;
    else
        return true;
}


bool Packing::checkNeighbourhoodConstraints(int branchingVertex, BucketGraph* G)
{
    int tmp = 0;
    auto startTime = std::chrono::high_resolution_clock::now();
    Vertex* vertexObj = G->getVertex(branchingVertex); return true;
    std::unordered_map<int, bool>* vertexAdjMap = vertexObj->getAdj();
    auto endTime = std::chrono::high_resolution_clock::now();
    
    bool currentConstraintsHold = true;
    std::vector<int> updatedConstraintVertices = std::vector<int>();
    int numActiveNeighbours = 0; //need this later for adding new constraint
    for(auto it = vertexAdjMap->begin(); it != vertexAdjMap->end(); ++it)
    {
        int neighbourOfVertex = it->first;
        if(G->isActive(neighbourOfVertex))
        {
            numActiveNeighbours++;
            if(existsNeighbourhoodConstraintForVertex(neighbourOfVertex))
            {
                currentConstraintsHold = currentConstraintsHold & incrementNeighbourhoodConstraintForVertex(neighbourOfVertex);
                updatedConstraintVertices.push_back(neighbourOfVertex);
                if(!currentConstraintsHold)
                {
                    break;
                }
            }
        }
    }
    if(!currentConstraintsHold) 
    {
        //revert changes made to constraints 
        for(auto it = updatedConstraintVertices.begin(); it != updatedConstraintVertices.end(); ++it)
        {
            decrementNeighbourhoodConstraintForVertex(*it);
        }
    }
    else
    {
        addNeighbourhoodConstraintForVertex(branchingVertex, numActiveNeighbours);
    }
    return currentConstraintsHold;
}

void Packing::undoNeighbourhoodConstraints(int branchingVertex, BucketGraph* G)
{
    deleteNeighbourhoodConstraintForVertex(branchingVertex);

    Vertex* vertexObj = G->getVertex(branchingVertex);
    std::unordered_map<int, bool>* vertexAdjMap = vertexObj->getAdj();
    for(auto it = vertexAdjMap->begin(); it != vertexAdjMap->end(); ++it)
    {
        int neighbourOfVertex = it->first;
        if(G->isActive(neighbourOfVertex))
        {
            if(existsNeighbourhoodConstraintForVertex(neighbourOfVertex))
            {
                decrementNeighbourhoodConstraintForVertex(neighbourOfVertex);
            }
        }
    }
}

