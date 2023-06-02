#include "Reductions.h"
#include "boost/intrusive/list.hpp"
#include "BucketGraph.h"


using namespace boost::intrusive;

bool Reductions::rule_HighDegree(BucketGraph* G, int* k)
{
    if(!(G->getMaxDegree() > *k)) return false; //cannot apply rule
    
    Reduction* reduction = new Reduction(RULE::HIGH_DEGREE, 0, new std::vector<int>());
    appliedRules->push_back(reduction);

    //delete vertices that have to be in the vertex cover
    while(G->getMaxDegree() > *k)
    {
        if(*k == 0) return true; //cannot delete more vertices, no possible vertex cover exists
        int maxDegVertex = G->getMaxDegreeVertex();
        reduction->kDecrement++;
        reduction->deletedVertices->push_back(maxDegVertex);
        *k = *k - 1;
        G->setInactive(maxDegVertex);
    }
    return true;
}

bool Reductions::rule_DegreeZero(BucketGraph* G)
{
    list<BucketVertex>* degZeroBucket = G->getVerticesOfDegree(0);

    if(degZeroBucket == nullptr || degZeroBucket->empty()) return false;

    Reduction* reduction = new Reduction(RULE::DEGREE_ZERO, 0, new std::vector<int>());
    appliedRules->push_back(reduction);

    for(auto it = degZeroBucket->begin(); it != degZeroBucket->end(); ++it)
    {
        reduction->deletedVertices->push_back(it->index);
        G->setInactive(it->index);
    }
    return true;
}

bool Reductions::rule_Buss(BucketGraph* G, int* k, int numVertices, int numEdges)
{
    int k_square = std::pow((*k), 2);

    if((numVertices > k_square + *k) || (numEdges > k_square))
        return true;
    return false;
}

bool Reductions::rule_DegreeOne(BucketGraph* G, int* k)
{
    list<BucketVertex>* degOneBucket = G->getVerticesOfDegree(1);

    if(degOneBucket == nullptr || degOneBucket->empty()) return false;

    Reduction* reduction = new Reduction(RULE::DEGREE_ONE, 0, new std::vector<int>());
    appliedRules->push_back(reduction);

    for(auto it = degOneBucket->begin(); it != degOneBucket->end(); ++it)
    {
        if(*k == 0) return true; //cannot delete more vertices, no possible vertex cover exists
        int neighbourToDelete = G->getVertex(it->index)->adj->front();
        reduction->deletedVertices->push_back(it->index);
        reduction->kDecrement++;
        *k = *k - 1;
        G->setInactive(it->index);
    }
    return true;
}