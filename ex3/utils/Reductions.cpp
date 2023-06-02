#include "Reductions.h"

bool Reductions::rule_HighDegree(BucketGraph* G, int k)
{
    ReductionVertices* rVertices = nullptr;
    if(G->getMaxDegree() > k)
    {
        ReductionVertices* rVertices = new ReductionVertices(RULE::HIGH_DEGREE, 0, new std::vector<int>());
        appliedRules->push_back(rVertices);
    }
    //delete vertices that have to be in the vertex cover
    while(G->getMaxDegree() > k)
    {
        int maxDegVertex = G->getMaxDegreeVertex();
        rVertices->kDecrement++;
        rVertices->deletedVertices->push_back(maxDegVertex);
        G->setInactive(maxDegVertex);
        k--;
    }

    if(rVertices == nullptr)
        return false;
    else
        return true;
}

bool Reductions::rule_DegreeZero()
{

}