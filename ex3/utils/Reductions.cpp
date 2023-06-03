#include "Reductions.h"
#include "boost/intrusive/list.hpp"
#include "BucketGraph.h"
#include <iostream>
#include <chrono>


using namespace boost::intrusive;

RULE_APPLICATION_RESULT Reductions::rule_HighDegree(BucketGraph* G, int* k)
{
    if(!(G->getMaxDegree() > *k)) return INAPPLICABLE; //cannot apply rule
    
    Reduction* reduction = new Reduction(RULE::HIGH_DEGREE, 0, nullptr, new std::vector<int>());
    appliedRules->push_back(reduction);

    //std::cout << "HIGHDEG Culling vertices: {";
    //delete vertices that have to be in the vertex cover
    while(G->getMaxDegree() > *k)
    {
        if(*k == 0) return INSUFFICIENT_BUDGET; //cannot delete more vertices, no possible vertex cover exists
        int maxDegVertex = G->getMaxDegreeVertex();
        reduction->kDecrement++;
        reduction->deletedVCVertices->push_back(maxDegVertex);
        *k = *k - 1;
        G->setInactive(maxDegVertex);
        //std::cout << maxDegVertex << ", ";
    }
    //std::cout << "}" << std::endl;
    return APPLICABLE;
}

RULE_APPLICATION_RESULT Reductions::rule_DegreeZero(BucketGraph* G)
{
    list<BucketVertex>* degZeroBucket = G->getVerticesOfDegree(0);

    if(degZeroBucket == nullptr || degZeroBucket->empty()) return INAPPLICABLE;

    //G->print();

    Reduction* reduction = new Reduction(RULE::DEGREE_ZERO, 0, new std::vector<int>());
    appliedRules->push_back(reduction);

    //std::cout << "DEGZERO Culling vertices: {";
    while(!degZeroBucket->empty())
    {
        auto first = degZeroBucket->begin();
        reduction->deletedVertices->push_back(first->index);
        G->setInactive(first->index);
        //std::cout << first->index << ", ";
    }
    //std::cout << "}" << std::endl;
    return APPLICABLE;
}

RULE_APPLICATION_RESULT Reductions::rule_Buss(BucketGraph* G, int* k, int numVertices, int numEdges)
{
    int k_square = std::pow((*k), 2);

    if((numVertices > k_square + *k) || (numEdges > k_square))
        return INSUFFICIENT_BUDGET;
    return INAPPLICABLE;
}

RULE_APPLICATION_RESULT Reductions::rule_DegreeOne(BucketGraph* G, int* k)
{
    list<BucketVertex>* degOneBucket = G->getVerticesOfDegree(1);

    if(degOneBucket == nullptr || degOneBucket->empty()) return INAPPLICABLE;

    Reduction* reduction = new Reduction(RULE::DEGREE_ONE, 0, nullptr, new std::vector<int>());
    appliedRules->push_back(reduction);

    while(!degOneBucket->empty())
    {
        if(*k == 0) return INSUFFICIENT_BUDGET; //cannot delete more vertices, no possible vertex cover exists
        auto it = degOneBucket->begin();
        int neighbourToDelete = G->getVertex(it->index)->adj->front();
        reduction->deletedVCVertices->push_back(it->index);
        reduction->kDecrement++;
        *k = *k - 1;
        G->setInactive(it->index);
    }
    return APPLICABLE;
}

RULE_APPLICATION_RESULT Reductions::rule_LPFlow(BucketGraph* G, int* k)
{
    std::vector<int>* delVertices = new std::vector<int>();
    std::vector<int>* delVCVertices = new std::vector<int>();
    //G->hopcroftKarpMatchingSize();
    //std::cout << "LPFlow" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    G->edmondsKarpFlow();
    //std::cout << "executed edmondsKarp" << std::endl;
    //auto start = std::chrono::high_resolution_clock::now();
    G->setBipartMatchingFlowComponentsInactive(delVertices, delVCVertices);
    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "Computed flow reduction with vc_verts=" << delVCVertices->size() << ", verts=" << delVertices->size() << " in " << (std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() / 1000) / (double) 1000 << " seconds" << std::endl;
    //std::cout << "determined strongly connected components" << std::endl;
    Reduction* reduction = new Reduction(RULE::LPFLOW, delVCVertices->size(), delVertices, delVCVertices);
    appliedRules->push_back(reduction);
    if((int) delVCVertices->size() == 0 && (int) delVertices->size() == 0)
    {
        //std::cout << "LPFlow: INAPPLICABLE" << std::endl;
        return INAPPLICABLE;
    }
    if((int) delVCVertices->size() > *k)
    {
        //std::cout << "LPFlow: INSUFFICENT_BUDGET" << std::endl;
        return INSUFFICIENT_BUDGET;
    }
    *k = *k - delVCVertices->size();
    return APPLICABLE;
}