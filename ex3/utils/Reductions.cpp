#include "ColorPrint.h"
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

    //std::cout << "DEGONE Culling vertices: {";
    while(!degOneBucket->empty())
    {
        if(*k == 0) return INSUFFICIENT_BUDGET; //cannot delete more vertices, no possible vertex cover exists
        auto it = degOneBucket->begin();
        int neighbourToDelete = G->getNthActiveNeighbour(it->index, 0);
        reduction->deletedVCVertices->push_back(neighbourToDelete);
        reduction->kDecrement++;
        *k = *k - 1;
        G->setInactive(neighbourToDelete);
        //std::cout << neighbourToDelete << ", ";
    }
    //std::cout << "}" << std::endl;
    return APPLICABLE;
}

RULE_APPLICATION_RESULT Reductions::rule_LPFlow(BucketGraph* G, int* k)
{
    std::vector<int>* delVertices = new std::vector<int>();
    std::vector<int>* delVCVertices = new std::vector<int>();
    G->hopcroftKarpMatchingSize();
    //std::cout << "LPFlow" << std::endl;
    //auto start = std::chrono::high_resolution_clock::now();
    G->edmondsKarpFlow();
    //std::cout << "executed edmondsKarp" << std::endl;
    //auto start = std::chrono::high_resolution_clock::now();
    G->setBipartMatchingFlowComponentsInactive(delVertices, delVCVertices);
    //auto stop = std::chrono::high_resolution_clock::now();
    //std::cout << "Computed flow reduction with vc_verts=" << delVCVertices->size() << ", verts=" << delVertices->size() << " in " << (std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() / 1000) / (double) 1000 << " seconds" << std::endl;
    //std::cout << "determined strongly connected components" << std::endl;
    Reduction* reduction = new Reduction(RULE::LPFLOW, delVCVertices->size(), delVertices, delVCVertices);
    appliedRules->push_back(reduction);
    if((int) delVCVertices->size() == 0 && (int) delVertices->size() == 0)
    {
        //std::cout << "LPFlow: INAPPLICABLE" << std::endl;
        appliedRules->pop_back();
        return INAPPLICABLE;
    }
    if((int) delVCVertices->size() > *k)
    {
        //std::cout << "LPFlow: INSUFFICENT_BUDGET" << std::endl;
        /* G->setActive(delVertices);
        G->setActive(delVCVertices); */
        appliedRules->pop_back();
        return INSUFFICIENT_BUDGET;
    }
    *k = *k - delVCVertices->size();    // TODO: im getting core dumps, when restoring
    return APPLICABLE;                  // it should be correct to set vertices back to active, but this doesnt work
}

/*----------------------------------------------------------*/
/*------------------   Reduction Rules   -------------------*/
/*----------------------------------------------------------*/
//void Reductions::initRuleCounter()
//{
//    rule_0 = 0;
//    rule_1 = 0;
//    rule_2 = 0;
//    rule_3 = 0;
//    rule_4 = 0;
//    rule_5 = 0;
//}

void Reductions::printReductionRules()
{
    if(appliedRules->empty())
        return;

    for (auto i = 0; i < (int) appliedRules->size(); ++i) {
        Reduction* reductionR = appliedRules->at(i);
        std::string toPrint =  "Applied Reduction Rule Nr " + std::to_string(reductionR->rule) + ":";
        std::cout << ColorPrint::dye(toPrint, 'r') << std::endl ;

        std::string deleteV =  "Deleted Vertices are: ";
        for (int j = 0; j < (int) reductionR->deletedVertices->size(); ++j) {
            deleteV += std::to_string(reductionR->deletedVertices->at(j)) + ", ";
        }
        std::cout << ColorPrint::dye(deleteV, 'c') << std::endl ;

        if(reductionR->rule == DEGREE_TWO)
        {
            std::string savedAdj = "Saved Adjacency List: ";
            for (int k = 0; k < (int) std::get<1>(*reductionR->mergeVertexInfo)->size(); ++k)
            {
                savedAdj += std::to_string(std::get<1>(*reductionR->mergeVertexInfo)->at(k)) + ", ";
            }
            std::cout << ColorPrint::dye(savedAdj, 'c') << std::endl;
        }
        std::cout << std::endl;
    }
}

RULE_APPLICATION_RESULT Reductions::rule_DegreeTwo(BucketGraph* G, int* k)
{
    list<BucketVertex>* degTwoBucket = G->getVerticesOfDegree(2);

    if(degTwoBucket == nullptr || degTwoBucket->empty()) return INAPPLICABLE;

    std::cout << "DEGTWO Culling vertices: ";
    while(!degTwoBucket->empty())
    {
        if(*k == 0) return INSUFFICIENT_BUDGET; //cannot delete more vertices, no possible vertex cover exists

        auto it = degTwoBucket->begin();
        std::pair<int, int>* neighbours = G->getFirstTwoActiveNeighbours(it->index); //should always return valid vertices

        // save deleted vertex
        Reduction* delVer = new Reduction(RULE::DEGREE_TWO, 0, new std::vector<int>());

        //if no merge, take neighbours, otherwise case destinction whether merged vertex is in vc
        delVer->deletedVCVertices->push_back(neighbours->first);
        delVer->deletedVCVertices->push_back(neighbours->second);
        delVer->deletedVertices->push_back(it->index);

        std::cout << it->index  << ", " << neighbours->first << ", " << neighbours->second;

        delVer->mergeVertexInfo = nullptr;
        // CASE The neighbours know each other, take them into vc
        if(G->vertexHasEdgeTo(neighbours->first, neighbours->second))
        {
            delVer->kDecrement = 2;
            G->setInactive(delVer->deletedVertices);
            G->setInactive(delVer->deletedVCVertices);
            (*k) = (*k) - 2;
        }
        // CASE Neighbours don't know each other => setInactive or do that in addReducedVertices?
        else
        {
            delVer->kDecrement = 1;
            delVer->mergeVertexInfo = G->merge(it->index, neighbours->first, neighbours->second); //sets merged vertices inactive
            (*k) = (*k) - 1;
            std::cout << " | merging into " << std::get<0>(*delVer->mergeVertexInfo);
        }
        appliedRules->push_back(delVer);
    }
    return APPLICABLE;
}