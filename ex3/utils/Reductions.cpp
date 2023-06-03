#include "ColorPrint.h"
#include "Reductions.h"
#include "boost/intrusive/list.hpp"
#include "BucketGraph.h"
#include <iostream>

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
        if(*k == 0) return INSUFFIENT_BUDGET; //cannot delete more vertices, no possible vertex cover exists
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
    //G->print();
    return APPLICABLE;
}

RULE_APPLICATION_RESULT Reductions::rule_Buss(BucketGraph* G, int* k, int numVertices, int numEdges)
{
    int k_square = std::pow((*k), 2);

    if((numVertices > k_square + *k) || (numEdges > k_square))
        return INSUFFIENT_BUDGET;
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
        if(*k == 0) return INSUFFIENT_BUDGET; //cannot delete more vertices, no possible vertex cover exists
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

    G->edmondsKarpFlow();
    G->getBipartMatchingFlowComponents(delVertices, delVCVertices);
    if((int) delVCVertices->size() > *k)
    {
        return INSUFFIENT_BUDGET;
    }
    if((int) delVCVertices->size() == 0 && (int) delVertices->size() == 0)
    {
        return INAPPLICABLE;
    }
    G->setInactive(delVertices);
    G->setInactive(delVCVertices);
    Reduction* reduction = new Reduction(RULE::LPFLOW, delVCVertices->size(), delVertices, delVCVertices);
    appliedRules->push_back(reduction);
    *k = *k - delVCVertices->size();
    return APPLICABLE;
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
            for (int k = 0; k < (int) reductionR->savedAdjacency->size(); ++k)
            {
                savedAdj += std::to_string(reductionR->savedAdjacency->at(k)) + ", ";
            }
            std::cout << ColorPrint::dye(savedAdj, 'c') << std::endl;
        }
        std::cout << std::endl;
    }

}

RULE_APPLICATION_RESULT Reductions::rule_DegreeTwo(BucketGraph* G, int* k)
{
    list<BucketVertex>* degTwoBucket = G->getVerticesOfDegree(0);

    if(degTwoBucket == nullptr || degTwoBucket->empty()) return INAPPLICABLE;

    Reduction* reduction = new Reduction(RULE::DEGREE_TWO, 0, nullptr, new std::vector<int>());
    appliedRules->push_back(reduction);

    while(!degTwoBucket->empty())
    {
        if(*k == 0) return INSUFFIENT_BUDGET; //cannot delete more vertices, no possible vertex cover exists

        auto it = degTwoBucket->begin();
        int vertDegTwo = G->getVertex(it->index)->adj->front();

        std::vector<int>* neighbours = G->getNeighbours(vertDegTwo);

        int neighbourOne = neighbours->at(0);
        int neighbourTwo = neighbours->at(1);

        // save deleted vertex
        Reduction* delVer = new Reduction(RULE::DEGREE_TWO, 0, new std::vector<int>());
        delVer->deletedVertices->push_back(neighbourOne);            // at (0)
        delVer->deletedVertices->push_back(neighbourTwo);     // at (1)
        delVer->deletedVertices->push_back(vertDegTwo);                     // at (2)

        delVer->savedAdjacency = neighbours;

        // CASE The neighbours know each other
        if(G->vertexHasEdgeTo(neighbourOne,neighbourTwo))
        {
            delVer->kDecrement = 2;

            G->setInactive(vertDegTwo);
            (*k)-=2;
        }
        // CASE Neighbours don't know each other => setInactive or do that in addReducedVertices?
        else
        {
            delVer->kDecrement = 1;

//            merge()
//            setVertexAdjacency(vertexID, putAdjacencyTogether(shortestNeighbourhood, otherNeighbour));

            (*k)--;
        }
        G->setInactive(neighbourOne);
        G->setInactive(neighbourTwo);

        appliedRules->push_back(delVer);

    }
    return APPLICABLE;
}

