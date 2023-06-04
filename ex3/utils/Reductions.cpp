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

    if(enterDebug){
        std::string enterPrint = "IN BUSS RULE";
        std::cout << ColorPrint::dye(enterPrint, 'p') << std::endl ;
        std::cout << std::endl;
    }

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

/*----------------------------------------------------------*/
/*------------------   Domination Rule   -------------------*/
/*----------------------------------------------------------*/
void Reductions::initRuleCounter()
{
    rule_0 = 0;
    rule_1 = 0;
    rule_2 = 0;
    rule_High = 0;
    rule_LPF = 0;
    rule_Dom = 0;
}

void Reductions::printCounters()
{
    std::string r0 =  "Rule 0 = " + std::to_string(rule_0);
    std::cout << ColorPrint::dye(r0, 'r') << std::endl ;
    std::string r1 =  "Rule 1 = " + std::to_string(rule_1);
    std::cout << ColorPrint::dye(r1, 'r') << std::endl ;
    std::string r2 =  "Rule 2 = " + std::to_string(rule_2);
    std::cout << ColorPrint::dye(r2, 'r') << std::endl ;
    std::string rH =  "Rule High Degree = " + std::to_string(rule_High);
    std::cout << ColorPrint::dye(rH, 'r') << std::endl ;
    std::string rD =  "Rule Domination = " + std::to_string(rule_Dom);
    std::cout << ColorPrint::dye(rD, 'r') << std::endl ;
    std::string rL =  "Rule LP-Flow = " + std::to_string(rule_LPF);
    std::cout << ColorPrint::dye(rL, 'r') << std::endl ;
    std::cout << std::endl ;
}

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

//        if(reductionR->rule == DEGREE_TWO)
//        {
//            std::string savedAdj = "Saved Adjacency List: ";
//            for (int k = 0; k < (int) reductionR->savedAdjacency->size(); ++k)
//            {
//                savedAdj += std::to_string(reductionR->savedAdjacency->at(k)) + ", ";
//            }
//            std::cout << ColorPrint::dye(savedAdj, 'c') << std::endl;
//        }
//        std::cout << std::endl;
    }

}

void Reductions::printDominationSets()
{
    if(!isThereDomination) return;

    bool printDebug = true;

    int cnt = 0;

    // Take first out of sorted Dominator
    std::vector<int>* dominationSet;
    int maxSubsetIdx = (int)dominationSets->size();
    std::cout << "There are " << maxSubsetIdx << " dominators" << std::endl;

    while(cnt < maxSubsetIdx) {
        dominationSet = dominationSets->at(cnt);

        std::string enterPrint = "Domination Vertex Nr: "; //+ std::to_string(cnt);
        enterPrint += " of Vertex " + std::to_string(dominationSet->back());
        std::cout << ColorPrint::dye(enterPrint, 'r') << std::endl;

        std::cout << "It dominates " << (int) dominationSet->size()-1 << std::endl;
        std::string iterPrint;
        for (int i = 0; i < (int) dominationSet->size() - 1; ++i) {
//            std::cout << i  << ", ";
            iterPrint += std::to_string(dominationSet->at(i)) + ", ";
        }
        std::cout << ColorPrint::dye(iterPrint, 'c') << std::endl;
        cnt++;
    }
}

bool Reductions::isDominated(BucketGraph* G, int dom, bool printDebug)
{
    if (printDebug)
    {
        std::string enterPrint =  "For neighbour: " + std::to_string(dom);
        std::cout << ColorPrint::dye(enterPrint, 'r') << std::endl ;
    }

    std::vector<int> * neighbour = G->getNeighbours(dom);

    if((int) neighbour->size() == 1)
    {
        return false;
    }

    for (int i = 0; i < (int) neighbour->size(); ++i) {
        int n = neighbour->at(i);
        if(G->isActive(n) && G->dominationHelper->at(n) == 0)
        {
            if (printDebug)
            {
                std::string notNeig =  "No subset because of " + std::to_string(n);
                std::cout << ColorPrint::dye(notNeig, 'c') << std::endl ;
            }
            return false;
        }
    }
    return true;
}

RULE_APPLICATION_RESULT Reductions::rule_Domination(BucketGraph* G, int* k)
{
    int del = 0;

    int maxDeg = G->getMaxDegree();
    if(maxDeg < 3)
        return INAPPLICABLE;

    Reduction* reduction = new Reduction(RULE::DOMINATION, 0, nullptr, new std::vector<int>());

    while(maxDeg > 2) {
        if(*k == 0) return INSUFFICIENT_BUDGET; //cannot delete more vertices, no possible vertex cover exists

        list<BucketVertex>* degDom = G->getVerticesOfDegree(maxDeg);

        if( degDom->empty()) {
            maxDeg--;
            continue;
        }

        // Go through all vertices of degree Bucket
        for(auto it = degDom->begin(); it != degDom->end(); it++)
        {
            int v = it->index;

            // Get neighbours
            std::vector<int>* neighbours = G->getNeighbours(v);

            // Set flags to neighbours O(deg(v))
            G->dominationHelper->at(v) = 1;
            for (int i = 0; i < (int) neighbours->size(); ++i) {
                G->dominationHelper->at(neighbours->at(i)) = 1;
            }

            // Check if v has one neighbour that is dominated by it
            for (int i = 0; i < (int) neighbours->size(); ++i) {

                // Inactive or cannot be dominated
                if(!G->isActive(neighbours->at(i)) || G->getVertexDegree(neighbours->at(i)) > maxDeg) {
                    continue;
                }

                if (isDominated(G, neighbours->at(i), false)){
                    reduction->kDecrement++;
                    reduction->deletedVCVertices->push_back(v);
                    (*k) = (*k) - 1;
                    break;
                }

            }

            // Reset flags from neighbours
            G->dominationHelper->at(v) = 0;
            for (int i = 0; i < (int) neighbours->size(); ++i) {
                G->dominationHelper->at(neighbours->at(i)) = 0;
            }

            if(*k == 0){
                for (int i = del; i < reduction->kDecrement; ++i) {
                    G->setInactive(reduction->deletedVCVertices->at(i));
                }
                appliedRules->push_back(reduction);
                return INSUFFICIENT_BUDGET;
            }
        }

        // Currently setting Inactive here because of looping through the buckets
        // because it will corrupt the buckets pointer
        for (int i = del; i < reduction->kDecrement; ++i) {
            G->setInactive(reduction->deletedVCVertices->at(i));
        }

        del = reduction->kDecrement;
        maxDeg--;
    }

    rule_Dom += reduction->kDecrement;
    if(reduction->kDecrement == 0){
        delete reduction;
        return INAPPLICABLE;
    }

    appliedRules->push_back(reduction);

    return APPLICABLE;
}

void Reductions::initDominationVector(BucketGraph *G)
{
    int cnt = 0;

    int maxDeg = G->getMaxDegree();
    if(maxDeg < 3)
    {
        isThereDomination = false;
        return;
    }

    while(maxDeg > 2) {

        list<BucketVertex>* degDom = G->getVerticesOfDegree(maxDeg);

        if( degDom->empty()) {
            maxDeg--;
            continue;
        }

        // Go through all vertices of degree Bucket
        for(auto it = degDom->begin(); it != degDom->end(); it++)
        {
            std::vector<int>* curSubsets = new std::vector<int>();
            int v = it->index;

            // Get neighbours
            std::vector<int>* neighbours = G->getNeighbours(v);

            // Set flags to neighbours O(deg(v))
            G->dominationHelper->at(v) = 1;
            for (int i = 0; i < (int) neighbours->size(); ++i) {
                G->dominationHelper->at(neighbours->at(i)) = 1;
            }

            // Check if v has one neighbour that is dominated by it
            for (int i = 0; i < (int) neighbours->size(); ++i) {
                // Cannot be dominated
                if(G->getVertexDegree(neighbours->at(i)) > maxDeg) {
                    continue;
                }
                if (isDominated(G, neighbours->at(i), false)){
                    curSubsets->push_back(neighbours->at(i));
                    cnt++;
                }
            }

            // Reset flags from neighbours
            G->dominationHelper->at(v) = 0;
            for (int i = 0; i < (int) neighbours->size(); ++i) {
                G->dominationHelper->at(neighbours->at(i)) = 0;
            }
            // Put the dominant Vertex at the end
            if(!curSubsets->empty()) {
                curSubsets->push_back(v);
                dominationSets->push_back(curSubsets);
            }
        }

        maxDeg--;
    }

    if(cnt == 0)
        isThereDomination = false;
    else
        isThereDomination = true;
}

RULE_APPLICATION_RESULT Reductions::rule_DominationMitInit(BucketGraph* G, int* k)
{
    if(!isThereDomination) return INAPPLICABLE;

    bool printDebug = true;
    bool enterDebug = true;

    if(enterDebug){
        std::string enterPrint = "IN DOMINATION RULE";
        std::cout << ColorPrint::dye(enterPrint, 'y') << std::endl ;
        std::cout << std::endl;
    }

    int cnt = 0;
    int del = 0;

    // Take first out of sorted Dominator
    std::vector<int>* dominationSet;
    int maxSubsetIdx = (int)dominationSets->size() -1;
    int maxDeg = G->getVertexDegree(dominationSets->at(0)->back());
    if(maxDeg < 3)
        return INAPPLICABLE;

    Reduction* reduction = new Reduction(RULE::DOMINATION, 0, nullptr, new std::vector<int>());

    if (printDebug)
    {
        std::string domPrint =  "There are " + std::to_string(maxSubsetIdx) + " dominators";
        std::cout << ColorPrint::dye(domPrint, 'c') << std::endl ;
    }

    while(cnt < maxSubsetIdx)
    {
        if(*k == 0) return INSUFFICIENT_BUDGET; //cannot delete more vertices, no possible vertex cover exists

        // if current Vertex inactive
        if(!G->isActive(dominationSets->at(cnt)->back()))
        {
//            std::cout << "Am I here breaking?" << std::endl ;
            cnt++;
            continue;
        }

        dominationSet = dominationSets->at(cnt);
        int dominatorIdx = dominationSet->back();

        if (printDebug)
        {
            std::string enterPrint =  "Domination Nr: " + std::to_string(cnt+1);
            enterPrint += " of Vertex " + std::to_string(dominatorIdx);
            std::cout << ColorPrint::dye(enterPrint, 'r') << std::endl ;
        }


        std::string iterPrint;
        for (int i = 0; i < (int) dominationSet->size() - 1; ++i) {

            iterPrint +=  std::to_string(dominationSet->at(i)) + ", ";


            if(G->isActive(dominationSet->at(i)) && G->getVertexDegree(dominationSet->at(i)) > 0)
            {
//                std::cout << "Am I here?" << std::endl ;
                reduction->kDecrement++;
                reduction->deletedVCVertices->push_back(dominatorIdx);
                (*k) = (*k) - 1;
                G->setInactive(dominatorIdx);
//                G->printBucketQueue();
//                G->print();
            }
            else
            {
//                std::cout << "Or There ?" << std::endl ;
                continue;
            }
            if (printDebug)
            {
                std::cout << ColorPrint::dye(iterPrint, 'r') << std::endl ;
            }
            break;

        }
//        std::cout << "Am I here though?" << std::endl ;
        cnt++;
    }

    rule_Dom += reduction->kDecrement;
    appliedRules->push_back(reduction);

    if(rule_Dom == 0)
        return INAPPLICABLE;
    return APPLICABLE;
}