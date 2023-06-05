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

    /* if(enterDebug){
        std::string enterPrint = "IN BUSS RULE";
        std::cout << ColorPrint::dye(enterPrint, 'p') << std::endl ;
        std::cout << std::endl;
    } */

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
    //auto startHop = std::chrono::high_resolution_clock::now();
    G->hopcroftKarpMatchingSize();
    //auto stopHop = std::chrono::high_resolution_clock::now();
    //std::cout << "LPFlow" << std::endl;
    //auto startEd = std::chrono::high_resolution_clock::now();
    G->edmondsKarpFlow();
    //auto stopEd = std::chrono::high_resolution_clock::now();
    //std::cout << "executed edmondsKarp" << std::endl;
    //auto startCom = std::chrono::high_resolution_clock::now();
    G->setBipartMatchingFlowComponentsInactive(delVertices, delVCVertices, *k, 0.5);
    /* G->getBipartMatchingFlowComponents(delVertices, delVCVertices);
    for (int i=0; i < (int) delVertices->size(); i++)
    {
        G->setInactive(delVertices->at(i));
        G->setInactive(delVCVertices->at(i));
    } */
    auto stopCom = std::chrono::high_resolution_clock::now();
    /* double Hop = (std::chrono::duration_cast<std::chrono::microseconds>(stopHop - startHop).count() /  1000) / (double) 1000;
    double Ed = (std::chrono::duration_cast<std::chrono::microseconds>(stopEd - startEd).count() /  1000) / (double) 1000;
    double Com = (std::chrono::duration_cast<std::chrono::microseconds>(stopCom - startCom).count() /  1000) / (double) 1000;
    std::cout << "Computed flow reduction with vc_verts=" << delVCVertices->size() << ", verts=" << delVertices->size() << ", k=" << *k << " in " << Hop + Ed + Com << " seconds (HopcroftKarp: " << Hop << " + EdmondsKarp: " << Ed << " + Connected Components: " << Com << ")" << std::endl; */
    Reduction* reduction = new Reduction(RULE::LPFLOW, delVCVertices->size(), delVertices, delVCVertices);
    appliedRules->push_back(reduction);
    *k = *k - delVCVertices->size();
    if((int) delVCVertices->size() == 0 && (int) delVertices->size() == 0)
    {
        //std::cout << "LPFlow: INAPPLICABLE" << std::endl;
        return INAPPLICABLE;
    }
    /* std::cout << "Deleted component: {";
    for (int j=0; j<(int) delVertices->size(); j++)
    {
        std::cout << delVertices->at(j) << ", ";
    }
    std::cout << "} / ";
    std::cout << "{";
    for (int j=0; j<(int) delVCVertices->size(); j++)
    {
        std::cout << delVCVertices->at(j) << ", ";
    }
    std::cout << "}" << std::endl; */
    if(*k < 0)
    {
        //std::cout << "LPFlow: INSUFFICIENT_BUDGET" << std::endl;
        return INSUFFICIENT_BUDGET;
    }
    return APPLICABLE;
}

RULE_APPLICATION_RESULT Reductions::rule_DegreeTwo(BucketGraph* G, int* k)
{
    list<BucketVertex>* degTwoBucket = G->getVerticesOfDegree(2);

    if(degTwoBucket == nullptr || degTwoBucket->empty()) return INAPPLICABLE;

    //std::cout << "DEGTWO Culling vertices: " << std::endl;
    while(!degTwoBucket->empty())
    {
        //std::cout << "start of while loop" << std::endl;
        if(*k == 0) return INSUFFICIENT_BUDGET; //cannot delete more vertices, no possible vertex cover exists

        auto it = degTwoBucket->begin();

        std::pair<int, int>* neighbours = G->getFirstTwoActiveNeighbours(it->index); //should always return valid vertices
        if(neighbours->first == -1 || neighbours->second == -1)
        {
            G->print();
            throw std::invalid_argument("rule_DegreeTwo: deg2 vertex " + std::to_string(it->index) + " doesn't have two active neighbours\n");
        }


        //std::cout << "degree of deg2 vertex: " << G->getVertexDegree(it->index) << std::endl;
        //G->print();
        // save deleted vertex
        Reduction* delVer = new Reduction(RULE::DEGREE_TWO, 0, new std::vector<int>(), new std::vector<int>());

        //if no merge, take neighbours, otherwise case destinction whether merged vertex is in vc
        delVer->deletedVCVertices->push_back(neighbours->first);
        delVer->deletedVCVertices->push_back(neighbours->second);
        delVer->deletedVertices->push_back(it->index);

        //std::cout << it->index  << ", " << neighbours->first << ", " << neighbours->second << std::endl;

        delVer->mergeVertexInfo = nullptr;

        // CASE The neighbours know each other, take them into vc
        if(G->vertexHasEdgeTo(neighbours->first, neighbours->second))
        {
            //std::cout << " before nomerge" << std::endl;
            delVer->kDecrement = 2;
            //G->print();
            //std::cout << "before setting delVertices inactive: " << it->index << std::endl;
            G->setInactive(delVer->deletedVertices);
            //std::cout << "before setting delVCVertices inactive: " << neighbours->first << ", " << neighbours->second << std::endl;
            //G->print();
            G->setInactive(delVer->deletedVCVertices);
            //std::cout << " after setting inactive" << std::endl;
            //G->print();
            (*k) = (*k) - 2;
        }
        // CASE Neighbours don't know each other => setInactive or do that in addReducedVertices?
        else
        {
            //G->print();
            //G->printBucketQueue();
            //std::cout << " before  merge " << std::endl;
            delVer->kDecrement = 1;
            delVer->mergeVertexInfo = G->merge(it->index, neighbours->first, neighbours->second); //sets merged vertices inactive
            (*k) = (*k) - 1;
            //std::cout << " | merging into " << std::get<0>(*delVer->mergeVertexInfo) << std::endl;
            //G->print();
            //G->printBucketQueue();
        }
        appliedRules->push_back(delVer);
        //std::cout << "---------" << std::endl;
    }
    //std::cout << "----end----" << std::endl;
    return APPLICABLE;
}

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
//    if(!isThereDomination) {
//        std::cout << "No Domination Possible!" << std::endl;
//        return;
//    }

    bool printDebug = true;

    int cnt = 0;

    // Take first out of sorted Dominator
    std::vector<int>* dominationSet;
    int maxSubsetIdx = (int)dominationSets->size();
    std::cout << "There are " << maxSubsetIdx << " dominators" << std::endl;

    while(cnt < maxSubsetIdx) {
        dominationSet = dominationSets->at(cnt);

        std::string enterPrint = "Domination Vertex Nr: "; //+ std::to_string(cnt);
        enterPrint += std::to_string(dominationSet->back());
        enterPrint += " dominating " + std::to_string((int) dominationSet->size()-1);
        std::cout << ColorPrint::dye(enterPrint, 'r') << std::endl;

        std::string iterPrint;
        for (int i = 0; i < (int) dominationSet->size() - 1; ++i) {
//            std::cout << i  << ", ";
            iterPrint += std::to_string(dominationSet->at(i)) + ", ";
        }
        std::cout << ColorPrint::dye(iterPrint, 'c') << std::endl;
        cnt++;
    }
}

bool Reductions::isDominated(BucketGraph* G, int dom, std::vector<bool>* pendingDeletions, bool printDebug)
{
    if (printDebug)
    {
        std::string enterPrint =  "For neighbour: " + std::to_string(dom);
        std::cout << ColorPrint::dye(enterPrint, 'r') << std::endl ;
    }

    if(G->getVertexDegree(dom) == 1)
    {
        return true;
    }
    std::vector<int> * neighbour = G->getNeighbours(dom);

    for (int i = 0; i < (int) neighbour->size(); ++i) {
        int n = neighbour->at(i);
        if(G->isActive(n) && G->dominationHelper->at(n) == 0 && !pendingDeletions->at(n))
//        if(G->isActive(n) && !G->dominationHelperBool->at(n))
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

/* RULE_APPLICATION_RESULT Reductions::rule_Domination(BucketGraph* G, int* k)
{
    int maxDeg = G->getMaxDegree();
    if(maxDeg < 3)
        return INAPPLICABLE;

    // Heuristic
    if(dominationHeuristic && cntDom > 2)
        if(maxDeg > arbitraryDegreeLimiter)
            maxDeg = arbitraryDegreeLimiter;

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration;

    if(printDebug){
        std::string enterPrint = "IN DOMINATION RULE";
        std::cout << ColorPrint::dye(enterPrint, 'p') << std::endl ;
        std::cout << std::endl;
    }

    Reduction* reduction = new Reduction(RULE::DOMINATION, 0, nullptr, new std::vector<int>());

    // While there are Nodes with at least Deg 3
    while(maxDeg > 2) {
        if(*k - reduction->kDecrement == 0) {
            delete reduction;
            return INSUFFICIENT_BUDGET; //cannot delete more vertices, no possible vertex cover exists
        }
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
                    break;
                }
            }

            // Reset flags from neighbours
            G->dominationHelper->at(v) = 0;
            for (int i = 0; i < (int) neighbours->size(); ++i) {
                G->dominationHelper->at(neighbours->at(i)) = 0;
            }

            if(*k - reduction->kDecrement == 0){
                delete reduction;
                return INSUFFICIENT_BUDGET;
            }
        }

        auto stop = std::chrono::high_resolution_clock::now();
        duration = stop - start;

        // Dauer länger als 20ms
        if(duration.count() > 25.0)
        {
            start = std::chrono::high_resolution_clock::now();
            if(maxDeg > 20)
            {
                maxDeg = 20;
                continue;
            }
        }
        maxDeg--;
    }


    rule_Dom += reduction->kDecrement;
    if(reduction->kDecrement == 0){
        delete reduction;
        return INAPPLICABLE;
    }

    if(printDebug){
        std::string cntapp = "DOMINATION RULE APPLIED " + std::to_string(reduction->kDecrement);
        std::cout << ColorPrint::dye(cntapp, 'c') << std::endl ;
//        std::cout << std::endl;
    }
    if(printDebug){
        std::string kPrint = "k =" + std::to_string(*k);
        std::cout << ColorPrint::dye(kPrint, 'r') << std::endl ;
        std::cout << std::endl;
    }

    (*k) = (*k) - reduction->kDecrement;
    G->setInactive(reduction->deletedVCVertices);

    appliedRules->push_back(reduction);

    if(printTimer)
        if(cntDom < 100)
            std::cout << "This " << cntDom << "-nth reducing of " << reduction->kDecrement << " Vertices took " << duration.count() << "ms." << std::endl;

    cntDom++;
    return APPLICABLE;
} */

RULE_APPLICATION_RESULT Reductions::rule_Domination(BucketGraph* G, int* k)
{
    //std::cout << "----------Domination Rule Start-----------" << std::endl;
    int maxDeg = G->getMaxDegree();
    if(maxDeg < 3)
        return INAPPLICABLE;

    // Heuristic
    if(dominationHeuristic && cntDom > 2)
        if(maxDeg > arbitraryDegreeLimiter)
            maxDeg = arbitraryDegreeLimiter;

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration;

    std::vector<bool> pendingDeletions = std::vector<bool>(G->dominationHelper->size());
    for (int i=0; i<(int) pendingDeletions.size(); i++) { pendingDeletions[i] = false; }

    Reduction* reduction = new Reduction(RULE::DOMINATION, 0, nullptr, new std::vector<int>());


    // While there are Nodes with at least Deg 3
    while(maxDeg > 2) {
        if(*k - reduction->kDecrement == 0) {
            (*k) = (*k) + reduction->kDecrement;
            reduction->deletedVCVertices->clear();
            return INSUFFICIENT_BUDGET; //cannot delete more vertices, no possible vertex cover exists
        }
        list<BucketVertex>* degDom = G->getVerticesOfDegree(maxDeg);

        if(degDom->empty()) {
            maxDeg--;
            continue;
        }

        //std::cout << "Dom: starting iteration for bucket " << maxDeg << std::endl;
        // Go through all vertices of degree Bucket
        for(auto it = degDom->begin(); it != degDom->end(); it++)
        {
            if (!G->isActive(it->index) || pendingDeletions[it->index]) { continue; };
            int v = it->index;
            //std::cout << ColorPrint::dye("checking vertex ", 'y') << v << ColorPrint::dye(" for neighbour domination", 'y') << std::endl;

            // Get neighbours
            std::vector<int>* neighbours = G->getNeighbours(v);

            // Set flags to neighbours O(deg(v))
            G->dominationHelper->at(v) = 1;
            for (int i = 0; i < (int) neighbours->size(); ++i) {
                G->dominationHelper->at(neighbours->at(i)) = 1;
            }
            //std::cout << ColorPrint::dye("set domination helper", 'y') << std::endl;

            // Check if v has one neighbour that is dominated by it
            for (int i = 0; i < (int) neighbours->size(); ++i) {

                // Inactive or cannot be dominated
                if(!G->isActive(neighbours->at(i)) || G->getVertexDegree(neighbours->at(i)) > maxDeg || pendingDeletions[neighbours->at(i)]) {
                    continue;
                }
                //std::cout << ColorPrint::dye("checking neighbour ", 'y') << neighbours->at(i) << ColorPrint::dye(" is dominated", 'y') << std::endl;

                if (isDominated(G, neighbours->at(i), &pendingDeletions, false)){
                    for(int j = 0; j < (int) reduction->deletedVCVertices->size(); j++)
                    {
                        if(reduction->deletedVCVertices->at(j) == v)
                        {
                            std::cout << ColorPrint::dye("domination duplicate", 'r') << std::endl;
                        }
                    }
                    reduction->kDecrement++;
                    //std::cout << ColorPrint::dye("---> neighbour ", 'g') << neighbours->at(i) << ColorPrint::dye(" is dominated", 'g') << std::endl;
                    reduction->deletedVCVertices->push_back(v);
//                    (*k) = (*k) - 1;
                    //G->setInactive(v);
                    pendingDeletions[v] = true;
                    //std::cout << ColorPrint::dye("---> set ", 'g') << v << ColorPrint::dye(" inactive", 'g') << std::endl;
                    break;
                }
            }
            
            // Reset flags from neighbours
            G->dominationHelper->at(v) = 0;
            for (int i = 0; i < (int) neighbours->size(); ++i) {
                G->dominationHelper->at(neighbours->at(i)) = 0;
            }
            //std::cout << ColorPrint::dye("reset domination helper", 'y') << std::endl;

            if(*k - reduction->kDecrement == 0){
                delete reduction;
                return INSUFFICIENT_BUDGET;
            }
        }
        //std::cout << ColorPrint::dye("*******terminated examining bucket********", 'y') << std::endl;

        auto stop = std::chrono::high_resolution_clock::now();
        duration = stop - start;

        // Dauer länger als 20ms
        if(duration.count() > 25.0)
        {
            start = std::chrono::high_resolution_clock::now();
            if(maxDeg > 20)
            {
                maxDeg = 20;
                continue;
            }
        }
        maxDeg--;
    }

    if(reduction->kDecrement == 0){
        //reduction->deletedVCVertices = {};
        return INAPPLICABLE;
    }

    G->setInactive(reduction->deletedVCVertices);
    //appliedRules->push_back(reduction);

    if(printDebug){
        std::string cntapp = "DOMINATION RULE APPLIED " + std::to_string(reduction->kDecrement);
        std::cout << ColorPrint::dye(cntapp, 'c') << std::endl ;
//        std::cout << std::endl;
    }
    if(printDebug){
        std::string kPrint = "k =" + std::to_string(*k);
        std::cout << ColorPrint::dye(kPrint, 'r') << std::endl ;
        std::cout << std::endl;
    }

    if(printTimer)
        if(cntDom < 100)
            std::cout << "This " << cntDom << "-nth reducing of " << reduction->kDecrement << " Vertices took " << duration.count() << "ms." << std::endl;

    (*k) = (*k) - reduction->kDecrement;
    appliedRules->push_back(reduction);

    cntDom++;
    return APPLICABLE;
}

/* RULE_APPLICATION_RESULT Reductions::rule_DegreeOne(BucketGraph* G, int* k)
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
} */

//RULE_APPLICATION_RESULT Reductions::rule_Domination_Upgraded(BucketGraph* G, int* k)
//{
//    if(G->deactivatedNodes->empty())
//        return INAPPLICABLE;
//
//    int del = 0;
//
//    Reduction* reduction = new Reduction(RULE::DOMINATION, 0, nullptr, new std::vector<int>());
//
//    // Recheck for Domination Sets
//    for (int i = 0; i < (int) G->deactivatedNodes->size(); ++i) {
//        if(*k == 0) {
//            G->deactivatedNodes->clear();
////            appliedRules->push_back(reduction);
//            return INSUFFICIENT_BUDGET; //cannot delete more vertices, no possible vertex cover exists
//        }
//        std::vector<int>* neighborHoodToCheck = G->getNeighbours(G->deactivatedNodes->at(i));
//        int neighborHoodSize = (int) neighborHoodToCheck->size();
//
//        if(neighborHoodToCheck->empty())
//        {
//            continue;
//        }
//
//        // Check if new neighbourhood has a Dominator
//        for(auto it = neighborHoodToCheck->begin(); it != neighborHoodToCheck->end(); it++)
//        {
//            int v = *it;
//
//            // Get neighbours
//            std::vector<int>* neighbours = G->getNeighbours(v);
//
//            if(neighbours->empty() || neighborHoodSize < 3)
//                continue;
//
//            // Set flags to neighbours O(deg(v))
//            G->dominationHelper->at(v) = 1;
//            for (int i = 0; i < (int) neighbours->size(); ++i) {
//                G->dominationHelper->at(neighbours->at(i)) = 1;
//            }
//
//            // Check if v has one neighbour that is dominated by it
//            for (int i = 0; i < (int) neighbours->size(); ++i) {
//
//                // Inactive or cannot be dominated
//                if(!G->isActive(neighbours->at(i)) || G->getVertexDegree(neighbours->at(i)) > neighborHoodSize) {
//                    continue;
//                }
//
//                if (isDominated(G, neighbours->at(i), false)){
////                if (isDominated2(G, dominationHelper, neighbours->at(i), false)){
//                    reduction->kDecrement++;
//                    reduction->deletedVCVertices->push_back(v);
////                    (*k) = (*k) - 1;
//                    break;
//                }
//            }
//
//            // Reset flags from neighbours
//            G->dominationHelper->at(v) = 0;
//            for (int i = 0; i < (int) neighbours->size(); ++i) {
//                G->dominationHelper->at(neighbours->at(i)) = 0;
//            }
//
//            delete neighbours;
//
//            if(*k == 0){
//                G->deactivatedNodes->clear();
////                appliedRules->push_back(reduction);
//                return INSUFFICIENT_BUDGET;
//            }
//        }
//
////        for (int i = del; i < reduction->kDecrement; ++i) {
////            G->setInactive(reduction->deletedVCVertices->at(i));
////        }
////        del = reduction->kDecrement;
//
//        delete neighborHoodToCheck;
//    }
//
//    G->deactivatedNodes->clear();
//
//    rule_Dom += reduction->kDecrement;
//    if(reduction->kDecrement == 0){
//        delete reduction;
//        return INAPPLICABLE;
//    }
//
//    (*k) = (*k) - reduction->kDecrement;
//    G->setInactive(reduction->deletedVCVertices);
//
//    appliedRules->push_back(reduction);
//
//    return APPLICABLE;
//}

//void Reductions::initDominationVector(BucketGraph *G)
//{
//    int cnt = 0;
//
//    int maxDeg = G->getMaxDegree();
//    if(maxDeg < 3)
//    {
//        isThereDomination = false;
//        return;
//    }
//
//    while(maxDeg > 2) {
//
//        list<BucketVertex>* degDom = G->getVerticesOfDegree(maxDeg);
//
//        if( degDom->empty()) {
//            maxDeg--;
//            continue;
//        }
//
//        // Go through all vertices of degree Bucket
//        for(auto it = degDom->begin(); it != degDom->end(); it++)
//        {
//            std::vector<int>* curSubsets = new std::vector<int>();
//            int v = it->index;
//
//            // Get neighbours
//            std::vector<int>* neighbours = G->getNeighbours(v);
//
//            // Set flags to neighbours O(deg(v))
//            G->dominationHelper->at(v) = 1;
//            for (int i = 0; i < (int) neighbours->size(); ++i) {
//                G->dominationHelper->at(neighbours->at(i)) = 1;
//            }
//
//            // Check if v has one neighbour that is dominated by it
//            for (int i = 0; i < (int) neighbours->size(); ++i) {
//                // Cannot be dominated
//                if(G->getVertexDegree(neighbours->at(i)) > maxDeg) {
//                    continue;
//                }
//                if (isDominated(G, neighbours->at(i), false)){
//                    curSubsets->push_back(neighbours->at(i));
//                    cnt++;
//                }
//            }
//
//            // Reset flags from neighbours
//            G->dominationHelper->at(v) = 0;
//            for (int i = 0; i < (int) neighbours->size(); ++i) {
//                G->dominationHelper->at(neighbours->at(i)) = 0;
//            }
//            // Put the dominant Vertex at the end
//            if(!curSubsets->empty()) {
//                curSubsets->push_back(v);
//                dominationSets->push_back(curSubsets);
//            }
//        }
//
//        maxDeg--;
//    }
//
//    if(cnt == 0)
//        isThereDomination = false;
//    else
//        isThereDomination = true;
//
//    numberOfDominators = cnt;
//}

//RULE_APPLICATION_RESULT Reductions::rule_DominationMitInit(BucketGraph* G, int* k)
//{
//    if(!isThereDomination) return INAPPLICABLE;
//
//    bool printDebug = false;
//    bool enterDebug = printDebug;
//
//    if(enterDebug){
//        std::string enterPrint = "IN DOMINATION RULE";
//        std::cout << ColorPrint::dye(enterPrint, 'y') << std::endl ;
//        std::cout << std::endl;
//    }
//
//    int cnt = 0;
//    int del = 0;
//
//    // Take first out of sorted Dominator in Deg
//    std::vector<int>* dominationSet;
//    int maxSubsetIdx = (int)dominationSets->size();
//    int maxDeg = G->getVertexDegree(dominationSets->at(0)->back());
//    if(maxDeg < 3)
//        return INAPPLICABLE;
//
//    Reduction* reduction = new Reduction(RULE::DOMINATION, 0, nullptr, new std::vector<int>());
//
//    numberOfDominators = 0;
//    for(int l = 0; l < maxSubsetIdx; l++)
//    {
//        // if current Vertex inactive
//        if(G->isActive(dominationSets->at(l)->back()))
//            numberOfDominators ++;
//    }
//
//    if (numberOfDominators == 0){
//        return INAPPLICABLE;
//    }
//    if (printDebug)
//    {
//        std::string domPrint =  "There are " + std::to_string(numberOfDominators) + " active dominators";
//        std::cout << ColorPrint::dye(domPrint, 'c') << std::endl ;
//    }
//
//    for(int l = 0; l < maxSubsetIdx; l++)
//    {
//        if(*k == 0){
//            appliedRules->push_back(reduction);
//            return INSUFFICIENT_BUDGET;
//        }
//
//        // if current Vertex inactive
//        if(!G->isActive(dominationSets->at(l)->back()))
//        {
////            std::cout << "Am I here breaking?" << std::endl ;
////            cnt++;
//            continue;
//        }
//
//        dominationSet = dominationSets->at(l);
//        int dominatorIdx = dominationSet->back();
//
//        if (printDebug)
//        {
//            std::string enterPrint =  "Domination Vertex: " + std::to_string(dominatorIdx);
//            std::cout << ColorPrint::dye(enterPrint, 'r') << std::endl ;
//        }
//
//        std::string iterPrint;
//        for (int i = 0; i < (int) dominationSet->size() - 1; ++i) {
//
//            iterPrint +=  std::to_string(dominationSet->at(i)) + ", ";
//
//            if(G->isActive(dominationSet->at(i)) && G->getVertexDegree(dominationSet->at(i)) > 1)
//            {
//                reduction->kDecrement++;
//                reduction->deletedVCVertices->push_back(dominatorIdx);
//                (*k) = (*k) - 1;
//                G->setInactive(dominatorIdx);
//            }
//            else
//            {
//                continue;
//            }
//            if (printDebug)
//            {
//                std::cout << ColorPrint::dye(iterPrint, 'r') << std::endl ;
//            }
//            break;
//        }
//    }
//
//    rule_Dom += reduction->kDecrement;
//
//    if(reduction->kDecrement == 0)
//    {
//        delete reduction;
//        return INAPPLICABLE;
//    }
//
//    if(printDebug){
//        std::string cntapp = "DOMINATION RULE APPLIED " + std::to_string(reduction->kDecrement);
//        std::cout << ColorPrint::dye(cntapp, 'p') << std::endl ;
//        std::cout << std::endl;
//    }
//
//    appliedRules->push_back(reduction);
//
////    if(rule_Dom == 0)
////        return INAPPLICABLE;
//    return APPLICABLE;
//}
