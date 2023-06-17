#include "ColorPrint.h"
#include "Reductions.h"
#include "BucketGraph.h"
#include <iostream>
#include <chrono>
#include "boost/intrusive/list.hpp"

using namespace boost::intrusive;

RULE_APPLICATION_RESULT Reductions::rule_HighDegree(BucketGraph* G, int* k)
{
    if(*k < 0) return INSUFFICIENT_BUDGET; //cannot delete more vertices, no possible vertex cover exists
    //std::cout << "HIGHDEG with maxdeg: " << G->getMaxDegree() << " and " << *k << '\n';
    if(!(G->getMaxDegree() > *k) ||*k == 0) return INAPPLICABLE; //cannot apply rule
    
    Reduction* reduction = new Reduction(RULE::HIGH_DEGREE, 0, nullptr, new std::vector<int>());
    appliedRules->push_back(reduction);

    //std::cout << "HIGHDEG Culling vertices: {" << '\n';
    //delete vertices that have to be in the vertex cover
    while(G->getMaxDegree() > *k)
    {
        if(*k <= 0) return INSUFFICIENT_BUDGET; //cannot delete more vertices, no possible vertex cover exists
        int maxDegVertex = G->getMaxDegreeVertex();
        reduction->kDecrement++;
        reduction->deletedVCVertices->push_back(maxDegVertex);
        *k = *k - 1;
        G->setInactive(maxDegVertex);
        //std::cout << maxDegVertex << ", " << '\n';
    }
    //std::cout << "}" << '\n';
    return APPLICABLE;
}

RULE_APPLICATION_RESULT Reductions::rule_Buss(BucketGraph* G, int* k, int numVertices, int numEdges)
{

    /* if(enterDebug){
        std::string enterPrint = "IN BUSS RULE";
        std::cout << ColorPrint::dye(enterPrint, 'p') << '\n' ;
        std::cout << '\n';
    } */

    int k_square = (*k)*(*k);

    if((numVertices > k_square + *k) || (numEdges > k_square))
        return INSUFFICIENT_BUDGET;
    return INAPPLICABLE;
}

RULE_APPLICATION_RESULT Reductions::rule_DegreeOne(BucketGraph* G, int* k, bool checkBudget)
{
    list<BucketVertex>* degOneBucket = G->getVerticesOfDegree(1);

    if(degOneBucket == nullptr || degOneBucket->empty()) return INAPPLICABLE;

    Reduction* reduction = new Reduction(RULE::DEGREE_ONE, 0, nullptr, new std::vector<int>());
    appliedRules->push_back(reduction);

    //std::cout << "DEGONE Culling vertices: {";
    while(!degOneBucket->empty())
    {
        if(*k == 0 && checkBudget) return INSUFFICIENT_BUDGET; //cannot delete more vertices, no possible vertex cover exists
        auto it = degOneBucket->begin();
        int neighbourToDelete = G->getNthActiveNeighbour(it->index, 0);
        reduction->deletedVCVertices->push_back(neighbourToDelete);
        reduction->kDecrement++;
        *k = *k - 1;
        G->setInactive(neighbourToDelete);
        //std::cout << neighbourToDelete << ", ";
    }
    //std::cout << "}" << '\n';
    return APPLICABLE;
}

RULE_APPLICATION_RESULT Reductions::rule_DegreeTwo(BucketGraph* G, int* k, bool checkBudget)
{
    list<BucketVertex>* degTwoBucket = G->getVerticesOfDegree(2);

    if(degTwoBucket == nullptr || degTwoBucket->empty()) return INAPPLICABLE;

    //std::cout << "DEGTWO Culling vertices: " << '\n';
    while(!degTwoBucket->empty())
    {
        //std::cout << "start of while loop" << '\n';
        if(*k == 0 && checkBudget) return INSUFFICIENT_BUDGET; //cannot delete more vertices, no possible vertex cover exists

        auto it = degTwoBucket->begin();

        std::pair<int, int>* neighbours = G->getFirstTwoActiveNeighbours(it->index); //should always return valid vertices
        if(neighbours->first == -1 || neighbours->second == -1)
        {
            G->print();
            throw std::invalid_argument("rule_DegreeTwo: deg2 vertex " + std::to_string(it->index) + " doesn't have two active neighbours\n");
        }


        //std::cout << "degree of deg2 vertex: " << G->getVertexDegree(it->index) << '\n';
        //G->print();
        // save deleted vertex
        Reduction* delVer = new Reduction(RULE::DEGREE_TWO, 0, new std::vector<int>(), new std::vector<int>());

        //if no merge, take neighbours, otherwise case destinction whether merged vertex is in vc
        delVer->deletedVCVertices->push_back(neighbours->first);
        delVer->deletedVCVertices->push_back(neighbours->second);
        delVer->deletedVertices->push_back(it->index);

        //std::cout << it->index  << ", " << neighbours->first << ", " << neighbours->second << '\n';

        delVer->mergeVertexInfo = nullptr;

        // CASE The neighbours know each other, take them into vc
        if(G->vertexHasEdgeTo(neighbours->first, neighbours->second))
        {
            //std::cout << " before nomerge" << '\n';
            delVer->kDecrement = 2;
            //G->print();
            //std::cout << "before setting delVertices inactive: " << it->index << '\n';
            G->setInactive(delVer->deletedVertices);
            //std::cout << "before setting delVCVertices inactive: " << neighbours->first << ", " << neighbours->second << '\n';
            //G->print();
            G->setInactive(delVer->deletedVCVertices);
            //std::cout << " after setting inactive" << '\n';
            //G->print();
            (*k) = (*k) - 2;
        }
        // CASE Neighbours don't know each other => setInactive or do that in addReducedVertices?
        else
        {
            //G->print();
            //G->printBucketQueue();
            //std::cout << " before  merge " << '\n';
            delVer->kDecrement = 1;
            delVer->mergeVertexInfo = G->merge(it->index, neighbours->first, neighbours->second); //sets merged vertices inactive
            (*k) = (*k) - 1;
            //std::cout << " | merging into " << std::get<0>(*delVer->mergeVertexInfo) << '\n';
            //G->print();
            //G->printBucketQueue();
        }
        appliedRules->push_back(delVer);
        //std::cout << "---------" << '\n';
    }
    //std::cout << "----end----" << '\n';
    return APPLICABLE;
}

RULE_APPLICATION_RESULT Reductions::rule_LPFlow(BucketGraph* G, int* k, bool checkBudget)
{
    if(!G->LP_INITIALISED)
        throw std::invalid_argument("LP data structures not initialised");
    std::vector<int>* delVertices = new std::vector<int>();
    std::vector<int>* delVCVertices = new std::vector<int>();
    //auto startHop = std::chrono::high_resolution_clock::now();
    G->hopcroftKarpMatchingSize();
    //auto stopHop = std::chrono::high_resolution_clock::now();
    //auto startCom = std::chrono::high_resolution_clock::now();
    G->setBipartMatchingFlowComponentsInactive(delVertices, delVCVertices, *k, 0.5);
    /* auto stopCom = std::chrono::high_resolution_clock::now();
    double Hop = (std::chrono::duration_cast<std::chrono::microseconds>(stopHop - startHop).count() /  1000) / (double) 1000;
    double Com = (std::chrono::duration_cast<std::chrono::microseconds>(stopCom - startCom).count() /  1000) / (double) 1000; */
    //std::cout << "Computed flow reduction with vc_verts=" << delVCVertices->size() << ", verts=" << delVertices->size() << ", k=" << *k << " in " << Hop + Ed + Com << " seconds (HopcroftKarp: " << Hop << " + EdmondsKarp: " << Ed << " + Connected Components: " << Com << ")" << '\n';
    Reduction* reduction = new Reduction(RULE::LPFLOW, delVCVertices->size(), delVertices, delVCVertices);
    appliedRules->push_back(reduction);
    *k = *k - delVCVertices->size();
    if((int) delVCVertices->size() == 0 && (int) delVertices->size() == 0)
    {
        //std::cout << "LPFlow: INAPPLICABLE" << '\n';
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
    std::cout << "}" << '\n'; */
    if(*k < 0 && checkBudget)
    {
        //std::cout << "LPFlow: INSUFFICIENT_BUDGET" << '\n';
        return INSUFFICIENT_BUDGET;
    }
    return APPLICABLE;
}

RULE_APPLICATION_RESULT Reductions::rule_Domination_BE(BucketGraph* G, int* k, bool checkBudget)
{
    //std::cout << "Domination start" << '\n';
    //G->print();
    int maxDeg = G->getMaxDegree();
    if(maxDeg <= 2) { return INAPPLICABLE; }

    Reduction* reduction = new Reduction(RULE::DEGREE_ONE, 0, nullptr, new std::vector<int>());

    //loop through buckets in descending order
    for(int curDeg = maxDeg; curDeg > 2; curDeg--) //in each while loop, we look at --it, because we started with end()
    {
        //std::cout << "Iterating through bucket " << curDeg << '\n';
        //G->printBucketQueue();
        //list<BucketVertex>* bucket = G->getVerticesOfDegree(curDeg);
        Bucket* bucket = G->getBucket(curDeg);
        if(bucket->vertices.empty())
        {
            continue;
        }

        //loop through vertices of bucket
        auto u_it = bucket->getStableIterator(); // ++bucket->begin(); // TODO: sometimes randomly u_it is an inactive vertex
        while((*u_it)->index != bucket->vertices.end()->index)
        {
            //std::cout << "Checking iterator for vertex u=" << u_it->index << '\n';
            //if no budget left
            if(k == 0 && G->getMaxDegree() > 0 && checkBudget)
            {
                //revert all changes made
                *k = *k + reduction->kDecrement;
                G->setActive(reduction->deletedVCVertices);
                return INSUFFICIENT_BUDGET;
            }
            bool dominates = true;
            //loop through u's neighbours
            Vertex* u = G->getVertex((*u_it)->index);
            //std::cout << "Checking vertex u=" << u->getIndex() << '\n';
            for(int i = 0; i < (int) u->getAdj()->size(); i++)
            {
                Vertex* v = G->getVertex(u->getAdj()->at(i));
                if(!v->getActive()) { continue; }
                if(u->getDegree() < v->getDegree()) { continue; }
                //std::cout << "Checking vertex v=" << v->getIndex() << '\n';

                //loop through v's neighbours
                for(int j = 0; j < (int) v->getAdj()->size(); j++)
                {
                    Vertex* vn = G->getVertex(v->getAdj()->at(j));
                    if(!vn->getActive() || vn->getIndex() == u->getIndex()) { continue; }
                    //std::cout << "Checking vertex v's neighbour vn=" << vn->getIndex() << '\n';
                    //if u does not have edge (u, vn) to neighbour of v
                    if(!G->vertexHasEdgeTo(u->getIndex(), vn->getIndex()))
                    {
                        //std::cout << "Vertex u=" << u->getIndex() << " is missing edge: (" << u->getIndex() << ", " << vn->getIndex() << ")" << '\n';
                        dominates = false;
                        break;
                    }
                    //std::cout << "Vertex u=" << u->getIndex() << " can match edge: (" << v->getIndex() << ", " << vn->getIndex() << ")" << " with edge: (" << u->getIndex() << ", " << vn->getIndex() << ")" << '\n';
                }

                //u dominates neighbour v
                if(dominates)
                {
                    //std::cout << "Vertex u=" << u->getIndex() << " dominates vertex v=" << v->getIndex() << '\n';
                    break;
                }
            }

            //u dominates any neighbour
            if(dominates)
            {
                (*u_it) = ++(*u_it);
                //std::cout << "Vertex u=" << u->getIndex() << " is active " << u->getActive() << " and has degree: " << u->getDegree() << '\n';
                //std::cout << "Setting vertex u=" << u->getIndex() << " inactive " << '\n';
                reduction->kDecrement++;
                *k = *k - 1;
                reduction->deletedVCVertices->push_back(u->getIndex());
                G->setInactive(u->getIndex());
                continue;
            }
            (*u_it) = ++(*u_it);
        }
    }
    if(reduction->deletedVCVertices->size() > 0)
    {
        appliedRules->push_back(reduction);
        return APPLICABLE;
    }
    else
    {
        return INAPPLICABLE;
    }
}

RULE_APPLICATION_RESULT Reductions::rule_Unconfined(BucketGraph* G, int* k, bool checkBudget)
{
    /* list<int> S;
    list<int> neighbours;
    list<int> uniqueNeighbours;
    // TODO: concurrent modification issue, when setting verts inactive?
    // TODO: maybe iterate over vertexreferences
    for (auto vertex = G->getActiveList()->begin(); vertex != G->getActiveList()->end(); ++vertex)
    {
        // How to get neighbourhood of S fast? ---> keep track of it
        // how to check neighbourhood intersection fast? ---> for loop with counter & vertexHasEdgeTo check

        // S and neighbours of S list, that are kept up to date
        S = list<int>();
        neighbours = list<int>();
        uniqueNeighbours = list<int>();
        S.push_back(vertex->getIndex());
        for (auto neighbour = vertex->adj->begin(); neighbour != vertex->adj->end(); ++neighbour)
        {
            neighbours.push_back(*neighbour);
        }
        // search continuation loop
        while(true)
        {
            int best = -1;
            int bestUniqueNeighboursSize = -1;
            uniqueNeighbours.clear();
            // find best neighbour
            for (auto u = neighbours.begin(); u != neighbours.end(); ++u)
            {
                if(!u->getActive()) { continue; }
                bool valid = true;
                int SNeighbours = 0;
                for (auto it=u->adj->begin(); it != u->adj->end(); ++it)
                {
                    if(!u->getActive()) { continue; }
                    if(it->is_linked(neighbours))
                    {
                        uniqueNeighbours.push_back(it->index);
                    }
                    // if u has more than one neighbour in S, u does not fit the criteria
                    if (it->is_linked(S)) { // TODO: check if is_linked works
                        SNeighbours++;
                        if(SNeighbours > 1) { valid = false; break; }
                    }
                }
                if(!valid) { continue; }
                // TODO: update best/bestNeighbourhoodSize if u is better than the current best
                if(best != -1 || uniqueNeighbours.size() < bestUniqueNeighboursSize) { continue; }
                best = u->index;
                bestUniqueNeighboursSize = uniqueNeighbours.size();
            }
            // if no expansion vertex found, vertex "vertex" is not unconfined (continue with next vertex)
            if(best == -1) { break; }
            if(bestUniqueNeighboursSize == 0)
            {
                // TODO: delete vertex and take it into the vc
                break;
            }
            if(bestUniqueNeighboursSize > 0)
            {
                // TODO: push vertex neighbours into S and add non-intersecting neighbourhood to neighbourhood
            }
        }
    } */
    // TODO: if there were any applications return APPLICABLE
    return INAPPLICABLE;
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
    std::cout << ColorPrint::dye(r0, 'r') << '\n' ;
    std::string r1 =  "Rule 1 = " + std::to_string(rule_1);
    std::cout << ColorPrint::dye(r1, 'r') << '\n' ;
    std::string r2 =  "Rule 2 = " + std::to_string(rule_2);
    std::cout << ColorPrint::dye(r2, 'r') << '\n' ;
    std::string rH =  "Rule High Degree = " + std::to_string(rule_High);
    std::cout << ColorPrint::dye(rH, 'r') << '\n' ;
    std::string rD =  "Rule Domination = " + std::to_string(rule_Dom);
    std::cout << ColorPrint::dye(rD, 'r') << '\n' ;
    std::string rL =  "Rule LP-Flow = " + std::to_string(rule_LPF);
    std::cout << ColorPrint::dye(rL, 'r') << '\n' ;
    std::cout << '\n' ;
}

void Reductions::printReductionRules()
{
    if(appliedRules->empty())
        return;

    for (auto i = 0; i < (int) appliedRules->size(); ++i) {
        Reduction* reductionR = appliedRules->at(i);
        std::string toPrint =  "Applied Reduction Rule Nr " + std::to_string(reductionR->rule) + ":";
        std::cout << ColorPrint::dye(toPrint, 'r') << '\n' ;

        std::string deleteV =  "Deleted Vertices are: ";
        for (int j = 0; j < (int) reductionR->deletedVertices->size(); ++j) {
            deleteV += std::to_string(reductionR->deletedVertices->at(j)) + ", ";
        }
        std::cout << ColorPrint::dye(deleteV, 'c') << '\n' ;

//        if(reductionR->rule == DEGREE_TWO)
//        {
//            std::string savedAdj = "Saved Adjacency List: ";
//            for (int k = 0; k < (int) reductionR->savedAdjacency->size(); ++k)
//            {
//                savedAdj += std::to_string(reductionR->savedAdjacency->at(k)) + ", ";
//            }
//            std::cout << ColorPrint::dye(savedAdj, 'c') << '\n';
//        }
//        std::cout << '\n';
    }

}

void Reductions::printDominationSets()
{
//    if(!isThereDomination) {
//        std::cout << "No Domination Possible!" << '\n';
//        return;
//    }

    bool printDebug = true;

    int cnt = 0;

    // Take first out of sorted Dominator
    std::vector<int>* dominationSet;
    int maxSubsetIdx = (int)dominationSets->size();
    std::cout << "There are " << maxSubsetIdx << " dominators" << '\n';

    while(cnt < maxSubsetIdx) {
        dominationSet = dominationSets->at(cnt);

        std::string enterPrint = "Domination Vertex Nr: "; //+ std::to_string(cnt);
        enterPrint += std::to_string(dominationSet->back());
        enterPrint += " dominating " + std::to_string((int) dominationSet->size()-1);
        std::cout << ColorPrint::dye(enterPrint, 'r') << '\n';

        std::string iterPrint;
        for (int i = 0; i < (int) dominationSet->size() - 1; ++i) {
//            std::cout << i  << ", ";
            iterPrint += std::to_string(dominationSet->at(i)) + ", ";
        }
        std::cout << ColorPrint::dye(iterPrint, 'c') << '\n';
        cnt++;
    }
}

bool Reductions::isDominated(BucketGraph* G, int dom, std::vector<bool>* pendingDeletions, bool printDebug)
{
    if (printDebug)
    {
        std::string enterPrint =  "For neighbour: " + std::to_string(dom);
        std::cout << ColorPrint::dye(enterPrint, 'r') << '\n' ;
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
                std::cout << ColorPrint::dye(notNeig, 'c') << '\n' ;
            }
            return false;
        }
    }
    return true;
}

RULE_APPLICATION_RESULT Reductions::rule_Domination(BucketGraph* G, int* k, bool checkBudget)
{
    //std::cout << "----------Domination Rule Start-----------" << '\n';
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
        if(*k - reduction->kDecrement == 0 && checkBudget) {
            (*k) = (*k) + reduction->kDecrement;
            reduction->deletedVCVertices->clear();
            delete reduction->deletedVCVertices;
            delete reduction;
            return INSUFFICIENT_BUDGET; //cannot delete more vertices, no possible vertex cover exists
        }
        list<BucketVertex>* degDom = G->getVerticesOfDegree(maxDeg);

        if(degDom->empty()) {
            maxDeg--;
            continue;
        }

        //std::cout << "Dom: starting iteration for bucket " << maxDeg << '\n';
        // Go through all vertices of degree Bucket
        for(auto it = degDom->begin(); it != degDom->end(); it++)
        {
            if (!G->isActive(it->index) || pendingDeletions[it->index]) { continue; };
            int v = it->index;
            //std::cout << ColorPrint::dye("checking vertex ", 'y') << v << ColorPrint::dye(" for neighbour domination", 'y') << '\n';

            // Get neighbours
            std::vector<int>* neighbours = G->getNeighbours(v);

            // Set flags to neighbours O(deg(v))
            G->dominationHelper->at(v) = 1;
            for (int i = 0; i < (int) neighbours->size(); ++i) {
                G->dominationHelper->at(neighbours->at(i)) = 1;
            }
            //std::cout << ColorPrint::dye("set domination helper", 'y') << '\n';

            // Check if v has one neighbour that is dominated by it
            for (int i = 0; i < (int) neighbours->size(); ++i) {

                // Inactive or cannot be dominated
                if(!G->isActive(neighbours->at(i)) || G->getVertexDegree(neighbours->at(i)) > maxDeg || pendingDeletions[neighbours->at(i)]) {
                    continue;
                }
                //std::cout << ColorPrint::dye("checking neighbour ", 'y') << neighbours->at(i) << ColorPrint::dye(" is dominated", 'y') << '\n';

                if (isDominated(G, neighbours->at(i), &pendingDeletions, false)){
                    for(int j = 0; j < (int) reduction->deletedVCVertices->size(); j++)
                    {
                        if(reduction->deletedVCVertices->at(j) == v)
                        {
                            std::cout << ColorPrint::dye("domination duplicate", 'r') << '\n';
                        }
                    }
                    reduction->kDecrement++;
                    //std::cout << ColorPrint::dye("---> neighbour ", 'g') << neighbours->at(i) << ColorPrint::dye(" is dominated", 'g') << '\n';
                    reduction->deletedVCVertices->push_back(v);
//                    (*k) = (*k) - 1;
                    //G->setInactive(v);
                    pendingDeletions[v] = true;
                    //std::cout << ColorPrint::dye("---> set ", 'g') << v << ColorPrint::dye(" inactive", 'g') << '\n';
                    break;
                }
            }
            
            // Reset flags from neighbours
            G->dominationHelper->at(v) = 0;
            for (int i = 0; i < (int) neighbours->size(); ++i) {
                G->dominationHelper->at(neighbours->at(i)) = 0;
            }
            // free neighbours
            delete neighbours;
            //std::cout << ColorPrint::dye("reset domination helper", 'y') << '\n';

            if(*k - reduction->kDecrement == 0 && checkBudget){
                delete reduction->deletedVCVertices;
                delete reduction;
                return INSUFFICIENT_BUDGET;
            }
        }
        //std::cout << ColorPrint::dye("*******terminated examining bucket********", 'y') << '\n';

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
        //reduction->deletedVCVertices->clear();
        delete reduction->deletedVCVertices;
        delete reduction;
        return INAPPLICABLE;
    }

    G->setInactive(reduction->deletedVCVertices);
    //appliedRules->push_back(reduction);

    if(printDebug){
        std::string cntapp = "DOMINATION RULE APPLIED " + std::to_string(reduction->kDecrement);
        std::cout << ColorPrint::dye(cntapp, 'c') << '\n' ;
//        std::cout << '\n';
    }
    if(printDebug){
        std::string kPrint = "k =" + std::to_string(*k);
        std::cout << ColorPrint::dye(kPrint, 'r') << '\n' ;
        std::cout << '\n';
    }

    if(printTimer)
        if(cntDom < 100)
            std::cout << "This " << cntDom << "-nth reducing of " << reduction->kDecrement << " Vertices took " << duration.count() << "ms." << '\n';

    (*k) = (*k) - reduction->kDecrement;
    appliedRules->push_back(reduction);

    cntDom++;
    return APPLICABLE;
}
