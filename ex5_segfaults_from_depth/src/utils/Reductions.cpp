#include "ColorPrint.h"
#include "Reductions.h"
#include "BucketGraph.h"
#include <iostream>
#include <chrono>
#include <random>
#include "boost/intrusive/list.hpp"

using namespace boost::intrusive;
typedef ColorPrint cp;

RULE_APPLICATION_RESULT Reductions::rule_HighDegree(BucketGraph* G, int* k, int depth)
{
    if(*k < 0) return INSUFFICIENT_BUDGET; //cannot delete more vertices, no possible vertex cover exists
    //std::cout << "HIGHDEG with maxdeg: " << G->getMaxDegree() << " and " << *k << '\n';
    if(!(G->getMaxDegree() > *k) ||*k == 0) return INAPPLICABLE; //cannot apply rule
    
    Reduction* reduction = new Reduction(RULE::HIGH_DEGREE, 0, nullptr, new std::vector<int>());
    reduction->rDepth = depth;
    appliedRules->push_back(reduction);

    //std::cout << "HIGHDEG Culling vertices: {" << '\n';
    //delete vertices that have to be in the vertex cover
    while(G->getMaxDegree() > *k)
    {
        if(*k <= 0) return INSUFFICIENT_BUDGET; //cannot delete more vertices, no possible vertex cover exists
        int maxDegVertex = G->getMaxDegreeVertex();
        reduction->kDecrement++;
        reduction->deletedVCVertices->push_back(maxDegVertex);

//        std::cout << "REDUCE: High Degree" << std::endl;
//        std::cout << "Vertex: " << maxDegVertex << '\n';
//        std::cout << "-----------" << '\n';

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

RULE_APPLICATION_RESULT Reductions::rule_DegreeOne(BucketGraph* G, int* k, int depth, bool checkBudget, bool printDebug)
{
    list<BucketVertex>* degOneBucket = G->getVerticesOfDegree(1);

    if(degOneBucket == nullptr || degOneBucket->empty()) return INAPPLICABLE;

    Reduction* reduction = new Reduction(RULE::DEGREE_ONE, 0, nullptr, new std::vector<int>());
    reduction->rDepth = depth;

    //std::cout << "DEGONE Culling vertices: {";
    auto startDeg1 = std::chrono::high_resolution_clock::now();
    while(!degOneBucket->empty())
    {
        //revert if no budget
        if(*k == 0 && checkBudget)
        {
            G->setActive(reduction->deletedVCVertices);
            *k = *k + reduction->kDecrement;
            freeReductionRule(reduction, false);
            return INSUFFICIENT_BUDGET;
        } //cannot delete more vertices, no possible vertex cover exists



        auto it = degOneBucket->begin();
        int neighbourToDelete = G->getNthActiveNeighbour(it->index, 0);
        reduction->deletedVCVertices->push_back(neighbourToDelete);
        reduction->kDecrement++;
        *k = *k - 1;
        G->setInactive(neighbourToDelete);

//        std::cout << "REDUCE: Deg 1" << std::endl;
//        std::cout << "Degree 1 vertex: " << it->index << '\n';
//        std::cout << neighbourToDelete << '\n';
//        std::cout << "-----------" << '\n';
        //std::cout << neighbourToDelete << ", ";
    }
    auto stopDeg1 = std::chrono::high_resolution_clock::now();
    double Deg1 = (std::chrono::duration_cast<std::chrono::microseconds>(stopDeg1 - startDeg1).count() /  1000) / (double) 1000;
    if (printDebug)
        std::cout << "#Reduced " << reduction->deletedVCVertices->size() << " Deg1 neighbour vertices in " << Deg1 << " seconds" << std::endl;
    appliedRules->push_back(reduction);
    return APPLICABLE;
}

RULE_APPLICATION_RESULT Reductions::rule_DegreeTwo(BucketGraph* G, int* k, int depth, bool checkBudget, bool printDebug)
{
    list<BucketVertex>* degTwoBucket = G->getVerticesOfDegree(2);

    if(degTwoBucket == nullptr || degTwoBucket->empty()) return INAPPLICABLE;

//    std::cout << "DEGTWO Culling vertices: " << '\n';
    int numberOfReducedVertices = 0;
    int numberOfReducedVCVertices = 0;
    auto startDeg2 = std::chrono::high_resolution_clock::now();
    while(!degTwoBucket->empty())
    {
//        std::cout << "start of while loop" << '\n';
        if(*k == 0 && checkBudget) return INSUFFICIENT_BUDGET; //cannot delete more vertices, no possible vertex cover exists

        auto it = degTwoBucket->begin();

        std::pair<int, int> neighbours = G->getFirstTwoActiveNeighbours(it->index); //should always return valid vertices

//        std::cout << "REDUCE: Deg 2" << std::endl;
//        std::cout << "Degree 2 vertex: " << it->index << '\n';
//        std::cout << neighbours.first << ", " << neighbours.second << '\n';
//        std::cout << "-----------" << '\n';

        if(neighbours.first == -1 || neighbours.second == -1)
        {
            G->print();
            throw std::invalid_argument("rule_DegreeTwo: deg2 vertex " + std::to_string(it->index) + " doesn't have two active neighbours\n");
        }


//        std::cout << "degree of deg2 vertex: " << G->getVertexDegree(it->index) << '\n';
        //G->print();
        // save deleted vertex
        Reduction* delVer = new Reduction(RULE::DEGREE_TWO, 0, new std::vector<int>(), new std::vector<int>());
        delVer->rDepth = depth;

        //if no merge, take neighbours, otherwise case destinction whether merged vertex is in vc
        delVer->deletedVCVertices->push_back(neighbours.first);
        delVer->deletedVCVertices->push_back(neighbours.second);
        delVer->deletedVertices->push_back(it->index);

//        std::cout << it->index  << ", " << neighbours.first << ", " << neighbours.second << '\n';

        delVer->mergeVertexInfo = nullptr;

        // CASE The neighbours know each other, take them into vc
        if(G->vertexHasEdgeTo(neighbours.first, neighbours.second))
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
            numberOfReducedVertices++;
            numberOfReducedVCVertices += 2;
        }
        // CASE Neighbours don't know each other => setInactive or do that in addReducedVertices?
        else
        {
            //G->print();
            //G->printBucketQueue();
//            std::cout << " before  merge " << '\n';
            delVer->kDecrement = 1;
            delVer->mergeVertexInfo = G->merge(it->index, neighbours.first, neighbours.second); //sets merged vertices inactive
            (*k) = (*k) - 1;
            numberOfReducedVertices++;
            numberOfReducedVCVertices++;
//            std::cout << " | merging into " << std::get<0>(*delVer->mergeVertexInfo) << '\n';
            //G->print();
            //G->printBucketQueue();
            /* if(std::get<0>(*delVer->mergeVertexInfo) == 33) {
                G->print();
                std::cout << "reducing mergeVertex " << 33 << std::endl;
            } */
        }
        appliedRules->push_back(delVer);
        //std::cout << "---------" << '\n';
    }
    auto stopDeg2 = std::chrono::high_resolution_clock::now();
    double Deg2 = (std::chrono::duration_cast<std::chrono::microseconds>(stopDeg2 - startDeg2).count() /  1000) / (double) 1000;
    if (printDebug)
        std::cout << "#Reduced " << numberOfReducedVertices << " Deg2 vertices in " << Deg2 << " seconds" << std::endl;
//    std::cout << "----end----" << '\n';
    return APPLICABLE;
}

RULE_APPLICATION_RESULT Reductions::rule_DegreeTwo_Simple_Case(BucketGraph* G, int* k, int depth, bool checkBudget, bool printDebug)
{
    list<BucketVertex>* degTwoBucket = G->getVerticesOfDegree(2);

    if(degTwoBucket == nullptr || degTwoBucket->empty()) return INAPPLICABLE;

    //std::cout << "DEGTWO Culling vertices: " << '\n';
    int numberOfReducedVertices = 0;
    int numberOfReducedVCVertices = 0;
    auto startDeg2 = std::chrono::high_resolution_clock::now();
//    while(!degTwoBucket->empty())

    int cnt = 0;
    std::vector<int> tempDeleted;
    std::unordered_map<int, int> alreadyInactive;

    for(auto it = degTwoBucket->begin(); it != degTwoBucket->end(); it++) {
        //std::cout << "start of while loop" << '\n';
        if (*k == 0 && checkBudget)
            return INSUFFICIENT_BUDGET; //cannot delete more vertices, no possible vertex cover exists

//        auto it = degTwoBucket->begin();

        std::pair<int, int> neighbours = G->getFirstTwoActiveNeighbours(
                it->index); //should always return valid vertices
        if (neighbours.first == -1 || neighbours.second == -1) {
            G->print();
            throw std::invalid_argument("rule_DegreeTwo: deg2 vertex " + std::to_string(it->index) +
                                        " doesn't have two active neighbours\n");
        }
        if(alreadyInactive[neighbours.first] == 1 || alreadyInactive[neighbours.second] == 1)
            continue;

        if (!G->vertexHasEdgeTo(neighbours.first, neighbours.second)) {
            continue;
        }

        tempDeleted.push_back(it->index);
        tempDeleted.push_back(neighbours.first);
        tempDeleted.push_back(neighbours.second);

        alreadyInactive[it->index] = 1;
        alreadyInactive[neighbours.first] = 1;
        alreadyInactive[neighbours.second] = 1;

        cnt++;
    }

    if(tempDeleted.empty())
        return INAPPLICABLE;

    for (int i = 0; i < cnt; ++i) {
        Reduction* delVer = new Reduction(RULE::DEGREE_TWO, 0, new std::vector<int>(), new std::vector<int>());
        delVer->rDepth = depth;

        int v = tempDeleted.at(3*i+0);
        int n1 = tempDeleted.at(3*i+1);
        int n2 = tempDeleted.at(3*i+2);
        //if no merge, take neighbours, otherwise case destinction whether merged vertex is in vc
        delVer->deletedVCVertices->push_back(n1);
        delVer->deletedVCVertices->push_back(n2);
        delVer->deletedVertices->push_back(v);

        delVer->mergeVertexInfo = nullptr;

        delVer->kDecrement = 2;

        G->setInactive(v);
        G->setInactive(n1);
        G->setInactive(n2);

        (*k) = (*k) - 2;
        numberOfReducedVertices++;
        numberOfReducedVCVertices += 2;
        appliedRules->push_back(delVer);
        //std::cout << "---------" << '\n';
    }
    auto stopDeg2 = std::chrono::high_resolution_clock::now();
    double Deg2 = (std::chrono::duration_cast<std::chrono::microseconds>(stopDeg2 - startDeg2).count() /  1000) / (double) 1000;
    if (true || printDebug)
        std::cout << "#Reduced " << numberOfReducedVertices << " Deg2 vertices in " << Deg2 << " seconds" << std::endl;

    //std::cout << "----end----" << '\n';
    return APPLICABLE;
}

RULE_APPLICATION_RESULT Reductions::rule_DegreeTwo_Secure(BucketGraph* G, int* k)
{
    list<BucketVertex>* degTwoBucket = G->getVerticesOfDegree(2);

    if(degTwoBucket == nullptr || degTwoBucket->empty()) return INAPPLICABLE;

    int cnt = 0;
    std::vector<int> tempVC;
    std::vector<int> tempDeleted;
    std::unordered_map<int, int> alreadyInactive;
//
    for(auto it = degTwoBucket->begin(); it != degTwoBucket->end(); it++)
    {
//        std::cout << "start of while loop" << '\n';
        if(!G->isActive(it->index))
            continue;
//        std::cout << "Degree 2 vertex: " << it->index << '\n';

//        std::pair<int, int> neighbours = G->getActiveTwoNeighbours(it->index);
        std::pair<int, int> neighbours = G->getFirstTwoActiveNeighbours(it->index);

        if(alreadyInactive[it->index] == 1 || alreadyInactive[neighbours.first] == 1 || alreadyInactive[neighbours.second] == 1)
        {
//            std::cout << "Nope not this vertex: " << it->index << '\n';
//            std::cout << neighbours.first << ", " << neighbours.second << '\n';
            continue;
        }
//
//        std::cout << "degree of deg2 vertex: " << it->index << '\n';
        if(!G->isActive(neighbours.first) || !G->isActive(neighbours.second)) {
            continue;
        }

//        std::cout << neighbours.first << ", " << neighbours.second << '\n';
//
//
//        // CASE The neighbours know each other, take them into vc
        if(G->vertexHasEdgeTo(neighbours.first, neighbours.second))
        {
//            std::cout << "Neighbours know each other! " << '\n';

            int n1 = neighbours.first;
            int n2 = neighbours.second;
            tempDeleted.push_back(n1);
            tempDeleted.push_back(n2);
            tempDeleted.push_back(it->index);

            cnt++;

            alreadyInactive[it->index] = 1;
            alreadyInactive[neighbours.first] = 1;
            alreadyInactive[neighbours.second] = 1;

        }
//        std::cout << "---------" << '\n';
    }

    if(tempDeleted.empty() || cnt == 0) {
        return INAPPLICABLE;
    }
//    std::cout << "I have " << cnt<< " elements\n";

    for (int i = 0; i < cnt; ++i) {
//        std::cout << "I deleted Vertex: " << tempDeleted.at(3*i) << tempDeleted.at(3*i+1) << tempDeleted.at(3*i+2) << '\n';
        Reduction* delVer = new Reduction(RULE::DEGREE_TWO, 0, new std::vector<int>(), new std::vector<int>());
        delVer->deletedVCVertices->push_back(tempDeleted.at(3*i + 0));
        delVer->deletedVCVertices->push_back(tempDeleted.at(3*i + 1));
        delVer->deletedVertices->push_back(tempDeleted.at(3*i + 2));

        G->setInactive(delVer->deletedVCVertices);
        G->setInactive(delVer->deletedVertices);

        appliedRules->push_back(delVer);
    }

    return APPLICABLE;
//    return INAPPLICABLE;
}

bool Reductions::isTouchable(BucketGraph* G, std::vector<int> * neighbours, std::unordered_map<int, int>& inactiveMap)
{
    bool touchable = true;
    for (int i = 0; i < (int)neighbours->size(); ++i) {
        if(G->isActive(neighbours->at(i)) && inactiveMap[neighbours->at(i)] == 1) {
            touchable = false;
            break;
        }
    }
    return touchable;
}

RULE_APPLICATION_RESULT Reductions::rule_DegreeThree_Independent(BucketGraph* G, int depth, bool printDebug)
{
    list<BucketVertex>* degThreeBucket = G->getVerticesOfDegree(3);



    if(degThreeBucket == nullptr || degThreeBucket->empty()) {
        return INAPPLICABLE; }

    if(G->deg3ind) {
        std::cout << "\nRule: Degree 3: Indenpendent Set" << '\n';
//        std::cout << "Deg 3 Bucket of size: " << (int)degThreeBucket->size() << '\n';
    }

    auto startDeg3Ind = std::chrono::high_resolution_clock::now();

    int cnt = 0;

    std::vector<int> tempDeleted;
    std::vector<int> tempNeighbours;
    std::unordered_map<int, int> alreadyInactive;
    std::vector<std::vector<int>> tempAddedEdges;

    for(auto it = degThreeBucket->begin(); it != degThreeBucket->end(); it++)
    {
        int v = it->index;

        if(!G->isActive(v)) {
            continue;
        }

        std::vector<int>* neighbours = G->getNeighbours(v);

        int a = neighbours->at(0);
        int b = neighbours->at(1);
        int c = neighbours->at(2);

        if(!G->isActive(a) || !G->isActive(b) || !G->isActive(c)) {
            delete neighbours;
            continue;
        }

        if(G->vertexHasEdgeTo(a, b) || G->vertexHasEdgeTo(a, c) || G->vertexHasEdgeTo(b, c)) {
            if(G->deg3ind) {
                std::cout << "Neighbours are not independent!\n";
                std::cout << a << ", " << b << ", " << c << '\n';
            }
            delete neighbours;
            continue;
        }

        if(alreadyInactive[v] == 1 || alreadyInactive[a] == 1 || alreadyInactive[b] == 1 || alreadyInactive[c] == 1)
        {
            if(G->deg3ind) {
                std::cout << "Already treated: " << v << ", maybe next round!\n";
                std::cout << a << ", " << b << ", " << c << '\n';
            }
            delete neighbours;
            continue;
        }

        std::vector<int> addedEdgesToA;
        std::vector<int> addedEdgesToB;
        std::vector<int> addedEdgesToC;

        // Edges {a, b} {b, c}
        addedEdgesToA.push_back(b);
        addedEdgesToB.push_back(c);

        std::vector<int>* nA = G->getNeighbours(a);
        std::vector<int>* nB = G->getNeighbours(b);
        std::vector<int>* nC = G->getNeighbours(c);

        if(!isTouchable(G, nA, alreadyInactive)
           || !isTouchable(G, nB, alreadyInactive)
           || !isTouchable(G, nC, alreadyInactive)) {
            delete nA;
            delete nB;
            delete nC;
            delete neighbours;
            continue;
        }

        // Edges {a, N(b)}
        for (int i = 0; i < (int)nB->size(); ++i) {
//            if(G->isActive(nB->at(i)) && !G->vertexHasEdgeTo(a, nB->at(i)) && alreadyInactive[nB->at(i)] != 1){
            if(G->isActive(nB->at(i)) && G->getVertexDegree(nB->at(i)) != 0 && !G->vertexHasEdgeTo(a, nB->at(i))){
                addedEdgesToA.push_back(nB->at(i));
            }
        }

        // Edges {b, N(c)}
        for (int i = 0; i < (int)nC->size(); ++i) {
            if(G->isActive(nC->at(i)) && G->getVertexDegree(nC->at(i)) != 0 && !G->vertexHasEdgeTo(b, nC->at(i))){
//            if(G->isActive(nC->at(i)) && !G->vertexHasEdgeTo(b, nC->at(i)) && alreadyInactive[nC->at(i)] != 1) {
                addedEdgesToB.push_back(nC->at(i));
            }
        }

        // Edges {c, N(a)}
        for (int i = 0; i < (int)nA->size(); ++i) {
            if(G->isActive(nA->at(i)) && G->getVertexDegree(nA->at(i)) != 0 && !G->vertexHasEdgeTo(c, nA->at(i))){
//            if(G->isActive(nA->at(i)) && !G->vertexHasEdgeTo(c, nA->at(i)) && alreadyInactive[nA->at(i)] != 1) {
                addedEdgesToC.push_back(nA->at(i));
            }
        }

        alreadyInactive[it->index] = 1;
        alreadyInactive[a] = 1;
        alreadyInactive[b] = 1;
        alreadyInactive[c] = 1;

        tempDeleted.push_back(v);
        tempNeighbours.push_back(a);
        tempNeighbours.push_back(b);
        tempNeighbours.push_back(c);

        tempAddedEdges.push_back(addedEdgesToA);
        tempAddedEdges.push_back(addedEdgesToB);
        tempAddedEdges.push_back(addedEdgesToC);

        cnt++;

        delete neighbours;
        delete nA;
        delete nB;
        delete nC;

        if(G->deg3ind) {
            std::cout << "VALID = " << cnt+1 << '\n';
            std::string foundV = "Found vertex: " + std::to_string(v) + '\n';
            std::cout << ColorPrint::dye(foundV, 'y');
            std::cout << a << ", " << b << ", " << c << '\n';
            std::cout << "INDEPENDENT SET!" << '\n';
            std::cout << ColorPrint::dye(lineDelimiter, 'y');
        }

    }

    if(tempDeleted.empty()) {
        if(G->deg3ind)
            std::cout << "Nothing was applied\n";
        return INAPPLICABLE;
    }

    int edgeCnt = 0;

    for (int i = 0; i < cnt; ++i) {
        Reduction* delVer = new Reduction(RULE::DEGREE_THREE_IND, 0, new std::vector<int>(), new std::vector<int>());
        delVer->addedEdges = new std::vector<std::vector<int>>;

        G->cntDeg3Ind++;

        delVer->rDepth = depth;
        delVer->deletedVertices->push_back(tempDeleted.at(i));

        int a = tempNeighbours.at(3*i + 0);
        int b = tempNeighbours.at(3*i + 1);
        int c = tempNeighbours.at(3*i + 2);

        std::vector<int> edgeToA = tempAddedEdges.at(3*i+0);
        std::vector<int> edgeToB = tempAddedEdges.at(3*i+1);
        std::vector<int> edgeToC = tempAddedEdges.at(3*i+2);

        for (int j = 0; j < (int)edgeToA.size(); ++j) {
            if(G->getVertexDegree(edgeToA.at(j)) != 0)
                G->addEdgeToVertex(a, edgeToA.at(j));
        }

        for (int j = 0; j < (int)edgeToB.size(); ++j) {
            if(G->getVertexDegree(edgeToB.at(j)) != 0)
                G->addEdgeToVertex(b, edgeToB.at(j));
        }

        for (int j = 0; j < (int)edgeToC.size(); ++j) {
            if(G->getVertexDegree(edgeToC.at(j)) != 0)
                G->addEdgeToVertex(c, edgeToC.at(j));
        }

        edgeCnt += (int)edgeToA.size();
        edgeCnt += (int)edgeToB.size();
        edgeCnt += (int)edgeToC.size();

        delVer->deletedVCVertices->push_back(a);
        delVer->deletedVCVertices->push_back(b);
        delVer->deletedVCVertices->push_back(c);
        delVer->addedEdges->push_back(edgeToA);
        delVer->addedEdges->push_back(edgeToB);
        delVer->addedEdges->push_back(edgeToC);

        G->setInactive(delVer->deletedVertices);
        appliedRules->push_back(delVer);

        if(G->deg3ind) {
            std::cout << "\nApplying rule to v = " << tempDeleted.at(i) << '\n';
            std::cout << "At recursion = " << delVer->rDepth << '\n';

            std::cout << "Edge 1:" << a << " with:" << '\n';
            for (int j = 0; j < (int)edgeToA.size(); ++j) {
                std::cout << edgeToA.at(j) << '\n';
            }
            std::cout << "Edge 2:" << b << " with:" << '\n';
            for (int j = 0; j < (int)edgeToB.size(); ++j) {
                std::cout << edgeToB.at(j) << '\n';
            }
            std::cout << "Edge 3:" << c << " with:" << '\n';
            for (int j = 0; j < (int)edgeToC.size(); ++j) {
                std::cout << edgeToC.at(j) << '\n';
            }
        }

    }
    auto stopDeg3Ind = std::chrono::high_resolution_clock::now();
    double Deg3Ind = (std::chrono::duration_cast<std::chrono::microseconds>(stopDeg3Ind - startDeg3Ind).count() /  1000) / (double) 1000;

    if (printDebug) {
        std::string deg3indred = "Reduced "+ std::to_string(cnt) + " Deg3: Independent Set: in ";
        deg3indred += std::to_string(Deg3Ind)+ " seconds adding " + std::to_string(edgeCnt);
        deg3indred += " edges at rec Depth " + std::to_string(depth) + '\n';
//        std::cout << '#' <<ColorPrint::dye(deg3indred, 'p');
        std::cout << '#' << deg3indred;
    }

    return APPLICABLE;
}

RULE_APPLICATION_RESULT Reductions::rule_DegreeThree_Clique(BucketGraph* G, int* k, int depth, bool checkBudget, bool printDebug)
{
    list<BucketVertex>* degThreeBucket = G->getVerticesOfDegree(3);

    if(degThreeBucket == nullptr || degThreeBucket->empty()) {
//        if(G->deg3clique)
//            std::cout << "\nNo degree 3 Vertex!" << '\n';
        return INAPPLICABLE; }

    if(G->deg3clique) {
        std::cout << "\nRule: Degree 3: 2-Clique-Neighborhood" << '\n';
//        std::cout << "Deg 3 Bucket of size: " << (int)degThreeBucket->size() << '\n';
        std::cout << "---------" << '\n';
    }

    auto startDeg3clique = std::chrono::high_resolution_clock::now();

    int cnt = 0;

    std::vector<int> tempDeleted;
    std::vector<int> tempNeighbours;
    std::unordered_map<int, int> alreadyInactive;
    std::vector<std::vector<int>> tempAddedEdges;

    for(auto it = degThreeBucket->begin(); it != degThreeBucket->end(); it++)
    {
        if(*k - cnt < 1 && checkBudget)
            return INSUFFICIENT_BUDGET;

        int v = it->index;

        if(!G->isActive(v)) {
            continue;
        }

        std::vector<int>* neighbours = G->getNeighbours(v);

        int a = neighbours->at(0);
        int b = neighbours->at(1);
        int c = neighbours->at(2);

        if(!G->isActive(a) || !G->isActive(b) || !G->isActive(c)) {
            delete neighbours;
            continue;
        }

        if(alreadyInactive[v] == 1 || alreadyInactive[a] == 1 || alreadyInactive[b] == 1 || alreadyInactive[c] == 1)
        {
//            if(G->deg3clique) {
//                std::cout << "Already treated: " << v << ", maybe next round!\n";
//                std::cout << a << ", " << b << ", " << c << '\n';
//                std::cout << "---------" << '\n';
//            }
            delete neighbours;
            continue;
        }

        int c11;
        int c12;
        int c2;

        bool connectionAB = G->vertexHasEdgeTo(a, b);
        bool connectionAC = G->vertexHasEdgeTo(a, c);
        bool connectionBC = G->vertexHasEdgeTo(b, c);

        // C1: a-b, C2: c
        if (connectionAB && !connectionAC && !connectionBC)
        {
            c11 = a;
            c12 = b;
            c2 = c;
        }// C1: a-c, C2: b
        else if (!connectionAB && connectionAC && !connectionBC)
        {
            c11 = a;
            c12 = c;
            c2 = b;
        }// C1: b-c, C2: a
        else if (!connectionAB && !connectionAC && connectionBC)
        {
            c11 = b;
            c12 = c;
            c2 = a;
        }
        else{
//            if(!connectionAB && !connectionAC && !connectionBC) {
//                if(G->deg3clique) {
//                    std::cout << "Neighbours are independent!\n";
//                    std::cout << a << ", " << b << ", " << c << '\n';
//                    std::cout << "---------" << '\n';
//                }
//            }
//            else{
//                if(G->deg3clique) {
//                    std::cout << "Neighbours dominate: " << v << '\n';
//                    std::cout << a << ", " << b << ", " << c << '\n';
//                    std::cout << "---------" << '\n';
//                }
//            }
            delete neighbours;
            continue;
        }

//        std::cout << "REDUCE: Deg 3 Clique" << std::endl;
//        std::cout << "Degree 3 vertex: " << v << '\n';
//        std::cout << c11 << ", " << c12 << ", " << c2 << '\n';
//        std::cout << "-----------" << '\n';

        if(G->deg3clique) {
            std::cout << ColorPrint::dye(lineDelimiter, 'y');
            std::cout << "VALID = " << cnt+1 << '\n';
            std::string foundV = "Found vertex: " + std::to_string(v) + '\n';
            std::cout << "2-CLIQUE!" << '\n';
            std::cout << ColorPrint::dye(foundV, 'y');
            std::cout << c11 << ", " << c12 << ", " << c2 << '\n';
            std::cout << ColorPrint::dye(lineDelimiter, 'y');
        }

        std::vector<int>* nC2 = G->getNeighbours(c2);

        std::vector<int> addedEdgesToC11;
        std::vector<int> addedEdgesToC12;

        if(!isTouchable(G, nC2, alreadyInactive))
        {
            delete nC2;
            delete neighbours;
            continue;
        }

//        std::cout << "C2: " << c2 << " with deg(c2): " << G->getVertexDegree(c2) << '\n';
        // Adding edges to N(C2)
        for (int i = 0; i < (int)nC2->size(); ++i) {

            if(G->isActive(nC2->at(i))){ // && G->getVertexDegree(nC2->at(i)) > 0){
                // Edges {C11, N(C2)}
                if(!G->vertexHasEdgeTo(c11, nC2->at(i)))
                    addedEdgesToC11.push_back(nC2->at(i));
                // Edges {C12, N(C2)}
                if(!G->vertexHasEdgeTo(c12, nC2->at(i)))
                    addedEdgesToC12.push_back(nC2->at(i));
            }
            else
                throw std::invalid_argument("There shouldn't be any inactive vertices here");

//            std::cout << "u in N(c2): " << nC2->at(i) << " with deg(u): " << G->getVertexDegree(nC2->at(i)) << '\n';
        }

        alreadyInactive[v] = 1;
        alreadyInactive[c11] = 1;
        alreadyInactive[c12] = 1;
        alreadyInactive[c2] = 1;

        tempDeleted.push_back(v);
        tempNeighbours.push_back(c11);
        tempNeighbours.push_back(c12);
        tempNeighbours.push_back(c2);

        tempAddedEdges.push_back(addedEdgesToC11);
        tempAddedEdges.push_back(addedEdgesToC12);

        cnt++;

        delete nC2;
        delete neighbours;
    }

    if(cnt == 0) {
        if(G->deg3clique)
            std::cout << "Nothing was applied\n";
        return INAPPLICABLE;
    }

    int edgeCnt = 0;

    for (int i = 0; i < cnt; ++i) {
        Reduction* delVer = new Reduction(RULE::DEGREE_THREE_CLIQ, 0, new std::vector<int>(), new std::vector<int>());
        delVer->addedEdges = new std::vector<std::vector<int>>;
        delVer->rDepth = depth;
        delVer->kDecrement = 1;

        (*k) = (*k) -1;

        G->cntDeg3Clique++;

        int v = tempDeleted.at(i);

        int c11 = tempNeighbours.at(3*i + 0);
        int c12 = tempNeighbours.at(3*i + 1);
        int c2 = tempNeighbours.at(3*i + 2);

        std::vector<int> edgeToC11 = tempAddedEdges.at(2*i+0);
        std::vector<int> edgeToC12 = tempAddedEdges.at(2*i+1);

        for (int j = 0; j < (int)edgeToC11.size(); ++j) {
            if(G->getVertexDegree(edgeToC11.at(j)) > 0)
                G->addEdgeToVertex(c11, edgeToC11.at(j));
        }
        for (int j = 0; j < (int)edgeToC12.size(); ++j) {
            if(G->getVertexDegree(edgeToC12.at(j)) > 0)
                G->addEdgeToVertex(c12, edgeToC12.at(j));
        }

        delVer->deletedVertices->push_back(v);
        delVer->deletedVCVertices->push_back(c11);
        delVer->deletedVCVertices->push_back(c12);
        delVer->deletedVCVertices->push_back(c2);
        delVer->addedEdges->push_back(edgeToC11);
        delVer->addedEdges->push_back(edgeToC12);
        appliedRules->push_back(delVer);

        G->setInactive(v);
        G->setInactive(c2);

        edgeCnt += (int)edgeToC11.size();
        edgeCnt += (int)edgeToC12.size();

//        std::cout << "REDUCE: Deg 3 Clique" << std::endl;
//        std::cout << "Degree 3 vertex: " << delVer->deletedVertices->at(0) << '\n';
//        std::cout << delVer->deletedVCVertices->at(0) << ", " << delVer->deletedVCVertices->at(1) << ", " << delVer->deletedVCVertices->at(2) << '\n';
//        std::cout << "-----------" << '\n';

        if(G->deg3clique) {
            std::cout << "\nApplying rule to v = " << v << '\n';
            std::cout << "i = " << i << '\n';
            std::cout << "At recursion = " << delVer->rDepth << '\n';

            std::cout << "Edge 1:" << c11 << " with Neighbour of C2:" << '\n';
            for (int j = 0; j < (int)edgeToC11.size(); ++j) {
                std::cout << edgeToC11.at(j) << '\n';
            }

            std::cout << "Edge 2:" << c12 << " with Neighbour of C2:" << '\n';
            for (int j = 0; j < (int)edgeToC12.size(); ++j) {
                std::cout << edgeToC12.at(j) << '\n';
            }
        }
    }

    auto stopDeg3clique = std::chrono::high_resolution_clock::now();
    double Deg3clique = (std::chrono::duration_cast<std::chrono::microseconds>(stopDeg3clique - startDeg3clique).count() /  1000) / (double) 1000;

    if (printDebug) {
        std::string deg3cliqred = "Reduced "+ std::to_string(cnt) + " Deg3: 2-Clique-Neighbourhood: in ";
        deg3cliqred += std::to_string(Deg3clique)+ " seconds adding " + std::to_string(edgeCnt);
        deg3cliqred += " edges at rec Depth " + std::to_string(depth) + '\n';
        std::cout << '#' << deg3cliqred;
    }

    return APPLICABLE;
}

RULE_APPLICATION_RESULT Reductions::rule_DegreeThree_Domination(BucketGraph* G, int* k, int depth, bool deg2inc, bool checkBudget, bool printDebug)
{
    list<BucketVertex>* degThreeBucket = G->getVerticesOfDegree(3);



    if(degThreeBucket == nullptr || degThreeBucket->empty()) {
        if(G->deg3dom)
            std::cout << "\nNo degree 3 Vertex!" << '\n';
        return INAPPLICABLE; }

    if(G->deg3dom) {
        std::cout << "\nRule: Degree 3: Domination" << '\n';
        std::cout << "degThreeBucket: " << degThreeBucket << '\n';
        std::cout << "Deg 3 Bucket of size: " << (int)degThreeBucket->size() << '\n';
        std::cout << "---------" << '\n';
    }

    auto startDeg3dom = std::chrono::high_resolution_clock::now();

    int cnt = 0;
    int kDecrement = 0;

    std::vector<int> tempDeleted;
    std::vector<int> tempNeighbours;
    std::unordered_map<int, int> alreadyInactive;
    std::vector<bool> tempClique;
    std::vector<std::vector<int>> tempAddedEdges;

    for(auto it = degThreeBucket->begin(); it != degThreeBucket->end(); it++)
    {
        int v = it->index;

        if(!G->isActive(v)) {
            continue;
        }

        if(*k - kDecrement < 1 && checkBudget)
            return INSUFFICIENT_BUDGET;

        std::vector<int>* neighbours = G->getNeighbours(v);

        int a = neighbours->at(0);
        int b = neighbours->at(1);
        int c = neighbours->at(2);

        if(!G->isActive(a) || !G->isActive(b) || !G->isActive(c)) {
            delete neighbours;
            continue;
        }

        if(!G->vertexHasEdgeTo(a, b) && !G->vertexHasEdgeTo(a, c) && !G->vertexHasEdgeTo(b, c)) {
            if(G->deg3dom) {
                std::cout << "Neighbours are independent!\n";
                std::cout << a << ", " << b << ", " << c << '\n';
                std::cout << "---------" << '\n';
            }
            delete neighbours;
            continue;
        }

        if(alreadyInactive[v] == 1 || alreadyInactive[a] == 1 || alreadyInactive[b] == 1 || alreadyInactive[c] == 1)
        {
            if(G->deg3dom) {
                std::cout << "Already treated: " << v << ", maybe next round!\n";
                std::cout << a << ", " << b << ", " << c << '\n';
                std::cout << "---------" << '\n';
            }
            delete neighbours;
            continue;
        }

        int dominant = a;
        int c1 = b;
        int c2 = c;

        bool clique = false;

        bool connectionAB = G->vertexHasEdgeTo(a, b);
        bool connectionAC = G->vertexHasEdgeTo(a, c);
        bool connectionBC = G->vertexHasEdgeTo(b, c);

        if(connectionAB && connectionAC && connectionBC)
        {
            clique = true;
        }
        // Dom: a, C: b-c
        else if (connectionAB && connectionAC && !connectionBC)
        {
            dominant = a;
            c1 = b;
            c2 = c;
        }// Dom: b, C: a-c
        else if (connectionAB && !connectionAC && connectionBC)
        {
            c1 = a;
            c2 = c;
            dominant = b;
        }// Dom: c, C: b-a
        else if (!connectionAB && connectionAC && connectionBC)
        {
            c2 = b;
            dominant = c;
            c1 = a;
        }
        else{
            if(G->deg3dom) {
                std::cout << "2-Clique Neighbourhood: " << v << '\n';
                std::cout << a << ", " << b << ", " << c << '\n';
                std::cout << "---------" << '\n';
            }
            delete neighbours;
            continue;
        }

        if(G->deg3dom) {
            std::cout << "VALID = " << cnt+1 << '\n';
            std::string foundV = "Found vertex: " + std::to_string(v) + '\n';
            std::cout << ColorPrint::dye(foundV, 'y');

            std::cout << dominant << ", " << c1 << ", " << c2 << '\n';
            if (clique)
                std::cout << "CLIQUE!" << '\n';
            else {
                std::cout << "No Clique, goes to deg 2!" << '\n';
                std::cout << "With Dominant Vertex: " << dominant << '\n';
            }
            std::cout << ColorPrint::dye(lineDelimiter, 'y');
        }

        if(clique) {
            kDecrement += 3;
//            if(deg2inc){
//                tempAddedEdges.push_back({});
//                tempAddedEdges.push_back({});
//            }
        }
        else {
            kDecrement += 1;

            {
//            if(deg2inc)
//            {
//                std::vector<int> addedEdgesToC1;
//                std::vector<int> addedEdgesToC2;
//
//                std::vector<int> *nC1 = G->getNeighbours(c1);
//                std::vector<int> *nC2 = G->getNeighbours(c2);
//
//                bool untouchable = false;
//
//                for (int i = 0; i < (int) nC1->size(); ++i) {
//                    if (alreadyInactive[nC1->at(i)] == 1) {
//                        untouchable = true;
//                        break;
//                    }
//                }
//                if (untouchable) {
//                    delete nC1;
//                    delete nC2;
//                    continue;
//                }
//                for (int i = 0; i < (int) nC2->size(); ++i) {
//                    if (alreadyInactive[nC2->at(i)] == 1) {
//                        untouchable = true;
//                        break;
//                    }
//                }
//                if (untouchable) {
//                    delete nC1;
//                    delete nC2;
//                    continue;
//                }
//
//                // Edges {c1, N(c2)}
//                for (int i = 0; i < (int) nC2->size(); ++i) {
//                    if (G->isActive(nC2->at(i)) && !G->vertexHasEdgeTo(c1, nC2->at(i))) {
//                        addedEdgesToC1.push_back(nC2->at(i));
//                    }
//                }
//
//                // Edges {c2, N(c1)}
//                for (int i = 0; i < (int) nC1->size(); ++i) {
//                    if (G->isActive(nC1->at(i)) && !G->vertexHasEdgeTo(c2, nC1->at(i))) {
//                        addedEdgesToC2.push_back(nC1->at(i));
//                    }
//                }
//
//                delete nC1;
//                delete nC2;
//
//                tempAddedEdges.push_back(addedEdgesToC1);
//                tempAddedEdges.push_back(addedEdgesToC2);
//            }
        }
        }

        alreadyInactive[it->index] = 1;
        alreadyInactive[dominant] = 1;
        alreadyInactive[c1] = 1;
        alreadyInactive[c2] = 1;

        tempDeleted.push_back(v);
        tempNeighbours.push_back(dominant);
        tempNeighbours.push_back(c1);
        tempNeighbours.push_back(c2);

        tempClique.push_back(clique);

        delete neighbours;

        cnt++;

    }

    if(tempDeleted.empty() || cnt == 0) {
        if(G->deg3dom)
            std::cout << "Nothing was applied\n";
        return INAPPLICABLE;
    }

    int kAdded = 0;
    int edgeCnt = 0;

    for (int i = 0; i < cnt; ++i) {
        Reduction* delVer = new Reduction(RULE::DEGREE_THREE_DOM, 0, new std::vector<int>(), new std::vector<int>());
        delVer->rDepth = depth;

        int v = tempDeleted.at(i);
        delVer->deletedVertices->push_back(v);

        int dom = tempNeighbours.at(3*i + 0);
        int c1 = tempNeighbours.at(3*i + 1);
        int c2 = tempNeighbours.at(3*i + 2);

        // Case 1: Clique
        if (tempClique.at(i))
        {
            G->cntDeg3Dom2++;
            delVer->kDecrement = 3;
            G->setInactive(c1);
            G->setInactive(c2);
        }
        else
        {
            delVer->kDecrement = 1;
            G->cntDeg3Dom1++;

            {//            if(deg2inc)
//            {
//                delVer->addedEdges = new std::vector <std::vector<int>>;
//
//                std::vector<int> edgeToC1 = tempAddedEdges.at(2 * i + 0);
//                std::vector<int> edgeToC2 = tempAddedEdges.at(2 * i + 1);
//
//                for (int j = 0; j < (int) edgeToC1.size(); ++j) {
//                    if (G->getVertexDegree(edgeToC1.at(j)) != 0)
//                        G->addEdgeToVertex(c1, edgeToC1.at(j));
//                }
//
//                for (int j = 0; j < (int) edgeToC2.size(); ++j) {
//                    if (G->getVertexDegree(edgeToC2.at(j)) != 0)
//                        G->addEdgeToVertex(c2, edgeToC2.at(j));
//                }
//
//                edgeCnt += (int) edgeToC1.size();
//                edgeCnt += (int) edgeToC2.size();
//
//                delVer->addedEdges->push_back(edgeToC1);
//                delVer->addedEdges->push_back(edgeToC2);
//            }
            }
        }
        G->setInactive(dom);

        delVer->deletedVCVertices->push_back(dom);
        delVer->deletedVCVertices->push_back(c1);
        delVer->deletedVCVertices->push_back(c2);
        appliedRules->push_back(delVer);

        (*k) = (*k) - delVer->kDecrement;
        kAdded += delVer->kDecrement;

        if(G->deg3dom) {
            std::cout << "\nApplying rule to v = " << v << '\n';
            std::cout << "At recursion = " << delVer->rDepth << '\n';
            std::cout << "k = " << *k << '\n';
            if (tempClique.at(i))
                std::cout << "Neighbours a clique " << '\n';
            else {
                std::cout << "Neighbours no clique " << '\n';
                std::cout << dom << ", " << c1 << ", " << c2 << '\n';
            }
            std::cout << "---------" << '\n';
        }

    }

    auto stopDeg3dom = std::chrono::high_resolution_clock::now();
    double Deg3dom = (std::chrono::duration_cast<std::chrono::microseconds>(stopDeg3dom - startDeg3dom).count() /  1000) / (double) 1000;


    if (printDebug) {
        std::string deg3domred = "Reduced "+ std::to_string(cnt) + " Deg3: Domination: in ";
        deg3domred += std::to_string(Deg3dom)+ " seconds adding " + std::to_string(kAdded);
        deg3domred += " to VC at rec depth: " + std::to_string(depth) + '\n';
//        std::cout << '#' << ColorPrint::dye(deg3domred, 'p');
        std::cout << '#' << deg3domred;
    }

    return APPLICABLE;
}

//using namespace std;

RULE_APPLICATION_RESULT Reductions::rule_DegreeFour_Clique(BucketGraph* G, int* k, int depth, bool checkBudget, bool printDebug)
{
    list<BucketVertex>* degFourBucket = G->getVerticesOfDegree(4);

    if(degFourBucket == nullptr || degFourBucket->empty()) {
//        if(G->deg3clique)
//            std::cout << "\nNo degree 3 Vertex!" << '\n';
        return INAPPLICABLE; }
//        return INAPPLICABLE;

    if(G->deg4clique) {
        std::cout << "\nRule: Degree 4: 2-Clique-Neighborhood" << '\n';
        std::cout << "Deg 4 Bucket of size: " << (int)degFourBucket->size() << '\n';
        std::cout << "---------" << '\n';
    }

    auto startDeg4clique = std::chrono::high_resolution_clock::now();

    int cnt = 0;
    int kDecrement = 0;

    std::vector<int> tempDeleted;
    std::vector<int> tempNeighbours;
    std::vector<int> tempDecrement;
    std::unordered_map<int, int> alreadyInactive;
    std::vector<std::vector<int>> tempAddedEdges;

    for(auto it = degFourBucket->begin(); it != degFourBucket->end(); it++)
    {
        if(*k - cnt < 1 && checkBudget)
            return INSUFFICIENT_BUDGET;

        int v = it->index;

        if(!G->isActive(v)) {
            continue;
        }

        std::vector<int>* neighbours = G->getNeighbours(v);

        int a = neighbours->at(0);
        int b = neighbours->at(1);
        int c = neighbours->at(2);
        int d = neighbours->at(3);

        if(!G->isActive(a) || !G->isActive(b) || !G->isActive(c) || !G->isActive(d)) {
            delete neighbours;
            continue;
        }

        if(alreadyInactive[v] == 1 || alreadyInactive[a] == 1 || alreadyInactive[b] == 1 || alreadyInactive[c] == 1 || alreadyInactive[d] == 1)
        {
//            if(G->deg4clique) {
//                std::cout << "Already treated: " << v << ", maybe next round!\n";
//                std::cout << a << ", " << b << ", " << c << '\n';
//                std::cout << "---------" << '\n';
//            }
            delete neighbours;
            continue;
        }

        int connections = 0;
        int nodesConnection[4] = {0, 0, 0, 0};

        bool connectionAB = G->vertexHasEdgeTo(a, b);
        bool connectionAC = G->vertexHasEdgeTo(a, c);
        bool connectionAD = G->vertexHasEdgeTo(a, d);
        bool connectionBC = G->vertexHasEdgeTo(b, c);
        bool connectionBD = G->vertexHasEdgeTo(b, d);
        bool connectionCD = G->vertexHasEdgeTo(c, d);

        if(connectionAB) {nodesConnection[0]++; nodesConnection[1]++; connections++; if(G->deg4clique)std::cout << "CONNECTION AB\n";}
        if(connectionAC) {nodesConnection[0]++; nodesConnection[2]++; connections++; if(G->deg4clique)std::cout << "CONNECTION AC\n";}
        if(connectionAD) {nodesConnection[0]++; nodesConnection[3]++; connections++; if(G->deg4clique)std::cout << "CONNECTION AD\n";}
        if(connectionBC) {nodesConnection[1]++; nodesConnection[2]++; connections++; if(G->deg4clique)std::cout << "CONNECTION BC\n";}
        if(connectionBD) {nodesConnection[1]++; nodesConnection[3]++; connections++; if(G->deg4clique)std::cout << "CONNECTION BD\n";}
        if(connectionCD) {nodesConnection[2]++; nodesConnection[3]++; connections++; if(G->deg4clique)std::cout << "CONNECTION CD\n";}

        // Case |C1| > |C2|
        if(connections == 3)
        {
            if(G->deg4clique)
                std::cout << "Case 3\n";

            bool cliquePresent = false;
            for (int j = 0; j < 4; ++j)
                if(nodesConnection[j] == 0)
                    cliquePresent = true;

            if(!cliquePresent)
            {
                if(G->deg4clique)
                std::cout << "No clique present!\n";
                delete neighbours;
                continue;
            }

            int c2;
            std::vector<int> C1;

            for (int j = 0; j < 4; ++j) {
                if(j == 0){
                    if(nodesConnection[j] > 0)
                        C1.push_back(a);
                    else
                        c2 = a;
                }else if(j == 1){
                    if(nodesConnection[j] > 0)
                        C1.push_back(b);
                    else
                        c2 = b;
                }else if(j == 2){
                    if(nodesConnection[j] > 0)
                        C1.push_back(c);
                    else
                        c2 = c;
                }else if(j == 3){
                    if(nodesConnection[j] > 0)
                        C1.push_back(d);
                    else
                        c2 = d;
                }
            }
            if((int)C1.size() != 3)
                throw std::invalid_argument("Wrong thought process");

            std::vector<int>* nC2 = G->getNeighbours(c2);

            std::vector<int> addedEdgesToC11;
            std::vector<int> addedEdgesToC12;
            std::vector<int> addedEdgesToC13;

            if(!isTouchable(G, nC2, alreadyInactive))
            {
                delete nC2;
                delete neighbours;
                continue;
            }

            // Adding edges to N(C2)
            for (int i = 0; i < (int)nC2->size(); ++i) {

                if(G->isActive(nC2->at(i))){ // && G->getVertexDegree(nC2->at(i)) > 0){
                    // Edges {C11, N(C2)}
                    if(!G->vertexHasEdgeTo(C1.at(0), nC2->at(i)))
                        addedEdgesToC11.push_back(nC2->at(i));
                    // Edges {C12, N(C2)}
                    if(!G->vertexHasEdgeTo(C1.at(1), nC2->at(i)))
                        addedEdgesToC12.push_back(nC2->at(i));
                    // Edges {C13, N(C2)}
                    if(!G->vertexHasEdgeTo(C1.at(2), nC2->at(i)))
                        addedEdgesToC13.push_back(nC2->at(i));
                }
                else
                    throw std::invalid_argument("There shouldn't be any inactive vertices here");

//            std::cout << "u in N(c2): " << nC2->at(i) << " with deg(u): " << G->getVertexDegree(nC2->at(i)) << '\n';
            }

            kDecrement += 1;

            tempNeighbours.push_back(C1.at(0));
            tempNeighbours.push_back(C1.at(1));
            tempNeighbours.push_back(C1.at(2));
            tempNeighbours.push_back(c2);

            tempAddedEdges.push_back(addedEdgesToC11);
            tempAddedEdges.push_back(addedEdgesToC12);
            tempAddedEdges.push_back(addedEdgesToC13);

            tempDecrement.push_back(1);

            delete nC2;
        }
        // Case |C1| = |C2|
        else if(connections == 4)
        {
            if(G->deg4clique)
                std::cout << "Case 4\n";

            bool cliquePresent = true;
            for (int j = 0; j < 4; ++j)
                if(nodesConnection[j] == 3)
                    cliquePresent = false;

            // Case domination
            if(!cliquePresent)
            {
                if(G->deg4clique)
                std::cout << "No clique present!\n";
                delete neighbours;
                continue;
            }

            for (int j = 0; j < 4; ++j)
                if(nodesConnection[j] != 2)
                    throw std::invalid_argument("Wrong interpretation");

            std::vector<int>* nA = G->getNeighbours(a);
            std::vector<int>* nB = G->getNeighbours(b);
            std::vector<int>* nC = G->getNeighbours(c);
            std::vector<int>* nD = G->getNeighbours(d);

            if(!isTouchable(G, nA, alreadyInactive)
            || !isTouchable(G, nB, alreadyInactive)
            || !isTouchable(G, nC, alreadyInactive)
            || !isTouchable(G, nD, alreadyInactive))
            {
                delete nA;
                delete nB;
                delete nC;
                delete nD;
                delete neighbours;
                continue;
            }

            std::vector<int> addedEdgesToC1_1;
            std::vector<int> addedEdgesToC1_2;

            int c1_1 = a;
            int c1_2;
            int c2_1;
            int c2_2;

            // Define cliques 1 and 2, c2_1 being the non edge to c1_1 and c2_2 to c1_2
            if(connectionAB && connectionAC) {
                c1_2 = b;
                c2_1 = d;
                c2_2 = c;
            }
            else if (connectionAB && connectionAD){
                c1_2 = b;
                c2_1 = c;
                c2_2 = d;
            }
            else if (connectionAC && connectionAD){
                c1_2 = c;
                c2_1 = b;
                c2_2 = d;
            }
            else
                throw std::invalid_argument("Wrong connections");

            if(G->deg4clique){
                std::cout << "a = " << a << ", b = " << b << ", c = " << c << ", d = " << d << '\n';
                std::cout << "c1_1 = " << c1_1 << '\n';
                std::cout << "c1_2 = " << c1_2 << '\n';
                std::cout << "c2_1 = " << c2_1 << '\n';
                std::cout << "c2_2 = " << c2_2 << '\n';
            }

            // Neighbours of the non-edge partners
            std::vector<int>* nC2_1;
            std::vector<int>* nC2_2;

            if(c2_1 == b) { nC2_1 = nB; nC2_2 = nD;}
            else if(c2_1 == c) { nC2_1 = nC; nC2_2 = nD;}
            else if(c2_1 == d) { nC2_1 = nD; nC2_2 = nC;}

            // Adding edges to N(C2_1)
            for (int i = 0; i < (int)nC2_1->size(); ++i) {
                if(G->isActive(nC2_1->at(i))){ // && G->getVertexDegree(nC2->at(i)) > 0){
                    // Edges {C11, N(C21)}
                    if(!G->vertexHasEdgeTo(c1_1, nC2_1->at(i)))
                        addedEdgesToC1_1.push_back(nC2_1->at(i));
                }
                else
                    throw std::invalid_argument("There shouldn't be any inactive vertices here");
            }
            // Adding edges to N(C2_2)
            for (int i = 0; i < (int)nC2_2->size(); ++i) {
                if(G->isActive(nC2_2->at(i))){ // && G->getVertexDegree(nC2->at(i)) > 0){
                    // Edges {C12, N(C22)}
                    if(!G->vertexHasEdgeTo(c1_2, nC2_2->at(i)))
                        addedEdgesToC1_2.push_back(nC2_2->at(i));
                }
                else
                    throw std::invalid_argument("There shouldn't be any inactive vertices here");
            }

            kDecrement += 2;

            tempNeighbours.push_back(c1_1);
            tempNeighbours.push_back(c1_2);
            tempNeighbours.push_back(c2_1);
            tempNeighbours.push_back(c2_2);

            tempAddedEdges.push_back(addedEdgesToC1_1);
            tempAddedEdges.push_back(addedEdgesToC1_2);
            tempAddedEdges.push_back({});

            tempDecrement.push_back(2);

            delete nA;
            delete nB;
            delete nC;
            delete nD;
        }
        else
        {
            if(G->deg4clique) {
                    std::cout << "This RULE doesn't apply, v = " << v << '\n';
                    std::cout << a << ", " << b << ", " << c << ", " << d << '\n';
                    std::cout << "---------" << '\n';
            }
            delete neighbours;
            continue;
        }

        if(G->deg4clique) {
            std::cout << ColorPrint::dye(lineDelimiter, 'y');
            std::cout << "VALID = " << cnt+1 << '\n';
            std::string foundV = "Found vertex: " + std::to_string(v) + '\n';
            if(connections == 3)
                std::cout << "CASE 3!" << '\n';
            else
                std::cout << "CASE 4!" << '\n';
            std::cout << ColorPrint::dye(foundV, 'y');
            std::cout << a << ", " << b << ", " << c << ", " << d << '\n';
            std::cout << ColorPrint::dye(lineDelimiter, 'y');
        }

        alreadyInactive[v] = 1;
        alreadyInactive[a] = 1;
        alreadyInactive[b] = 1;
        alreadyInactive[c] = 1;
        alreadyInactive[d] = 1;

        tempDeleted.push_back(v);

        cnt++;

        delete neighbours;
    }

    if(cnt == 0) {
        if(G->deg4clique)
            std::cout << "Nothing was applied\n";
        return INAPPLICABLE;
    }

    if(G->deg4clique)
        std::cout << "Found " << cnt << " instances!\n";


    int edgeCnt = 0;

    for (int i = 0; i < cnt; ++i) {
        Reduction* delVer = new Reduction(RULE::DEGREE_FOUR_CLIQUE, 0, new std::vector<int>(), new std::vector<int>());
        delVer->addedEdges = new std::vector<std::vector<int>>;
        delVer->rDepth = depth;

        if(G->deg4clique) {
            std::cout << "Adding to applied rules \n";
        }

        G->cntDeg4Clique++;

        int v = tempDeleted.at(i);

        // Case 3 Connections
        if(tempDecrement.at(i) == 1)
        {
            int c1_1 = tempNeighbours.at(4*i + 0);
            int c1_2 = tempNeighbours.at(4*i + 1);
            int c1_3 = tempNeighbours.at(4*i + 2);
            int c2 = tempNeighbours.at(4*i + 3);

            if(G->deg4clique) {
                std::cout << "Case 3: The order is: \n";
                std::cout << c1_1 << ", " << c1_2 << ", " << c1_3 << ", " << c2 << '\n';
            }

            std::vector<int> edgeToC1_1 = tempAddedEdges.at(3*i+0);
            std::vector<int> edgeToC1_2 = tempAddedEdges.at(3*i+1);
            std::vector<int> edgeToC1_3 = tempAddedEdges.at(3*i+2);

            for (int j = 0; j < (int)edgeToC1_1.size(); ++j) {
                if(G->getVertexDegree(edgeToC1_1.at(j)) > 0)
                    G->addEdgeToVertex(c1_1, edgeToC1_1.at(j));
            }
            for (int j = 0; j < (int)edgeToC1_2.size(); ++j) {
                if(G->getVertexDegree(edgeToC1_2.at(j)) > 0)
                    G->addEdgeToVertex(c1_2, edgeToC1_2.at(j));
            }
            for (int j = 0; j < (int)edgeToC1_3.size(); ++j) {
                if(G->getVertexDegree(edgeToC1_3.at(j)) > 0)
                    G->addEdgeToVertex(c1_3, edgeToC1_3.at(j));
            }

            delVer->deletedVCVertices->push_back(c1_1);
            delVer->deletedVCVertices->push_back(c1_2);
            delVer->deletedVCVertices->push_back(c1_3);
            delVer->deletedVCVertices->push_back(c2);
            delVer->addedEdges->push_back(edgeToC1_1);
            delVer->addedEdges->push_back(edgeToC1_2);
            delVer->addedEdges->push_back(edgeToC1_3);
            delVer->kDecrement = 1;
            (*k) = (*k) -1;

            G->setInactive(c2);
        }
        // Case 4 Connection
        else if(tempDecrement.at(i) == 2){
            int c1_1 = tempNeighbours.at(4*i + 0);
            int c1_2 = tempNeighbours.at(4*i + 1);
            int c2_1 = tempNeighbours.at(4*i + 2);
            int c2_2 = tempNeighbours.at(4*i + 3);

            if(G->deg4clique) {
                std::cout << "Case 4: The order is: \n";
                std::cout << c1_1 << ", " << c1_2 << ", " << c2_1 << ", " << c2_2 << '\n';
            }

            std::vector<int> edgeToC1_1 = tempAddedEdges.at(3*i+0);
            std::vector<int> edgeToC1_2 = tempAddedEdges.at(3*i+1);

            for (int j = 0; j < (int)edgeToC1_1.size(); ++j) {
                if(G->getVertexDegree(edgeToC1_1.at(j)) > 0) {
                    G->addEdgeToVertex(c1_1, edgeToC1_1.at(j));
                }
            }
            for (int j = 0; j < (int)edgeToC1_2.size(); ++j) {
                if(G->getVertexDegree(edgeToC1_2.at(j)) > 0) {
                    G->addEdgeToVertex(c1_2, edgeToC1_2.at(j));
                }
            }

            delVer->deletedVCVertices->push_back(c1_1);
            delVer->deletedVCVertices->push_back(c1_2);
            delVer->deletedVCVertices->push_back(c2_1);
            delVer->deletedVCVertices->push_back(c2_2);
            delVer->addedEdges->push_back(edgeToC1_1);
            delVer->addedEdges->push_back(edgeToC1_2);
//            delVer->addedEdges->push_back(edgeToC1_3);
            delVer->kDecrement = 2;
            (*k) = (*k) -2;


            G->setInactive(c2_1);
            G->setInactive(c2_2);
        }

        delVer->deletedVertices->push_back(v);
        appliedRules->push_back(delVer);

        G->setInactive(v);

//        if(G->deg4clique) {
//            std::cout << "\nApplying rule to v = " << v << '\n';
//            std::cout << "i = " << i << '\n';
//            std::cout << "At recursion = " << delVer->rDepth << '\n';
//
//            std::cout << "Edge 1:" << c11 << " with Neighbour of C2:" << '\n';
//            for (int j = 0; j < (int)edgeToC11.size(); ++j) {
//                std::cout << edgeToC11.at(j) << '\n';
//            }
//
//            std::cout << "Edge 2:" << c12 << " with Neighbour of C2:" << '\n';
//            for (int j = 0; j < (int)edgeToC12.size(); ++j) {
//                std::cout << edgeToC12.at(j) << '\n';
//            }
//        }
    }

    auto stopDeg4clique = std::chrono::high_resolution_clock::now();
    double Deg4clique = (std::chrono::duration_cast<std::chrono::microseconds>(stopDeg4clique - startDeg4clique).count() /  1000) / (double) 1000;

    if (printDebug) {
        std::string deg4clique = "Reduced "+ std::to_string(cnt) + " Deg3: 2-Clique-Neighbourhood: in ";
        deg4clique += std::to_string(Deg4clique)+ " seconds at rec Depth " + std::to_string(depth) + '\n';
        std::cout << '#' << deg4clique;
    }

    return APPLICABLE;
}

RULE_APPLICATION_RESULT Reductions::rule_DegreeFour_Domination(BucketGraph *G, int *k, int depth, bool checkBudget, bool printDebug) {
    list<BucketVertex>* degFourBucket = G->getVerticesOfDegree(4);

    if(degFourBucket == nullptr || degFourBucket->empty()) {
//        if(G->deg3clique)
//            std::cout << "\nNo degree 3 Vertex!" << '\n';
        return INAPPLICABLE; }
//        return INAPPLICABLE;

    if(G->deg4dom) {
        std::cout << "\nRule: Degree 4: 2-Clique-Neighborhood" << '\n';
        std::cout << "Deg 4 Bucket of size: " << (int)degFourBucket->size() << '\n';
        std::cout << "---------" << '\n';
    }

    auto startDeg4clique = std::chrono::high_resolution_clock::now();

    int cnt = 0;
    int kDecrement = 0;

    std::vector<int> tempDeleted;
    std::vector<int> tempNeighbours;
    std::vector<int> tempDecrement;
    std::vector<int> tempDoms;
    std::unordered_map<int, int> alreadyInactive;

    for(auto it = degFourBucket->begin(); it != degFourBucket->end(); it++)
    {
        if(*k - cnt < 1 && checkBudget)
            return INSUFFICIENT_BUDGET;

        int v = it->index;

        if(!G->isActive(v)) {
            continue;
        }

        std::vector<int>* neighbours = G->getNeighbours(v);

        int a = neighbours->at(0);
        int b = neighbours->at(1);
        int c = neighbours->at(2);
        int d = neighbours->at(3);

        if(!G->isActive(a) || !G->isActive(b) || !G->isActive(c) || !G->isActive(d)) {
            delete neighbours;
            continue;
        }

        if(alreadyInactive[v] == 1 || alreadyInactive[a] == 1 || alreadyInactive[b] == 1 || alreadyInactive[c] == 1 || alreadyInactive[d] == 1)
        {
//            if(G->deg4clique) {
//                std::cout << "Already treated: " << v << ", maybe next round!\n";
//                std::cout << a << ", " << b << ", " << c << '\n';
//                std::cout << "---------" << '\n';
//            }
            delete neighbours;
            continue;
        }

        int connections = 0;
        int nodesConnection[4] = {0, 0, 0, 0};

        bool connectionAB = G->vertexHasEdgeTo(a, b);
        bool connectionAC = G->vertexHasEdgeTo(a, c);
        bool connectionAD = G->vertexHasEdgeTo(a, d);
        bool connectionBC = G->vertexHasEdgeTo(b, c);
        bool connectionBD = G->vertexHasEdgeTo(b, d);
        bool connectionCD = G->vertexHasEdgeTo(c, d);

        if(connectionAB) {nodesConnection[0]++; nodesConnection[1]++; connections++; if(G->deg4dom)std::cout << "CONNECTION AB\n";}
        if(connectionAC) {nodesConnection[0]++; nodesConnection[2]++; connections++; if(G->deg4dom)std::cout << "CONNECTION AC\n";}
        if(connectionAD) {nodesConnection[0]++; nodesConnection[3]++; connections++; if(G->deg4dom)std::cout << "CONNECTION AD\n";}
        if(connectionBC) {nodesConnection[1]++; nodesConnection[2]++; connections++; if(G->deg4dom)std::cout << "CONNECTION BC\n";}
        if(connectionBD) {nodesConnection[1]++; nodesConnection[3]++; connections++; if(G->deg4dom)std::cout << "CONNECTION BD\n";}
        if(connectionCD) {nodesConnection[2]++; nodesConnection[3]++; connections++; if(G->deg4dom)std::cout << "CONNECTION CD\n";}

        if(connections < 3){
            if(G->deg4dom)
                std::cout << "No dominator present!\n";
            delete neighbours;
            continue;
        }

        bool dominatorPresent = false;
        std::vector<int> dominators;
        std::vector<int> subs;

        for (int j = 0; j < 4; ++j) {
            int domIdx;
            if(j == 0)
                domIdx = a;
            else if(j == 1)
                domIdx = b;
            else if(j == 2)
                domIdx = c;
            else if(j == 3)
                domIdx = d;
            else
                throw std::invalid_argument("Wrong dom id");
//            std::cout << domIdx << ": connections = " << nodesConnection[j] << "\n";

            if (nodesConnection[j] == 3) {
                dominatorPresent = true;
                dominators.push_back(domIdx);
            }
            else
                subs.push_back(domIdx);
        }
        if(!dominatorPresent)
        {
            if(G->deg4dom)
                std::cout << "No dominator present!\n";
            delete neighbours;
            continue;
        }

        int numDoms = (int)dominators.size();
        kDecrement += numDoms;


        for (int j = 0; j < numDoms; ++j) {
//            std::cout << dominators.at(j) <<  ": Dominator!\n";
            tempDoms.push_back(dominators.at(j));
        }
        for (int j = 0; j < 4-numDoms; ++j) {
//            std::cout << subs.at(j) <<  ": Sub!\n";
            tempDoms.push_back(subs.at(j));
        }

        tempNeighbours.push_back(a);
        tempNeighbours.push_back(b);
        tempNeighbours.push_back(c);
        tempNeighbours.push_back(d);

        tempDeleted.push_back(v);
        tempDecrement.push_back(numDoms);

        alreadyInactive[v] = 1;
        alreadyInactive[a] = 1;
        alreadyInactive[b] = 1;
        alreadyInactive[c] = 1;
        alreadyInactive[d] = 1;



        if(G->deg4dom) {
            std::cout << ColorPrint::dye(lineDelimiter, 'y');
            std::cout << "VALID = " << cnt+1 << '\n';
            std::string foundV = "Found " + std::to_string(numDoms) + " dominators: " + std::to_string(v) + "\n";
            for (int i = 0; i < numDoms; ++i) {
                foundV += std::to_string(tempDoms.at(i)) + " ";
            }
            foundV += '\n';
            std::cout << ColorPrint::dye(foundV, 'y');
            std::cout << a << ", " << b << ", " << c << ", " << d << '\n';
            std::cout << ColorPrint::dye(lineDelimiter, 'y');
        }

        cnt++;
        delete neighbours;
    }

    if(cnt == 0) {
        if(G->deg4clique)
            std::cout << "Nothing was applied\n";
        return INAPPLICABLE;
    }

    if(G->deg4clique)
        std::cout << "Found " << cnt << " instances!\n";


    int edgeCnt = 0;

    for (int i = 0; i < cnt; ++i) {
        Reduction* delVer = new Reduction(RULE::DEGREE_FOUR_DOM, 0, new std::vector<int>(), new std::vector<int>());
        delVer->addedEdges = new std::vector<std::vector<int>>;
        delVer->rDepth = depth;

        if(G->deg4dom) {
            std::cout << "Adding to applied rules \n";
        }

        G->cntDeg4Clique++;

        int v = tempDeleted.at(i);

        int numDoms = tempDecrement.at(i);

//        std::cout << numDoms <<  " Dominators!\n";
        for (int j = 0; j < 4; ++j) {
            if(j < numDoms) {
                G->setInactive(tempDoms.at(4*i + j));
//                std::cout << " Dominators: " << tempDoms.at(4*i + j) << "\n";
            }
//            else
//                std::cout << " Subs: " << tempDoms.at(4*i + j) << "\n";

            delVer->deletedVCVertices->push_back(tempDoms.at(4*i + j));
        }
        if(numDoms == 4)
            G->setInactive(v);

        delVer->kDecrement = numDoms;
        (*k) = (*k) -numDoms;

        delVer->deletedVertices->push_back(v);
        appliedRules->push_back(delVer);

//        if(G->deg4clique) {
//            std::cout << "\nApplying rule to v = " << v << '\n';
//            std::cout << "i = " << i << '\n';
//            std::cout << "At recursion = " << delVer->rDepth << '\n';
//
//            std::cout << "Edge 1:" << c11 << " with Neighbour of C2:" << '\n';
//            for (int j = 0; j < (int)edgeToC11.size(); ++j) {
//                std::cout << edgeToC11.at(j) << '\n';
//            }
//
//            std::cout << "Edge 2:" << c12 << " with Neighbour of C2:" << '\n';
//            for (int j = 0; j < (int)edgeToC12.size(); ++j) {
//                std::cout << edgeToC12.at(j) << '\n';
//            }
//        }
    }

    auto stopDeg4clique = std::chrono::high_resolution_clock::now();
    double Deg4clique = (std::chrono::duration_cast<std::chrono::microseconds>(stopDeg4clique - startDeg4clique).count() /  1000) / (double) 1000;

    if (printDebug) {
        std::string deg4clique = "Reduced "+ std::to_string(cnt) + " Deg4: Domination: in ";
        deg4clique += std::to_string(Deg4clique)+ " seconds at rec Depth " + std::to_string(depth) + '\n';
        std::cout << '#' << deg4clique;
    }

    return APPLICABLE;
}

RULE_APPLICATION_RESULT Reductions::rule_LPFlow(BucketGraph* G, int* k, int depth, bool checkBudget, bool printDebug)
{
    if(!G->LP_INITIALISED)
        throw std::invalid_argument("LP data structures not initialised");
    std::vector<int>* delVertices = new std::vector<int>();
    std::vector<int>* delVCVertices = new std::vector<int>();
    auto startHop = std::chrono::high_resolution_clock::now();
    G->hopcroftKarpMatchingSize();
    auto stopHop = std::chrono::high_resolution_clock::now();
    auto startCom = std::chrono::high_resolution_clock::now();
    G->setBipartMatchingFlowComponentsInactive(delVertices, delVCVertices, *k, INT32_MAX/* 0.5 */);
    auto stopCom = std::chrono::high_resolution_clock::now();
    double Hop = (std::chrono::duration_cast<std::chrono::microseconds>(stopHop - startHop).count() / 1000) / (double) 1000;
    double Com = (std::chrono::duration_cast<std::chrono::microseconds>(stopCom - startCom).count() / 1000) / (double) 1000;
    if (printDebug)
        std::cout << "#Reduced " << delVertices->size() + delVCVertices->size() << " LP vertices (" << delVCVertices->size() << " VC / " << delVertices->size() << " non-VC) in " << Hop + Com << " seconds (HopcroftKarp: " << Hop << " + Connected Components: " << Com << ")" << '\n';

    Reduction* reduction = new Reduction(RULE::LPFLOW, delVCVertices->size(), delVertices, delVCVertices);
    reduction->rDepth = depth;
    if((int) delVCVertices->size() == 0 && (int) delVertices->size() == 0)
    {
        //std::cout << "LPFlow: INAPPLICABLE" << '\n';
        delete delVertices;
        delete delVCVertices;
        delete reduction;
        return INAPPLICABLE;
    }
    if(*k - (int) delVCVertices->size() < 0 && checkBudget)
    {
        //std::cout << "LPFlow: INSUFFICIENT_BUDGET" << '\n';
        delete delVertices;
        delete delVCVertices;
        delete reduction;
        return INSUFFICIENT_BUDGET;
    }
    *k = *k - delVCVertices->size();
    appliedRules->push_back(reduction);
    return APPLICABLE;
}

RULE_APPLICATION_RESULT Reductions::rule_Domination_BE(BucketGraph* G, int* k, int depth, bool checkBudget)
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

/* RULE_APPLICATION_RESULT Reductions::rule_Unconfined(BucketGraph* G, int* k, bool checkBudget)
{
    int numAttempts = 10;
    bool applied = false;
    std::vector<bool> inS = std::vector<bool>(G->getTotalNumVertices());
    std::vector<bool> inSNeighbours = std::vector<bool>(G->getTotalNumVertices());
    for(int i=0; i<(int) inS.size(); i++) { inS[i] = false; inSNeighbours[i] = false; }

    std::vector<int> SNeighbours = std::vector<int>();

    Reduction* reduction = new Reduction(RULE::DEGREE_ONE, 0, nullptr, new std::vector<int>());
    for (int i = 0; i<numAttempts && i<(int) inS.size(); i++)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distr(0, G->getTotalNumVertices()-1);
        int v = (int) distr(gen);
        if(!G->isActive(v)) { continue; }
        Vertex* vertex = G->getVertex(v);

        // S and neighbours of S list, that are kept up to date
        for(int i=0; i<(int) inS.size(); i++) { inS[i] = false; inSNeighbours[i] = false; }
        SNeighbours.clear();
        inS[vertex->getIndex()] = true;
        for (auto neighbour = vertex->adj->begin(); neighbour != vertex->adj->end(); ++neighbour)
        {
            if(!G->isActive(*neighbour)) { continue; }
            inSNeighbours[*neighbour] = true;
            SNeighbours.push_back(*neighbour);
        }

        // search continuation loop
        while(true)
        {
            int best = -1;
            int bestOutsideNeighboursSize = -1;
            // find best neighbour
            for (auto u = SNeighbours.begin(); u != SNeighbours.end(); ++u)
            {
                if(!G->getVertex(*u)->getActive()) { continue; }
                if(inS[*u]) { continue; }
                bool valid = true;
                int sharedNeighbours = 0;
                int outsideNeighboursSize = 0;
                for (auto it=G->getVertex(*u)->adj->begin(); it != G->getVertex(*u)->adj->end(); ++it)
                {
                    if(!G->getVertex(*u)->getActive()) { continue; }
                    if(!inSNeighbours[*it] && !inS[*it])
                    {
                        outsideNeighboursSize++;
                    }
                    // if u has more than one neighbour in S, u does not fit the criteria
                    if (inS[*it]) {
                        sharedNeighbours++;
                        if(sharedNeighbours > 1) { valid = false; break; }
                    }
                }

                if(!valid) { continue; }
                // update best/bestNeighbourhoodSize if u is better than the current best
                //std::cout << cp::dye("Checking best=" + std::to_string(best) + " with " + std::to_string(bestUniqueNeighboursSize) + " vs. u=" + std::to_string(*u) + " with " + std::to_string(uniqueNeighbours.size()), 'y') << std::endl;
                if(best != -1 && outsideNeighboursSize >= bestOutsideNeighboursSize) { continue; }
                best = *u;
                bestOutsideNeighboursSize = outsideNeighboursSize;
            }

            // if no expansion vertex found, vertex "vertex" is not unconfined (continue with next vertex)
            if(best == -1) { break; }
            if(bestOutsideNeighboursSize == 0)
            {
                if((*k) <= 0 && checkBudget) {
                    (*k) = (*k) + reduction->kDecrement;
                    G->setActive(reduction->deletedVCVertices);
                    delete reduction->deletedVCVertices;
                    delete reduction;
                    return INSUFFICIENT_BUDGET;
                }
                applied = true;
                //std::cout << cp::dye("vertex is unconfined", 'g') << std::endl;
                reduction->deletedVCVertices->push_back(vertex->getIndex());
                reduction->kDecrement++;
                (*k) = (*k) - 1;
                G->setInactive(vertex->getIndex());
                break;
            }
            if(bestOutsideNeighboursSize > 0)
            {
                //std::cout << cp::dye("continuing search", 'y') << std::endl;
                // remove u from neighbours and add it to S
                inSNeighbours[best] = false;
                inS[best] = true;
                for (int f=0; f<(int) SNeighbours.size(); f++)
                {
                    if(SNeighbours[f] == best) { SNeighbours.erase(SNeighbours.begin()+f); break; }
                }
                // push u's neighbours into S and add non-intersecting neighbourhood to neighbours
                for(auto it=G->getVertex(best)->adj->begin(); it != G->getVertex(best)->adj->end(); ++it)
                {
                    inSNeighbours[*it] = false;
                    inS[*it] = true;
                }
                for(auto it=G->getVertex(best)->adj->begin(); it != G->getVertex(best)->adj->end(); ++it)
                {
                    for(auto itt=G->getVertex(*it)->adj->begin(); itt != G->getVertex(*it)->adj->end(); ++itt)
                    {
                        if(!inSNeighbours[*itt] && !inS[*itt]) { inSNeighbours[*itt] = true; SNeighbours.push_back(*itt); }
                    }
                }
                continue;
            }
            //std::cout << cp::dye("vertex is not unconfined", 'r') << std::endl;
            break;
        }
    }
    //G->print();
    //std::cout << cp::dye("checking if applied", 'r') << std::endl;
    // if there were any applications return APPLICABLE
    if (applied) {
        //std::cout << cp::dye("was applied", 'g') << std::endl;
        appliedRules->push_back(reduction);
        return APPLICABLE;
    }
    //std::cout << cp::dye("wasn't applied", 'r') << std::endl;
    delete reduction->deletedVCVertices;
    delete reduction;
    return INAPPLICABLE;
} */

RULE_APPLICATION_RESULT Reductions::rule_Unconfined(BucketGraph* G, int* k, int depth, bool checkBudget, bool printDebug)
{
    auto startUnconfined = std::chrono::high_resolution_clock::now();
    int numCheckedVertices = 0;

    bool applied = false;

    std::unordered_map<int, bool> S = std::unordered_map<int, bool>();
    std::unordered_map<int, bool> SNeighbours = std::unordered_map<int, bool>();

    Reduction* reduction = new Reduction(RULE::DEGREE_ONE, 0, nullptr, new std::vector<int>());
    reduction->rDepth = depth;
    auto bucketQueue = G->getBucketQueue();
    //G->printBucketQueue();
    //for(auto b = G->getStableBucketQueueIncIterator(); (*b)->degree != bucketQueue->end()->degree; (*b) = ++(*b))
    for(auto b = G->getStableBucketQueueDecIterator(); (*b)->degree > 0 && (*b)->degree != bucketQueue->end()->degree; (*b) = --(*b))
    {
        //std::cout << "Degree: " << (*b)->degree << std::endl;
        Bucket* bucket = G->getBucket((*b)->degree);
        //std::cout << "Degree: " << bucket->degree << std::endl;
        if(bucket->degree <= 0) { continue; }
        //std::cout << "Number of vertices: " << (*bucket).vertices.size() << std::endl;
        for(auto v = bucket->getStableIterator(); (*v)->index != (*bucket).vertices.end()->index; (*v) = ++(*v))
        {
            //if(!G->isVertexScheduledForUnconfined((*v)->index)) { continue; }
            //std::cout << "Vertex: " << (*v)->index << std::endl;
            Vertex* vertex = G->getVertex((*v)->index);
            //std::cout << "Got vertex: " << (*v)->index << std::endl;
            //std::cout << cp::dye("checking if vertex " + std::to_string(vertex->getIndex()) + " is unconfined", 'y') << std::endl;
            if(vertex == nullptr || !vertex->getActive()) { continue; }
            numCheckedVertices++;

            //std::cout << cp::dye("vertex " + std::to_string(vertex->getIndex()) + " is active", 'g') << std::endl;
            // S and neighbours of S, that are kept up to date
            S.clear();
            SNeighbours.clear();
            S.insert(std::pair<int, bool>(vertex->getIndex(), true));
            for (auto neighbourIt = vertex->adj_map->begin(); neighbourIt != vertex->adj_map->end(); ++neighbourIt)
            {
                int neighbour = neighbourIt->first;
                if(!G->isActive(neighbour)) { continue; }
                SNeighbours.insert(std::pair<int, bool>(neighbour, true));
            }
            /* G->print();
            std::cout << cp::dye("finished initialisation with S={", 'b');
            for(auto l=S.begin(); l!=S.end(); ++l)
            {
                std::cout << cp::dye(std::to_string(l->first)+", ", 'y');
            }
            std::cout << cp::dye("}", 'b');
            std::cout << cp::dye(" and SN={", 'b');
            for(auto l=SNeighbours.begin(); l!=SNeighbours.end(); ++l)
            {
                std::cout << cp::dye(std::to_string(l->first)+", ", 'y');
            }
            std::cout << cp::dye("}", 'b') << std::endl; */

            // search continuation loop
            while(true)
            {
                int best = -1;
                int bestOutsideNeighbour = -1;
                // find best neighbour
                int off = 0;
                for (auto u = SNeighbours.begin(); u != SNeighbours.end(); ++u)
                {
                    /* if(!G->getVertex(u->first)->getActive()) { continue; }
                    if(S.find(u->first) != S.end()) { continue; } */
                    bool valid = true;
                    int sharedNeighbours = 0;
                    int outsideNeighbour = -1;
                    //int outsideNeighboursSize = 0;
                    for (auto it=G->getVertex(u->first)->adj_map->begin(); it != G->getVertex(u->first)->adj_map->end(); ++it)
                    {
                        int vertex = it->first;
                        if(!G->getVertex(vertex)->getActive()) { continue; }
                        // if u has more than one neighbour in S, u does not fit the criteria
                        if (S.find(vertex) != S.end()) {
                            sharedNeighbours++;
                            //std::cout << cp::dye("vertex "+std::to_string(u->first)+" has S neighbour: "+std::to_string(*it), 'y') << std::endl;
                            if(sharedNeighbours > 1) {
                                //std::cout << cp::dye("vertex shares multiple neighbours", 'r') << std::endl;
                                // TODO: remove vertex from neighbourhood
                                //SNeighbours.erase()
                                valid = false;
                                break;
                            }
                        }
                        //if(!inSNeighbours[*it] && !inS[*it])
                        if(SNeighbours.find(vertex) == SNeighbours.end() && S.find(vertex) == S.end())
                        {
                            //std::cout << cp::dye("vertex has neighbour: "+std::to_string(*it), 'y') << std::endl;
                            if(outsideNeighbour != -1) { valid = false; break; }
                            outsideNeighbour = vertex;
                        }
                    }

                    if(!valid) { continue; }
                    // update best/bestNeighbourhoodSize if u is better than the current best
                    //std::cout << cp::dye("Checking best=" + std::to_string(best) + " with o=" + std::to_string(bestOutsideNeighbour) + " vs. u=" + std::to_string(u->first) + " with o=" + std::to_string(outsideNeighbour), 'y') << std::endl;
                    //if(outsideNeighbour == -1) { break; }
                    if(bestOutsideNeighbour != -1 && outsideNeighbour != -1) { continue; }
                    bestOutsideNeighbour = outsideNeighbour;
                    best = u->first;
                    if(bestOutsideNeighbour != -1) { continue; }
                    break;
                }

                //std::cout << cp::dye("Best neighbour: " + std::to_string(best), 'y') << std::endl;

                // if no expansion vertex found, vertex "vertex" is not unconfined (continue with next vertex)
                if(best == -1) {
                    //std::cout << cp::dye("vertex is not unconfined", 'r') << std::endl;
                    break;
                }
                if(bestOutsideNeighbour == -1)
                {
                    //G->unscheduleForUnconfined(vertex->getIndex());     // TODO: test if doing it here or at the end for all deleted vertices performs better
                    if((*k) <= 0 && checkBudget) {
                        if (printDebug)
                        {
                            auto endUnconfined = std::chrono::high_resolution_clock::now();
                            double unconfinedDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endUnconfined - startUnconfined).count() /  1000) / (double) 1000;
                            std::cout << "#INSUFFICIENT BUDGET k=" << reduction->deletedVCVertices->size() << " for unconfined rule in " << unconfinedDuration << " seconds" << " (checked " << numCheckedVertices << " vertices)" << std::endl;
                        }
                        (*k) = (*k) + reduction->kDecrement;
                        G->setActive(reduction->deletedVCVertices);
                        delete reduction->deletedVCVertices;
                        delete reduction;
                        return INSUFFICIENT_BUDGET;
                    }
                    applied = true;
                    //std::cout << cp::dye("vertex is unconfined", 'g') << std::endl;
                    reduction->deletedVCVertices->push_back(vertex->getIndex());
                    reduction->kDecrement++;
                    (*k) = (*k) - 1;
                    G->setInactive(vertex->getIndex());
                    /* if(vertex->getIndex() == 33) {
                        G->print();
                        std::cout << "unconfined reducing Vertex " << 33 << std::endl;
                    } */
                    break;
                }
                if(bestOutsideNeighbour != -1)
                {
                    // remove u from neighbours and add it to S
                    //SNeighbours.erase(best);
                    //S.insert(std::pair<int, bool>(best, true));
                    // add u's outside neighbour to S
                    S.insert(std::pair<int, bool>(bestOutsideNeighbour, true));
                    // add non-intersecting outsideNeighbour neighbourhood to neighbours
                    for(auto it = G->getVertex(bestOutsideNeighbour)->adj_map->begin(); it != G->getVertex(bestOutsideNeighbour)->adj_map->end(); ++it)
                    {
                        int vertex = it->first;
                        if(!G->isActive(vertex) || S.find(vertex) != S.end()) { continue; }
                        if(SNeighbours.find(vertex) == SNeighbours.end())
                        {
                            SNeighbours.insert(std::pair<int, bool>(vertex, true));
                        }
                        else
                        {
                            SNeighbours.erase(vertex);
                        }
                    }
                    /* std::cout << cp::dye("continuing search with S={", 'y');
                    for(auto l=S.begin(); l!=S.end(); ++l)
                    {
                        std::cout << cp::dye(std::to_string(l->first)+", ", 'y');
                    }
                    std::cout << cp::dye("}", 'y');
                    std::cout << cp::dye(" and SN={", 'y');
                    for(auto l=SNeighbours.begin(); l!=SNeighbours.end(); ++l)
                    {
                        std::cout << cp::dye(std::to_string(l->first)+", ", 'y');
                    }
                    std::cout << cp::dye("}", 'y') << std::endl; */
                    continue;
                }
                //std::cout << cp::dye("vertex is not unconfined", 'r') << std::endl;
                break;
            }
        }
    }
    //G->print();
    //std::cout << cp::dye("checking if applied", 'r') << std::endl;
    auto endUnconfined = std::chrono::high_resolution_clock::now();
    double unconfinedDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endUnconfined - startUnconfined).count() /  1000) / (double) 1000;
    if (printDebug)
        std::cout << "#Reduced " << reduction->deletedVCVertices->size() << " unconfined vertices in " << unconfinedDuration << " seconds" << " (checked " << numCheckedVertices << " vertices)" << std::endl;
    // if there were any applications return APPLICABLE
    if (applied) {
        //std::cout << cp::dye("was applied", 'g') << std::endl;
        appliedRules->push_back(reduction);
        /* std::cout << cp::dye("reduced unconfined vertices U={", 'g');
        for(int l=0; l<(int)reduction->deletedVCVertices->size();l++)
        {
            std::cout << cp::dye(std::to_string(reduction->deletedVCVertices->at(l))+", ", 'g');
        }
        std::cout << cp::dye("}", 'g') << std::endl; */
        return APPLICABLE;
    }
    //std::cout << cp::dye("wasn't applied", 'r') << std::endl;
    delete reduction->deletedVCVertices;
    delete reduction;
    return INAPPLICABLE;
}

/* RULE_APPLICATION_RESULT Reductions::rule_Unconfined(BucketGraph* G, int* k, bool checkBudget)
{
    auto startUnconfined = std::chrono::high_resolution_clock::now();

    bool applied = false;
    std::vector<bool> inS = std::vector<bool>(G->getTotalNumVertices());
    std::vector<bool> inSNeighbours = std::vector<bool>(G->getTotalNumVertices());
    for(int i=0; i<(int) inS.size(); i++) { inS[i] = false; inSNeighbours[i] = false; }

    std::vector<int> SNeighbours = std::vector<int>();

    Reduction* reduction = new Reduction(RULE::DEGREE_ONE, 0, nullptr, new std::vector<int>());
    auto bucketQueue = G->getBucketQueue();
    //G->printBucketQueue();
    //for(auto b = G->getStableBucketQueueIncIterator(); (*b)->degree != bucketQueue->end()->degree; (*b) = ++(*b))
    for(auto b = G->getStableBucketQueueDecIterator(); (*b)->degree > 0 && (*b)->degree != bucketQueue->end()->degree; (*b) = --(*b))
    {
        //std::cout << "Degree: " << (*b)->degree << std::endl;
        Bucket* bucket = G->getBucket((*b)->degree);
        //std::cout << "Degree: " << bucket->degree << std::endl;
        if(bucket->degree <= 0) { continue; }
        //std::cout << "Number of vertices: " << (*bucket).vertices.size() << std::endl;
        //auto v = (*bucket).getStableIterator();
        //std::cout << "Beginning vertex: " << (*v)->index << std::endl;
        //std::cout << "Last vertex: " << (------(*bucket).vertices.end())->index << std::endl;
        for(auto v = bucket->getStableIterator(); (*v)->index != (*bucket).vertices.end()->index; (*v) = ++(*v))
        {
            //std::cout << "Vertex: " << (*v)->index << std::endl;
            //std::cout << cp::dye("getting vertex", 'y') << std::endl;
            Vertex* vertex = G->getVertex((*v)->index);
            //std::cout << "Got vertex: " << (*v)->index << std::endl;
            //std::cout << cp::dye("checking if vertex " + std::to_string(vertex->getIndex()) + " is unconfined", 'y') << std::endl;
            if(vertex == nullptr || !vertex->getActive()) { continue; }
            //std::cout << cp::dye("vertex " + std::to_string(vertex->getIndex()) + " is active", 'g') << std::endl;
            // S and neighbours of S list, that are kept up to date
            for(int i=0; i<(int) inS.size(); i++) { inS[i] = false; inSNeighbours[i] = false; }
            SNeighbours.clear();
            inS[vertex->getIndex()] = true;
            for (auto neighbour = vertex->adj->begin(); neighbour != vertex->adj->end(); ++neighbour)
            {
                if(!G->isActive(*neighbour)) { continue; }
                inSNeighbours[*neighbour] = true;
                SNeighbours.push_back(*neighbour);
            }
            //G->print();
            //std::cout << cp::dye("checking if vertex " + std::to_string(vertex->getIndex()) + " is unconfined", 'b') << std::endl;

            //int index = vertex->getIndex();
            // search continuation loop
            while(true)
            {
                int best = -1;
                int bestOutsideNeighboursSize = -1;
                // find best neighbour
                int off = 0;
                for (auto u = SNeighbours.begin(); u != SNeighbours.end(); ++u)
                //for (int i = 0; i < (int) SNeighbours.size(); ++i)
                {
                    if(!G->getVertex(*u)->getActive()) { continue; }
                    if(inS[*u]) { continue; }
                    bool valid = true;
                    int sharedNeighbours = 0;
                    int outsideNeighboursSize = 0;
                    for (auto it=G->getVertex(*u)->adj->begin(); it != G->getVertex(*u)->adj->end(); ++it)
                    {
                        if(!G->getVertex(*u)->getActive()) { continue; }
                        if(!inSNeighbours[*it] && !inS[*it])
                        {
                            outsideNeighboursSize++;
                        }
                        // if u has more than one neighbour in S, u does not fit the criteria
                        if (inS[*it]) {
                            sharedNeighbours++;
                            if(sharedNeighbours > 1) {
                                valid = false;
                                break;
                            }
                        }
                    }

                    if(!valid) { continue; }
                    // update best/bestNeighbourhoodSize if u is better than the current best
                    //std::cout << cp::dye("Checking best=" + std::to_string(best) + " with " + std::to_string(bestUniqueNeighboursSize) + " vs. u=" + std::to_string(*u) + " with " + std::to_string(uniqueNeighbours.size()), 'y') << std::endl;
                    if(best != -1 && outsideNeighboursSize >= bestOutsideNeighboursSize) { continue; }
                    best = *u;
                    bestOutsideNeighboursSize = outsideNeighboursSize;
                }

                //std::cout << cp::dye("Best neighbour: " + std::to_string(best), 'y') << std::endl;

                // if no expansion vertex found, vertex "vertex" is not unconfined (continue with next vertex)
                if(best == -1) { break; }
                if(bestOutsideNeighboursSize == 0)
                {
                    if((*k) <= 0 && checkBudget) {
                        (*k) = (*k) + reduction->kDecrement;
                        G->setActive(reduction->deletedVCVertices);
                        delete reduction->deletedVCVertices;
                        delete reduction;
                        return INSUFFICIENT_BUDGET;
                    }
                    applied = true;
                    //std::cout << cp::dye("vertex is unconfined", 'g') << std::endl;
                    reduction->deletedVCVertices->push_back(vertex->getIndex());
                    reduction->kDecrement++;
                    (*k) = (*k) - 1;
                    G->setInactive(vertex->getIndex());
                    break;
                }
                if(bestOutsideNeighboursSize > 0)
                {
                    //std::cout << cp::dye("continuing search", 'y') << std::endl;
                    // remove u from neighbours and add it to S
                    inSNeighbours[best] = false;
                    inS[best] = true;
                    for (int f=0; f<(int) SNeighbours.size(); f++)
                    {
                        if(SNeighbours[f] == best) { SNeighbours.erase(SNeighbours.begin()+f); break; }
                    }
                    // push u's neighbours into S and add non-intersecting neighbourhood to neighbours
                    for(auto it=G->getVertex(best)->adj->begin(); it != G->getVertex(best)->adj->end(); ++it)
                    {
                        inSNeighbours[*it] = false;
                        inS[*it] = true;
                    }
                    for(auto it=G->getVertex(best)->adj->begin(); it != G->getVertex(best)->adj->end(); ++it)
                    {
                        for(auto itt=G->getVertex(*it)->adj->begin(); itt != G->getVertex(*it)->adj->end(); ++itt)
                        {
                            if(!inSNeighbours[*itt] && !inS[*itt]) { inSNeighbours[*itt] = true; SNeighbours.push_back(*itt); }
                        }
                    }
                    continue;
                }
                //std::cout << cp::dye("vertex is not unconfined", 'r') << std::endl;
                break;
            }
        }
    }
    //G->print();
    //std::cout << cp::dye("checking if applied", 'r') << std::endl;
    auto endUnconfined = std::chrono::high_resolution_clock::now();
    double unconfinedDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endUnconfined - startUnconfined).count() /  1000) / (double) 1000;
    std::cout << "Reduced " << reduction->deletedVCVertices->size() << " unconfined vertices in " << unconfinedDuration << " seconds" << std::endl;
    // if there were any applications return APPLICABLE
    if (applied) {
        //std::cout << cp::dye("was applied", 'g') << std::endl;
        appliedRules->push_back(reduction);
        return APPLICABLE;
    }
    //std::cout << cp::dye("wasn't applied", 'r') << std::endl;
    delete reduction->deletedVCVertices;
    delete reduction;
    return INAPPLICABLE;
} */

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

RULE_APPLICATION_RESULT Reductions::rule_Domination(BucketGraph* G, int* k, int depth, bool checkBudget)
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
    reduction->rDepth = depth;


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

void Reductions::printReductionStack()
{
    if(appliedRules == nullptr) { std::cout << "< nullptr >\n"; return; }
    std::cout << "--------- Reduction Stack ---------\n";
    if(appliedRules->empty()) { std::cout << "< empty >\n"; return; }

    for(auto it = appliedRules->begin(); it != appliedRules->end(); ++it)
    {
        std::cout << "rule: " << (*it)->rule << ", depth: " << (*it)->rDepth << ", kDecrement: " << (*it)->kDecrement << "\n";
    }
}

void Reductions::freeReductions()
{
    //free appliedRules
    if(appliedRules)
    {
        for(int i = 0; i < (int) appliedRules->size(); ++i)
        {
            Reduction* reduction = appliedRules->at(i);
            freeReductionRule(reduction, true);
        }
        delete appliedRules;
        appliedRules = NULL;
    }

    if(dominationSets)
    {
        for(int i = 0; i < (int) dominationSets->size(); ++i)
        {
            auto entry = dominationSets->at(i);
            if(entry) 
            { 
                delete entry;
                entry = NULL;
            }
        }
        delete dominationSets;
        dominationSets = NULL;
    }
}

/**
 * set @freeMergeVertexInfoData to true if you want to also free saved adjacency lists
*/
void Reductions::freeReductionRule(Reduction* reduction, bool freeMergeVertexInfoData)
{
    if(reduction)
    {
        if(reduction->deletedVertices) { delete reduction->deletedVertices; reduction->deletedVertices = NULL; }
        if(reduction->deletedVCVertices) { delete reduction->deletedVCVertices; reduction->deletedVCVertices = NULL; }
        if(reduction->mergeVertexInfo) {
            std::tuple<int, std::unordered_map<int, bool>*, std::vector<int>*>* mergeVertexInfo = reduction->mergeVertexInfo;
            if(freeMergeVertexInfoData)
            {
                auto original_adj_map = std::get<1>(*mergeVertexInfo);
                if(original_adj_map) { delete original_adj_map; original_adj_map = NULL; }
                auto added_vertices = std::get<2>(*mergeVertexInfo);
                if(added_vertices) { delete added_vertices; added_vertices = NULL; }
            }
            delete reduction->mergeVertexInfo;
            reduction->mergeVertexInfo = NULL;
        }
        delete reduction;
        reduction = NULL;
    }
}
