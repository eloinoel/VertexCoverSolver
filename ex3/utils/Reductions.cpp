#include "ColorPrint.h"
#include "Reductions.h"

#include "boost/intrusive/list.hpp"

//using namespace boost::intrusive;

/*bool Reductions::rule_HighDegree(BucketGraph* G, int* k)
{
    if(!(G->getMaxDegree() > *k)) return false; //cannot apply rule
    
    Reduction* reduction = new Reduction(RULE::HIGH_DEGREE, 0, new std::vector<int>());
    appliedRules->push_back(reduction);

    //delete vertices that have to be in the vertex cover
    while(G->getMaxDegree() > *k)
    {
        int maxDegVertex = G->getMaxDegreeVertex();
        reduction->kDecrement++;
        reduction->deletedVertices->push_back(maxDegVertex);
        k--;
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
        int neighbourToDelete = G->getVertex(it->index)->adj->front();
        reduction->deletedVertices->push_back(it->index);
        reduction->kDecrement++;
        *k--;
        G->setInactive(it->index);
    }
    return true;
}*/

/*----------------------------------------------------------*/
/*------------------   Reduction Rules   -------------------*/
/*----------------------------------------------------------*/
void reduce(int* k, Reductions* reductions)
{
    bool printDebug = true;
//    reductions->applyReductionRules(k, reductionArray, printDebug);
    reductions->initRuleCounter();

    if(printDebug){
        std::cout << std::endl;
        std::string befR = "Before Reduction Graph State";
        std::cout << ColorPrint::dye(befR, 'r') << std::endl ;
        printActiveList();
        printBucketQueue();
        std::cout << std::endl;
    }


    // if both rules are not applicable
    if(!reductions->rule_HighDegree(k, printDebug) && !reductions->rule_DegreeZero(printDebug))
        //if Buss rule == true => no vertex cover
        if(rule_Buss(k, printDebug))
            return false;


    // Apply rule 0, 1 & 2 at every deactivation?
    reductions->rule_DegreeOne(k, printDebug);
    reductions->rule_DegreeZero(printDebug);
//    reduction->rule_DegreeTwo(k, reductionArray);

    if (printDebug)
    {
        std::cout << std::endl;
        std::string kBef = "After all reductions";
        std::cout << ColorPrint::dye(kBef, 'r') << std::endl ;
        reductions->printReductionRules();
        std::cout << std::endl;
    }
    return true;

}

void BucketGraph::initRuleCounter()
{
    rule_0 = 0;
    rule_1 = 0;
    rule_2 = 0;
    rule_3 = 0;
    rule_4 = 0;
    rule_5 = 0;
}

void BucketGraph::printReductionRules()
{
    if(appliedRules->empty())
        return;

    for (auto i = 0; i < (int) reductionArray->size(); ++i) {
        Reduction reductionR = appliedRules->at(i);
        std::string toPrint =  "Applied Reduction Rule Nr " + std::to_string(reductionR.rule) + ":";
        std::cout << ColorPrint::dye(toPrint, 'r') << std::endl ;

        std::string deleteV =  "Deleted Vertices are: ";
        for (int j = 0; j < (int) reductionR.deletedVertices.size(); ++j) {
            deleteV += std::to_string(reductionR.deletedVertices.at(j)) + ", ";
        }
        std::cout << ColorPrint::dye(deleteV, 'c') << std::endl ;

        if(reductionR.rule == 2)
        {
            std::string savedAdj = "Saved Adjacency List: ";
            for (int k = 0; k < (int) reductionR.savedAdjacency->size(); ++k)
            {
                savedAdj += std::to_string(reductionR.savedAdjacency->at(k)) + ", ";
            }
            std::cout << ColorPrint::dye(savedAdj, 'c') << std::endl;
        }
        std::cout << std::endl;
    }

}

/*
bool BucketGraph::applyReductionRules(int* k, std::vector<ReductionVertices>* reductionArray, bool printDebug)
{
    initRuleCounter();

    if(printDebug){
        std::cout << std::endl;
        std::string befR = "Before Reduction Graph State";
        std::cout << ColorPrint::dye(befR, 'r') << std::endl ;
        printActiveList();
        printBucketQueue();
        std::cout << std::endl;
    }


    // if both rules are not applicable
    if(!rule_HighDegree(k, reductionArray, printDebug) && !rule_DegreeZero(reductionArray, printDebug))
        //if Buss rule == true => no vertex cover
        if(rule_Buss(k, printDebug))
            return false;


    // Apply rule 0, 1 & 2 at every deactivation?
    rule_DegreeOne(k, reductionArray, printDebug);
    rule_DegreeZero(reductionArray, printDebug);
//    rule_DegreeTwo(k, reductionArray);

    if (printDebug)
    {
        std::cout << std::endl;
        std::string kBef = "After all reductions";
        std::cout << ColorPrint::dye(kBef, 'r') << std::endl ;
        printReductionRules(reductionArray);
        std::cout << std::endl;
    }

    return true;
}
*/

bool BucketGraph::rule_HighDegree(int *k, bool printDebug)
{
    if (printDebug)
    {
        std::cout << std::endl;
        std::string toPrint =  "Arrived in Reduction Rule High Degree";
        std::cout << ColorPrint::dye(toPrint, 'r') << std::endl ;
        std::cout << std::endl;
    }

    int cnt = 0;

    int maxDegIdx = G->getMaxDegreeVertex();
    int maxDeg = G->getVertexDegree(maxDegIdx);

    if (maxDegIdx == -1 || maxDeg == -1)
        return false;

    while(maxDeg > (*k))
    {
        int tooHighToHandle = G->getFirstVertexOfDegree(maxDeg);

        while(tooHighToHandle != -1)
        {
            // save deleted vertex
            Reduction delVer;
            delVer.rule = 3;
            delVer.deletedVertices.push_back(tooHighToHandle);
            delVer.kDecrement = 1;
            delVer.savedAdjacency = getNeighbours(tooHighToHandle);

            appliedRules->push_back(delVer);

            G->setInactive(tooHighToHandle);
            rule_3++;
            (*k)--;

            tooHighToHandle = G->getFirstVertexOfDegree(maxDeg);
            cnt++;

            if (printDebug)
            {
                std::string appliedOn = "Rule High Degree on Vertex Nr: " + std::to_string(tooHighToHandle+1);
                std::cout << ColorPrint::dye(appliedOn, 'c') << std::endl ;
                G->printActiveList();
                G->printBucketQueue();
                std::cout << std::endl;
            }
        }
        maxDeg--;
    }

    if (printDebug)
    {
        std::string exitRule = "Rule High Degree applied " + std::to_string(rule_3) + " times.";
        std::cout << ColorPrint::dye(exitRule, 'y') << std::endl ;
        std::cout << std::endl;
    }

    if(cnt == 0)
        return false;
    return true;
}

bool BucketGraph::rule_DegreeZero(bool printDebug)
{
    if (printDebug)
    {
        std::cout << std::endl;
        std::string toPrint =  "Arrived in Reduction Rule Degree 0";
        std::cout << ColorPrint::dye(toPrint, 'r') << std::endl ;
        std::cout << std::endl;
    }

    int cnt = 0;

    int degreeZeroIdx = G->getFirstVertexOfDegree(0);

    while(degreeZeroIdx != -1)
    {
        // save deleted vertex
        Reduction delVer;
        delVer.rule = 0;
        delVer.kDecrement = 0;
        delVer.deletedVertices.push_back(degreeZeroIdx);

        appliedRules->push_back(delVer);

        G->setInactive(degreeZeroIdx);

        cnt++;
        rule_0++;

        degreeZeroIdx = G->getFirstVertexOfDegree(0);

        if (printDebug)
        {
            std::string appliedOn = "Rule Degree 0 on Vertex Nr: " + std::to_string(degreeZeroIdx+1);
            std::cout << ColorPrint::dye(appliedOn, 'c') << std::endl ;
            G->printActiveList();
            G->printBucketQueue();
            std::cout << std::endl;
        }


    }

    if (printDebug)
    {
        std::string exitRule = "Rule Degree Zero applied " + std::to_string(rule_0) + " times.";
        std::cout << ColorPrint::dye(exitRule, 'y') << std::endl ;
        std::cout << std::endl;
    }

    if(cnt == 0)
        return false;
    return true;
}

bool BucketGraph::rule_Buss(int* k, bool printDebug)
{
    if (printDebug)
    {
        std::cout << std::endl;
        std::string enterPrint =  "Just entered Buss Rule!";
        std::cout << ColorPrint::dye(enterPrint, 'r') << std::endl ;
        std::cout << std::endl;
    }

    int k_square = std::pow((*k), 2);

    // If |V| > k^2 + k || |E|>k^2 => no vertex cover
    if((G->getNumVertices() > k_square + (*k)) || (G->getNumEdges() > k_square))
    {
        if (printDebug)
        {
            std::string enterPrint =  "Applying Buss Rule";
            std::cout << ColorPrint::dye(enterPrint, 'c') << std::endl ;
            std::cout << std::endl;
        }
        return true;
    }
    if (printDebug)
    {
        std::string exitPrint =  "Buss Rule didn't apply";
        std::cout << ColorPrint::dye(exitPrint, 'c') << std::endl ;
        std::cout << std::endl;
    }

    return false;
}

void BucketGraph::rule_DegreeOne(int* k, bool printDebug)
{
    if (printDebug)
    {
        std::cout << std::endl;
        std::string enterPrint =  "Arrived in Reduction Rule Degree 1";
        std::cout << ColorPrint::dye(enterPrint, 'r') << std::endl ;
        std::cout << std::endl;
    }

    int degreeOneIdx = G->getFirstVertexOfDegree(1);

    while(degreeOneIdx != -1)
    {
        // THIS SHOULD HOLD THROUGH
        int neighbour = G->getNeighbours(degreeOneIdx)->at(0);

        // save deleted vertex
        Reduction delVer;
        delVer.rule = 1;
        delVer.kDecrement = 1;
        delVer.deletedVertices.push_back(neighbour);
        delVer.deletedVertices.push_back(degreeOneIdx);
        delVer.savedAdjacency = getNeighbours(neighbour);
//        ->push_back(getNeighbours(neighbour));
        appliedRules->push_back(delVer);

        G->setInactive(neighbour);
        G->setInactive(degreeOneIdx);
        (*k)--;

        if (printDebug)
        {
            std::string appliedOn = "Rule Degree 1 on Vertex Nr: " + std::to_string(degreeOneIdx+1);
            std::cout << ColorPrint::dye(appliedOn, 'c') << std::endl ;
            G->printActiveList();
            G->printBucketQueue();
            std::cout << std::endl;
        }

//        std::cout << "Applying reduction on Vertex Nr: " << degreeOneIdx+1 << std::endl ;
        rule_1++;
        degreeOneIdx = G->getFirstVertexOfDegree(1);
    }

    if (printDebug)
    {
        std::cout << std::endl;
        std::string exitPrint = "Rule Degree 1 applied " + std::to_string(rule_1) + " times.";
        std::cout << ColorPrint::dye(exitPrint, 'y') << std::endl ;
        std::cout << std::endl;
    }

}

/*
void BucketGraph::rule_DegreeTwo(int* k, std::vector<ReductionVertices>* reductionArray, bool printDebug)
{
 std::cout << "Arrived in Reduction Rule Degree 2" << std::endl;
    std::vector<int> degreeTwo = *getVerticesDegree(2);

    if(degreeTwo.empty())
        return;

    for (auto i = 0; i < (int) degreeTwo.size(); ++i) {
        int vertexID = degreeTwo.at(i);

        std::vector<int>* neighbours = getNeighbours(vertexID);
        int neighbourOne = neighbours->at(0);
        int neighbourTwo = neighbours->at(1);

        // Get shortest Neighbourhood to go through less to check for connection
        int shortestNeighbourhoodIdx;
        int otherNeighbourIdx;
        if(getVertexDegree(neighbourOne) > getVertexDegree(neighbourTwo) )
        {
            shortestNeighbourhoodIdx = neighbourTwo;
            otherNeighbourIdx = neighbourOne;
        }
        else
        {
            shortestNeighbourhoodIdx = neighbourOne;
            otherNeighbourIdx = neighbourTwo;
        }

        std::vector<int>* shortestNeighbourhood = getNeighbours(shortestNeighbourhoodIdx);
        std::vector<int>* otherNeighbour = getNeighbours(otherNeighbourIdx);

        // save deleted vertex
        ReductionVertices delVer;
        delVer.rule = 2;
        delVer.deletedVertices.push_back(otherNeighbourIdx);            // at (0)
        delVer.deletedVertices.push_back(shortestNeighbourhoodIdx);     // at (1)
        delVer.deletedVertices.push_back(vertexID);                     // at (2)

        delVer.savedAdjacency = getNeighbours(vertexID);
//        ->push_back(getNeighbours(vertexID));

        // CASE The neighbours know each other
        if(contains(shortestNeighbourhood,otherNeighbourIdx))
        {
            delVer.kDecrement = 2;
//            setInactive(otherNeighbourIdx);
//            setInactive(shortestNeighbourhoodIdx);
            setInactive(vertexID);
            (*k)-=2;
        }
        // CASE Neighbours don't know each other => setInactive or do that in addReducedVertices?
        else
        {
            delVer.kDecrement = 1;
            setVertexAdjacency(vertexID, putAdjacencyTogether(shortestNeighbourhood, otherNeighbour));
//            setInactive(otherNeighbourIdx);
//            setInactive(shortestNeighbourhoodIdx);
            (*k)--;
        }
        setInactive(otherNeighbourIdx);
        setInactive(shortestNeighbourhoodIdx);

        reductionArray->push_back(delVer);
    }
}

void BucketGraph::rule_Domination(int* k, std::vector<ReductionVertices>* reductionArray, bool printDebug)
{
    int maxDegIdx = getMaxDegreeVertex();
    int maxDeg = getVertexDegree(maxDegIdx);

    if (maxDegIdx == -1 || maxDeg == -1)
        return false;

    // for now degree 2
    while(maxDeg > 2) {
        std::vector<int> vertices = *getVerticesDegree(maxDeg);
        maxDeg--;
        if(vertices.empty())
            continue;

        for(auto v: vertices)
        {

        }

    }
}
*/

void BucketGraph::addReducedVertices(std::vector<int>* S, bool printDebug)
{
    if(reductionArray->empty())
        return;

    for (auto i = 0; i < (int) reductionArray->size(); ++i) {
        ReductionVertices curReduction = reductionArray->at(i);
        int curAddition = curReduction.deletedVertices.at(0);
        switch (curReduction.rule) {
            // Degree Zero Rule
            case 0:
                break;
                // Degree One Rule
            case 1:
                S->push_back(curAddition);
                break;
                // Degree Two Rule
            case 2:
                if (curReduction.kDecrement == 2)
                {
                    S->push_back(curReduction.deletedVertices.at(1));
                }
                S->push_back(curAddition);
                break;
                // High Degree Rule
            case 3:
                S->push_back(curAddition);
                break;
                // Domination Rule
            case 4:

                break;
            default:
                std::cout<< "There shouldn't be this rule";
                break;
        }
    }
}

void BucketGraph::addBackReducedVertices(int *k, bool printDebug)
{
    if(reductionArray->empty())
        return;

    for (auto i = 0; i < (int) reductionArray->size(); ++i) {
        ReductionVertices curReduction = reductionArray->at(i);
        int curVertex = curReduction.deletedVertices.at(0);
        switch (curReduction.rule) {
            // Degree Zero
            case 0:
                setActive(curVertex);
                break;
                // Degree One
            case 1:
                setActive(curVertex);
                setActive(curReduction.deletedVertices.at(1));
//                (*k)++;
                break;
            case 2:
                if (curReduction.kDecrement == 2)
                {
                    setActive(curReduction.deletedVertices.at(2));
                }
                setActive(curReduction.deletedVertices.at(1));
                setActive(curVertex);
                break;
                // High Degree Rule
            case 3:
//                setVertexAdjacencyBack(curVertex, curReduction.savedAdjacency);
                setActive(curVertex);
//                (*k)++;
                break;
            case 4:

                break;
            default:
                std::cout<< "There shouldn't be this rule";
                break;
        }
        (*k) += curReduction.kDecrement;
    }
}

/*
void BucketGraph::printReductionRules(std::vector<ReductionVertices>* reductionArray)
{
    if(reductionArray->empty())
        return;

    for (auto i = 0; i < (int) reductionArray->size(); ++i) {
        ReductionVertices reductionR = reductionArray->at(i);
        std::string toPrint =  "Applied Reduction Rule Nr " + std::to_string(reductionR.rule) + ":";
        std::cout << ColorPrint::dye(toPrint, 'r') << std::endl ;

        std::string deleteV =  "Deleted Vertices are: ";
        for (int j = 0; j < (int) reductionR.deletedVertices.size(); ++j) {
            deleteV += std::to_string(reductionR.deletedVertices.at(j)) + ", ";
        }
        std::cout << ColorPrint::dye(deleteV, 'c') << std::endl ;

        if(!reductionR.rule == 0)
        {
            std::string savedAdj = "Saved Adjacency List: ";
            for (int k = 0; k < (int) reductionR.savedAdjacency->size(); ++k)
            {
                savedAdj += std::to_string(reductionR.savedAdjacency->at(k)) + ", ";
            }
            std::cout << ColorPrint::dye(savedAdj, 'c') << std::endl;
        }
        std::cout << std::endl;
    }

}

/*
bool BucketGraph::applyReductionRules(int* k, std::vector<ReductionVertices>* reductionArray, bool printDebug)
{
    initRuleCounter();

    if(printDebug){
        std::cout << std::endl;
        std::string befR = "Before Reduction Graph State";
        std::cout << ColorPrint::dye(befR, 'r') << std::endl ;
        printActiveList();
        printBucketQueue();
        std::cout << std::endl;
    }


    // if both rules are not applicable
    if(!rule_HighDegree(k, reductionArray, printDebug) && !rule_DegreeZero(reductionArray, printDebug))
        //if Buss rule == true => no vertex cover
        if(rule_Buss(k, printDebug))
            return false;


    // Apply rule 0, 1 & 2 at every deactivation?
    rule_DegreeOne(k, reductionArray, printDebug);
    rule_DegreeZero(reductionArray, printDebug);
//    rule_DegreeTwo(k, reductionArray);

    if (printDebug)
    {
        std::cout << std::endl;
        std::string kBef = "After all reductions";
        std::cout << ColorPrint::dye(kBef, 'r') << std::endl ;
        printReductionRules(reductionArray);
        std::cout << std::endl;
    }

    return true;
}


bool BucketGraph::rule_HighDegree(int *k, std::vector<ReductionVertices>* reductionArray, bool printDebug)
{
    if (printDebug)
    {
        std::cout << std::endl;
        std::string toPrint =  "Arrived in Reduction Rule High Degree";
        std::cout << ColorPrint::dye(toPrint, 'r') << std::endl ;
        std::cout << std::endl;
    }

    int cnt = 0;

    int maxDegIdx = getMaxDegreeVertex();
    int maxDeg = getVertexDegree(maxDegIdx);

    if (maxDegIdx == -1 || maxDeg == -1)
        return false;

    while(maxDeg > (*k))
    {
        int tooHighToHandle = getFirstVertexOfDegree(maxDeg);

        while(tooHighToHandle != -1)
        {
            // save deleted vertex
            ReductionVertices delVer;
            delVer.rule = 3;
            delVer.deletedVertices.push_back(tooHighToHandle);
            delVer.kDecrement = 1;
            delVer.savedAdjacency = getNeighbours(tooHighToHandle);
            reductionArray->push_back(delVer);

            setInactive(tooHighToHandle);
            rule_3++;
            (*k)--;

            tooHighToHandle = getFirstVertexOfDegree(maxDeg);
            cnt++;

            if (printDebug)
            {
                std::string appliedOn = "Rule High Degree on Vertex Nr: " + std::to_string(tooHighToHandle+1);
                std::cout << ColorPrint::dye(appliedOn, 'c') << std::endl ;
                printActiveList();
                printBucketQueue();
                std::cout << std::endl;
            }
        }
        maxDeg--;
    }

    if (printDebug)
    {
        std::string exitRule = "Rule High Degree applied " + std::to_string(rule_3) + " times.";
        std::cout << ColorPrint::dye(exitRule, 'y') << std::endl ;
        std::cout << std::endl;
    }

    if(cnt == 0)
        return false;
    return true;
}

bool BucketGraph::rule_DegreeZero(std::vector<ReductionVertices>* reductionArray, bool printDebug)
{
    if (printDebug)
    {
        std::cout << std::endl;
        std::string toPrint =  "Arrived in Reduction Rule Degree 0";
        std::cout << ColorPrint::dye(toPrint, 'r') << std::endl ;
        std::cout << std::endl;
    }

    int cnt = 0;

    int degreeZeroIdx = getFirstVertexOfDegree(0);

    while(degreeZeroIdx != -1)
    {
        // save deleted vertex
        ReductionVertices delVer;
        delVer.rule = 0;
        delVer.kDecrement = 0;
        delVer.deletedVertices.push_back(degreeZeroIdx);
        reductionArray->push_back(delVer);

        setInactive(degreeZeroIdx);

        cnt++;
        rule_0++;

        degreeZeroIdx = getFirstVertexOfDegree(0);

        if (printDebug)
        {
            std::string appliedOn = "Rule Degree 0 on Vertex Nr: " + std::to_string(degreeZeroIdx+1);
            std::cout << ColorPrint::dye(appliedOn, 'c') << std::endl ;
            printActiveList();
            printBucketQueue();
            std::cout << std::endl;
        }


    }

    if (printDebug)
    {
        std::string exitRule = "Rule Degree Zero applied " + std::to_string(rule_0) + " times.";
        std::cout << ColorPrint::dye(exitRule, 'y') << std::endl ;
        std::cout << std::endl;
    }

    if(cnt == 0)
        return false;
    return true;
}

bool BucketGraph::rule_Buss(int* k, bool printDebug)
{
    if (printDebug)
    {
        std::cout << std::endl;
        std::string enterPrint =  "Just entered Buss Rule!";
        std::cout << ColorPrint::dye(enterPrint, 'r') << std::endl ;
        std::cout << std::endl;
    }

    int k_square = std::pow((*k), 2);

    // If |V| > k^2 + k || |E|>k^2 => no vertex cover
    if((getNumVertices() > k_square + (*k)) || (getNumEdges() > k_square))
    {
        if (printDebug)
        {
            std::string enterPrint =  "Applying Buss Rule";
            std::cout << ColorPrint::dye(enterPrint, 'c') << std::endl ;
            std::cout << std::endl;
        }
        return true;
    }
    if (printDebug)
    {
        std::string exitPrint =  "Buss Rule didn't apply";
        std::cout << ColorPrint::dye(exitPrint, 'c') << std::endl ;
        std::cout << std::endl;
    }

    return false;
}

void BucketGraph::rule_DegreeOne(int* k, std::vector<ReductionVertices>* reductionArray, bool printDebug)
{
    if (printDebug)
    {
        std::cout << std::endl;
        std::string enterPrint =  "Arrived in Reduction Rule Degree 1";
        std::cout << ColorPrint::dye(enterPrint, 'r') << std::endl ;
        std::cout << std::endl;
    }

    int degreeOneIdx = getFirstVertexOfDegree(1);

    while(degreeOneIdx != -1)
    {
        // THIS SHOULD HOLD THROUGH
        int neighbour = getNeighbours(degreeOneIdx)->at(0);

        // save deleted vertex
        ReductionVertices delVer;
        delVer.rule = 1;
        delVer.kDecrement = 1;
        delVer.deletedVertices.push_back(neighbour);
        delVer.deletedVertices.push_back(degreeOneIdx);
        delVer.savedAdjacency = getNeighbours(neighbour);
//        ->push_back(getNeighbours(neighbour));
        reductionArray->push_back(delVer);

        setInactive(neighbour);
        setInactive(degreeOneIdx);
        (*k)--;

        if (printDebug)
        {
            std::string appliedOn = "Rule Degree 1 on Vertex Nr: " + std::to_string(degreeOneIdx+1);
            std::cout << ColorPrint::dye(appliedOn, 'c') << std::endl ;
            printActiveList();
            printBucketQueue();
            std::cout << std::endl;
        }

//        std::cout << "Applying reduction on Vertex Nr: " << degreeOneIdx+1 << std::endl ;


        rule_1++;
        degreeOneIdx = getFirstVertexOfDegree(1);
    }

    if (printDebug)
    {
        std::cout << std::endl;
        std::string exitPrint = "Rule Degree 1 applied " + std::to_string(rule_1) + " times.";
        std::cout << ColorPrint::dye(exitPrint, 'y') << std::endl ;
        std::cout << std::endl;
    }

}

/*
void BucketGraph::rule_DegreeTwo(int* k, std::vector<ReductionVertices>* reductionArray, bool printDebug)
{
 std::cout << "Arrived in Reduction Rule Degree 2" << std::endl;
    std::vector<int> degreeTwo = *getVerticesDegree(2);

    if(degreeTwo.empty())
        return;

    for (auto i = 0; i < (int) degreeTwo.size(); ++i) {
        int vertexID = degreeTwo.at(i);

        std::vector<int>* neighbours = getNeighbours(vertexID);
        int neighbourOne = neighbours->at(0);
        int neighbourTwo = neighbours->at(1);

        // Get shortest Neighbourhood to go through less to check for connection
        int shortestNeighbourhoodIdx;
        int otherNeighbourIdx;
        if(getVertexDegree(neighbourOne) > getVertexDegree(neighbourTwo) )
        {
            shortestNeighbourhoodIdx = neighbourTwo;
            otherNeighbourIdx = neighbourOne;
        }
        else
        {
            shortestNeighbourhoodIdx = neighbourOne;
            otherNeighbourIdx = neighbourTwo;
        }

        std::vector<int>* shortestNeighbourhood = getNeighbours(shortestNeighbourhoodIdx);
        std::vector<int>* otherNeighbour = getNeighbours(otherNeighbourIdx);

        // save deleted vertex
        ReductionVertices delVer;
        delVer.rule = 2;
        delVer.deletedVertices.push_back(otherNeighbourIdx);            // at (0)
        delVer.deletedVertices.push_back(shortestNeighbourhoodIdx);     // at (1)
        delVer.deletedVertices.push_back(vertexID);                     // at (2)

        delVer.savedAdjacency = getNeighbours(vertexID);
//        ->push_back(getNeighbours(vertexID));

        // CASE The neighbours know each other
        if(contains(shortestNeighbourhood,otherNeighbourIdx))
        {
            delVer.kDecrement = 2;
//            setInactive(otherNeighbourIdx);
//            setInactive(shortestNeighbourhoodIdx);
            setInactive(vertexID);
            (*k)-=2;
        }
        // CASE Neighbours don't know each other => setInactive or do that in addReducedVertices?
        else
        {
            delVer.kDecrement = 1;
            setVertexAdjacency(vertexID, putAdjacencyTogether(shortestNeighbourhood, otherNeighbour));
//            setInactive(otherNeighbourIdx);
//            setInactive(shortestNeighbourhoodIdx);
            (*k)--;
        }
        setInactive(otherNeighbourIdx);
        setInactive(shortestNeighbourhoodIdx);

        reductionArray->push_back(delVer);
    }
}

void BucketGraph::rule_Domination(int* k, std::vector<ReductionVertices>* reductionArray, bool printDebug)
{
    int maxDegIdx = getMaxDegreeVertex();
    int maxDeg = getVertexDegree(maxDegIdx);

    if (maxDegIdx == -1 || maxDeg == -1)
        return false;

    // for now degree 2
    while(maxDeg > 2) {
        std::vector<int> vertices = *getVerticesDegree(maxDeg);
        maxDeg--;
        if(vertices.empty())
            continue;

        for(auto v: vertices)
        {

        }

    }
}


void BucketGraph::addReducedVertices(std::vector<int>* S, std::vector<ReductionVertices>* reductionArray, bool printDebug)
{
    if(reductionArray->empty())
        return;

    for (auto i = 0; i < (int) reductionArray->size(); ++i) {
        ReductionVertices curReduction = reductionArray->at(i);
        int curAddition = curReduction.deletedVertices.at(0);
        switch (curReduction.rule) {
            // Degree Zero Rule
            case 0:
                break;
                // Degree One Rule
            case 1:
                S->push_back(curAddition);
                break;
                // Degree Two Rule
            case 2:
                if (curReduction.kDecrement == 2)
                {
                    S->push_back(curReduction.deletedVertices.at(1));
                }
                S->push_back(curAddition);
                break;
                // High Degree Rule
            case 3:
                S->push_back(curAddition);
                break;
                // Domination Rule
            case 4:

                break;
            default:
                std::cout<< "There shouldn't be this rule";
                break;
        }
    }
}

void BucketGraph::addBackReducedVertices(int *k, std::vector<ReductionVertices>* reductionArray, bool printDebug)
{
    if(reductionArray->empty())
        return;

    for (auto i = 0; i < (int) reductionArray->size(); ++i) {
        ReductionVertices curReduction = reductionArray->at(i);
        int curVertex = curReduction.deletedVertices.at(0);
        switch (curReduction.rule) {
            // Degree Zero
            case 0:
                setActive(curVertex);
                break;
                // Degree One
            case 1:
                setActive(curVertex);
                setActive(curReduction.deletedVertices.at(1));
//                (*k)++;
                break;
            case 2:
                if (curReduction.kDecrement == 2)
                {
                    setActive(curReduction.deletedVertices.at(2));
                }
                setActive(curReduction.deletedVertices.at(1));
                setActive(curVertex);
                break;
                // High Degree Rule
            case 3:
//                setVertexAdjacencyBack(curVertex, curReduction.savedAdjacency);
                setActive(curVertex);
//                (*k)++;
                break;
            case 4:

                break;
            default:
                std::cout<< "There shouldn't be this rule";
                break;
        }
        (*k) += curReduction.kDecrement;
    }
}
*/