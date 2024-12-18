#ifndef REDUCTIONS_H
#define REDUCTIONS_H

#include <vector>
#include <list>
#include <unordered_map>
#include <iostream>

enum RULE
{
    DEGREE_ZERO,        // = 0
    DEGREE_ONE,         // = 1
    DEGREE_TWO,         // = 2
    HIGH_DEGREE,        // = 3
    DOMINATION,         // = 4
    UNCONFINED,         // = 5
    LPFLOW,             // = 6
    DEGREE_THREE_IND,    // = 7
    DEGREE_THREE_CLIQ,  // = 8
    DEGREE_THREE_DOM,  // = 9
    DEGREE_FOUR_CLIQUE,  // = 10
    DEGREE_FOUR_DOM  // = 11
};

enum RULE_APPLICATION_RESULT
{
    INAPPLICABLE,
    APPLICABLE,
    INSUFFICIENT_BUDGET //k doesn't allow for more vertex deletions -> no possible vertex cover of size k
};

class BucketGraph;

class Reduction
{
public:
    RULE rule;
    int kDecrement;
    int rDepth = -2; //-2 if not set, -1 if preprocessing, 0, 1, 2 ... otherwise if in branching
    std::vector<int>* deletedVertices = nullptr; // First idx is always to add in VC if(rule!=0)
    std::vector<int>* deletedVCVertices = nullptr;
    std::vector<std::vector<int>>* addedEdges = nullptr;
    /* mergeVertex, original adj_map, added vertices */
    std::tuple<int, std::unordered_map<int, bool>*, std::vector<int>*>* mergeVertexInfo = nullptr;

    Reduction() {};
    Reduction(RULE rule) { this->rule = rule; };
    Reduction(RULE rule, int kDecrement, std::vector<int>* deletedVertices)
    {
        this->rule = rule;
        this->kDecrement = kDecrement;
        this->deletedVertices = deletedVertices;
    };
    Reduction(RULE rule, int kDecrement, std::vector<int>* deletedVertices, std::vector<int>* deletedVCVertices)
    {
        this->rule = rule;
        this->kDecrement = kDecrement;
        this->deletedVertices = deletedVertices;
        this->deletedVCVertices = deletedVCVertices;
    };
    Reduction(RULE rule, int kDecrement, std::vector<int>* deletedVertices, std::vector<int>* deletedVCVertices, int recursionDepth)
    {
        this->rule = rule;
        this->kDecrement = kDecrement;
        this->deletedVertices = deletedVertices;
        this->deletedVCVertices = deletedVCVertices;
        this->rDepth = recursionDepth;
    };
};

class Reductions
{
public:

    int rule_0 , rule_1, rule_2 , rule_High, rule_Dom, rule_LPF, rule_B;

    bool dominationHeuristic = true;
    int cntDom = 0;

    std::string lineDelimiter = "---------\n";

    bool printDebug = false;
    bool printTimer = false;

    // max Degree to treat
    int arbitraryDegreeLimiter = 100;

    std::vector<Reduction*>* appliedRules;
    std::vector<std::vector<int>*>* dominationSets;

    Reductions()
    {
        appliedRules = new std::vector<Reduction*>();
        dominationSets = new std::vector<std::vector<int>*>();
    }

public:
    void printReductionStack();
    inline void printRule(Reduction* rule) { std::cout << "Rule: " <<  rule->rule << std::endl; }

    void freeReductions();
    void freeReductionRule(Reduction* reduction, bool freeMergeVertexInfoData);

    bool isDominated(BucketGraph* G, int dom, std::vector<bool>* pendingDeletions , bool printDebug);
    void initRuleCounter();

    void printReductionRules();
    void printCounters();
    void printDominationSets();

    // Degree 3
    bool isTouchable(BucketGraph* G, std::vector<int> * neighbours, std::unordered_map<int, int>& inactiveMap);

    RULE_APPLICATION_RESULT rule_DegreeThree_Independent(BucketGraph* G, int depth, bool printDebug = false);
    RULE_APPLICATION_RESULT rule_DegreeThree_Clique(BucketGraph* G, int* k, int depth, bool checkBudget, bool printDebug = false);
//    RULE_APPLICATION_RESULT rule_DegreeThree_Clique_Better(BucketGraph* G, int depth, bool printDebug = false);
    RULE_APPLICATION_RESULT rule_DegreeThree_Domination(BucketGraph* G, int* k, int depth, bool deg2inc, bool checkBudget, bool printDebug = false);

    RULE_APPLICATION_RESULT rule_DegreeFour_Clique(BucketGraph* G, int* k, int depth, bool checkBudget, bool printDebug = false);
    RULE_APPLICATION_RESULT rule_DegreeFour_Domination(BucketGraph* G, int* k, int depth, bool checkBudget, bool printDebug = false);

    //rules return true if they were applicable
    RULE_APPLICATION_RESULT rule_HighDegree(BucketGraph* G, int* k, int depth);
    /* only call if rule_HighDegree and rule_DegreeZero return false, returns true if no vertex cover of size k exists in graph */
    RULE_APPLICATION_RESULT rule_Buss(BucketGraph* G, int* k, int numVertices, int numEdges);
    RULE_APPLICATION_RESULT rule_DegreeOne(BucketGraph* G, int* k, int depth, bool checkBudget, bool printDebug = false);

    RULE_APPLICATION_RESULT rule_DegreeTwo(BucketGraph* G, int* k, int depth, bool checkBudget, bool printDebug = false);
    RULE_APPLICATION_RESULT rule_DegreeTwo_Simple_Case(BucketGraph* G, int* k, int depth, bool checkBudget, bool printDebug = false);
    RULE_APPLICATION_RESULT rule_DegreeTwo_Secure(BucketGraph* G, int* k);

    RULE_APPLICATION_RESULT rule_LPFlow(BucketGraph* G, int* k, int depth, bool checkBudget, bool printDebug = false);
    
    RULE_APPLICATION_RESULT rule_Domination(BucketGraph* G, int* k, int depth, bool checkBudget);
//    RULE_APPLICATION_RESULT rule_Domination(BucketGraph* G, int* k, int depth, bool checkBudget);

    RULE_APPLICATION_RESULT rule_Domination_BE(BucketGraph* G, int* k, int depth, bool checkBudget);

    RULE_APPLICATION_RESULT rule_Unconfined(BucketGraph* G, int* k, int depth, bool checkBudget, bool printDebug = false);
};

#endif
