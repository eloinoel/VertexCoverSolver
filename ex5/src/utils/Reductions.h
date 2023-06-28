#ifndef REDUCTIONS_H
#define REDUCTIONS_H

#include <vector>
#include <unordered_map>

enum RULE
{
    DEGREE_ZERO,        // = 0
    DEGREE_ONE,         // = 1
    DEGREE_TWO,         // = 2
    HIGH_DEGREE,        // = 3
    DOMINATION,         // = 4
    UNCONFINED,         // = 5
    LPFLOW              // = 6
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
    std::vector<int>* deletedVertices = nullptr; // First idx is always to add in VC if(rule!=0)
    std::vector<int>* deletedVCVertices = nullptr;
    /* mergeVertex, original adj, original adj_map, added vertices */
    std::tuple<int, std::vector<int>*, std::unordered_map<int, bool>*, std::vector<int>*>* mergeVertexInfo = nullptr;

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
};

class Reductions
{
public:
    /* Domination Rule */
    bool dominationHeuristic = true;
    int cntDom = 0;
    bool printDebug = false;
    bool printTimer = false;
    int arbitraryDegreeLimiter = 100;
    /* End Domination Rule */

    std::vector<Reduction*>* appliedRules;

    Reductions()
    {
        appliedRules = new std::vector<Reduction*>();
    }

public:

    void freeReductions();
    void freeReductionRule(Reduction* reduction, bool freeMergeVertexInfoData);

    bool isDominated(BucketGraph* G, int dom, std::vector<bool>* pendingDeletions , bool printDebug);
    void initRuleCounter();

    void printReductionRules();
    void printCounters();
    void printDominationSets();

    RULE_APPLICATION_RESULT rule_HighDegree(BucketGraph* G, int* k);
    /* only call if rule_HighDegree and rule_DegreeZero return false, returns true if no vertex cover of size k exists in graph */
    RULE_APPLICATION_RESULT rule_Buss(BucketGraph* G, int* k, int numVertices, int numEdges);
    RULE_APPLICATION_RESULT rule_DegreeOne(BucketGraph* G, int* k, bool checkBudget);

    RULE_APPLICATION_RESULT rule_DegreeTwo(BucketGraph* G, int* k, bool checkBudget);
    RULE_APPLICATION_RESULT rule_DegreeTwo_Secure(BucketGraph* G, int* k);

    RULE_APPLICATION_RESULT rule_LPFlow(BucketGraph* G, int* k, bool checkBudget);
    
    RULE_APPLICATION_RESULT rule_Domination(BucketGraph* G, int* k, bool checkBudget);

    RULE_APPLICATION_RESULT rule_Domination_BE(BucketGraph* G, int* k, bool checkBudget);

    RULE_APPLICATION_RESULT rule_Unconfined(BucketGraph* G, int* k, bool checkBudget);
};

#endif
