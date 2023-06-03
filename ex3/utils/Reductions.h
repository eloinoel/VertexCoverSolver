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
    LPFLOW              // = 5
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
    std::vector<int>* deletedVertices; // First idx is always to add in VC if(rule!=0)
    std::vector<int>* deletedVCVertices;
    std::tuple<int, std::vector<int>*, std::unordered_map<int, bool>*>* mergeVertexInfo;

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
    //int rule_0 , rule_1, rule_2 , rule_3,  rule_4, rule_5;

    std::vector<Reduction*>* appliedRules;

    Reductions()
    {
        appliedRules = new std::vector<Reduction*>();
    }

public:
//    void initRuleCounter();

    void printReductionRules();

    //rules return true if they were applicable
    RULE_APPLICATION_RESULT rule_HighDegree(BucketGraph* G, int* k);
    RULE_APPLICATION_RESULT rule_DegreeZero(BucketGraph* G);
    /* only call if rule_HighDegree and rule_DegreeZero return false, returns true if no vertex cover of size k exists in graph */
    RULE_APPLICATION_RESULT rule_Buss(BucketGraph* G, int* k, int numVertices, int numEdges);
    RULE_APPLICATION_RESULT rule_DegreeOne(BucketGraph* G, int* k);

    RULE_APPLICATION_RESULT rule_DegreeTwo(BucketGraph* G, int* k);

    RULE_APPLICATION_RESULT rule_LPFlow(BucketGraph* G, int* k);
};

#endif