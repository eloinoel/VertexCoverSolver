#ifndef REDUCTIONS_H
#define REDUCTIONS_H

#include <vector>

enum RULE{
    DEGREE_ZERO,        // = 0
    DEGREE_ONE,         // = 1
    DEGREE_TWO,         // = 2
    HIGH_DEGREE,        // = 3
    DOMINATION,         // = 4
    LPFLOW              // = 5
};

class BucketGraph;

class Reduction
{
public:
    RULE rule;
    int kDecrement;
    std::vector<int>* deletedVertices; // First idx is always to add in VC if(rule!=0)
    std::vector<int>* savedAdjacency;

    Reduction() {};
    Reduction(RULE rule) { this->rule = rule; };
    Reduction(RULE rule, int kDecrement, std::vector<int>* deletedVertices)
    {
        this->rule = rule;
        this->kDecrement = kDecrement;
        this->deletedVertices = deletedVertices;
    };
};

class Reductions
{
public:
    std::vector<Reduction*>* appliedRules;

public:
    //rules return true if they were applicable

    bool rule_HighDegree(BucketGraph* G, int* k);
    bool rule_DegreeZero(BucketGraph* G);
    /* only call if rule_HighDegree and rule_DegreeZero return false, returns true if no vertex cover of size k exists in graph */
    bool rule_Buss(BucketGraph* G, int* k, int numVertices, int numEdges);
    bool rule_DegreeOne(BucketGraph* G, int* k);

    //TODO:

    bool rule_DegreeTwo(int* k);

    bool rule_LPFlow(BucketGraph* G, int* k);
};

#endif