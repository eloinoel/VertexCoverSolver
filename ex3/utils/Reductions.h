#ifndef REDUCTIONS_H
#define REDUCTIONS_H

#include <vector>

#include "BucketGraph.h"

enum RULE{
    DEGREE_ZERO,        // = 0
    DEGREE_ONE,         // = 1
    DEGREE_TWO,         // = 2
    HIGH_DEGREE,        // = 3
    DOMINATION          // = 4
};

class ReductionVertices{
public:
    RULE rule;
    int kDecrement;
    std::vector<int>* deletedVertices; // First idx is always to add in VC if(rule!=0)
    std::vector<int>* savedAdjacency;

    ReductionVertices() {};
    ReductionVertices(RULE rule) { this->rule = rule; };
    ReductionVertices(RULE rule, int kDecrement, std::vector<int>* deletedVertices)
    {
        this->rule = rule;
        this->kDecrement = kDecrement;
        this->deletedVertices = deletedVertices;
    };
};

class Reductions
{
public:
    std::vector<ReductionVertices*>* appliedRules;

public:
    bool rule_HighDegree(BucketGraph* G, int k);
    bool rule_DegreeZero();

    bool rule_Buss(int* k);

    void rule_DegreeOne(int* k);
    void rule_DegreeTwo(int* k);
};

#endif