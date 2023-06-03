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

class Reduction{
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

    Reduction(RULE rule, int kDecrement, std::vector<int>* deletedVertices, std::vector<int>* saveAdjacency)
    {
        this->rule = rule;
        this->kDecrement = kDecrement;
        this->deletedVertices = deletedVertices;
        this->savedAdjacency = saveAdjacency;
    };
};

class Reductions
{
public:
    Reductions(BucketGraph* _G){G = _G;
        initRuleCounter();
    };

    BucketGraph* G;
    std::vector<Reduction*>* appliedRules;

    int rule_0;
    int rule_1;
    int rule_2;
    int rule_3;
    int rule_4;
    int rule_5;

//    bool printDebug;

public:
    // return bool indicating if no vertex cover possible
//    bool applyReductionRules(int* k, std::vector<ReductionVertices>* reductionArray, bool printDebug);

    // Adds the deleted vertices from the reduction rules to the vertex cover
    void addReducedVertices(std::vector<int>* S, bool printDebug);

    // Restores the initial kernel problem
    void addBackReducedVertices(int *k, bool printDebug);

    void initRuleCounter();

    void printReductionRules();

    bool rule_DegreeTwo(int* k);

    int size(){return (int)appliedRules->size();};

    //------------------------ Reduction Rules -----------------------
    //================================================================
    bool rule_HighDegree(int *k, bool printDebug);
    bool rule_DegreeZero(bool printDebug);

    bool rule_Buss(int* k, bool printDebug);

    void rule_DegreeOne(int* k, bool printDebug);
//    void rule_DegreeTwo(int* k, std::vector<ReductionVertices>* reductionArray, bool printDebug);
//    void rule_Domination(int* k, std::vector<ReductionVertices>* reductionArray, bool printDebug);
    //===============================================================

    /*
    // return bool indicating if no vertex cover possible
//    bool applyReductionRules(int* k, std::vector<ReductionVertices>* reductionArray, bool printDebug);

    // Adds the deleted vertices from the reduction rules to the vertex cover
    void addReducedVertices(std::vector<int>* S, std::vector<ReductionVertices>* reductionArray, bool printDebug);

    // Restores the initial kernel problem
    void addBackReducedVertices(int *k, std::vector<ReductionVertices>* reductionArray, bool printDebug);

    void initRuleCounter();

    void printReductionRules(std::vector<ReductionVertices>* reductionArray);

    bool rule_DegreeTwo(int* k);

    //------------------------ Reduction Rules -----------------------
    //================================================================
    bool rule_HighDegree(int *k, std::vector<ReductionVertices>* reductionVertices, bool printDebug);
    bool rule_DegreeZero(std::vector<ReductionVertices>* reductionArray, bool printDebug);

    bool rule_Buss(int* k, bool printDebug);

    void rule_DegreeOne(int* k, std::vector<ReductionVertices>* reductionArray, bool printDebug);
//    void rule_DegreeTwo(int* k, std::vector<ReductionVertices>* reductionArray, bool printDebug);
//    void rule_Domination(int* k, std::vector<ReductionVertices>* reductionArray, bool printDebug);
    *///================================================================
};

bool reduce(int* k, Reductions* reductions);

#endif
