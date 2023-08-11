#include <unordered_map>
#include <vector>

#include <chrono>
#include <stdio.h>
#include <iostream>

#include "Benchmark.h"


void benchmarkUnorderedMapVsVectorIterator()
{
    const long long int NUM_ELEMENTS = 100100100;
    const int NUM_RUNS = 100100;
    std::unordered_map<int, bool>* map = new std::unordered_map<int, bool>();
    std::vector<int>* vector = new std::vector<int>();

    //fill data structures
    for(int i = 0; i < NUM_ELEMENTS; i++)
    {
        map->insert({i, true});
        vector->push_back(i);
    }

    std::vector<double> unorderedMapDurations = std::vector<double>();
    std::vector<double> vectorDurations = std::vector<double>();
    //benchmark
    for(int i = 0; i < NUM_RUNS; ++i) 
    {
        auto startUnorderedMapIt = std::chrono::high_resolution_clock::now();
        for(auto it = map->begin(); it != map->end(); ++it)
        {

        }
        auto endUnorderedMapIt = std::chrono::high_resolution_clock::now();
        double unorderedMapItDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endUnorderedMapIt - startUnorderedMapIt).count());
        unorderedMapDurations.push_back(unorderedMapItDuration);

        auto startVectorIt = std::chrono::high_resolution_clock::now();
        for(auto it = vector->begin(); it != vector->end(); ++it)
        {

        }
        auto endVectorIt = std::chrono::high_resolution_clock::now();
        double vectorItDuration = (std::chrono::duration_cast<std::chrono::microseconds>(endVectorIt - startVectorIt).count());
        vectorDurations.push_back(vectorItDuration);
    }

    double meanUnorderedMapItDuration = 0.f;
    double meanVectorItDuration = 0.f;
    for(auto it = unorderedMapDurations.begin(); it != unorderedMapDurations.end(); ++it)
    {
        meanUnorderedMapItDuration += *it;
    }
    meanUnorderedMapItDuration = meanUnorderedMapItDuration / (double) unorderedMapDurations.size();
    for(auto it = vectorDurations.begin(); it != vectorDurations.end(); ++it)
    {
        meanVectorItDuration += *it;
    }
    meanVectorItDuration = meanVectorItDuration / (double) vectorDurations.size();

    std::cout << "Iterator benchmark: NUM_ELEMENTS = " << NUM_ELEMENTS << ", NUM_RUNS: " << NUM_RUNS << "\n";
    std::cout << "Mean unordered map duration: " << meanUnorderedMapItDuration << "\n";
    std::cout << "Mean vector duration: " << meanVectorItDuration << "\n";

    delete map;
    map = NULL;
    delete vector;
    vector = NULL;
}
