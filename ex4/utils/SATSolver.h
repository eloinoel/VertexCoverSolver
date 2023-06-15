#ifndef VERTEXCOVERSOLVER_SATSOLVER_H
#define VERTEXCOVERSOLVER_SATSOLVER_H

#include <stdexcept>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <array>

using namespace std;


class SATSolver {
public:
    int vertexIndex;
    int edgeCount;
    unordered_map<string, int> originalVertexNames;
    vector<string> namesList;
    unordered_map<int, string> indexToNames;
private:

public:
    SATSolver(){solver();};
private:
    string eraseLeadingTrailingWhitespacesFromString(string str);

    bool isVertexCharacter(char c);

    vector<pair<string, string>> readStandardInput();

    string writeOpbMinCond();

    string writeIntObpConstraint(pair<int, int> p);

    string writeStringObpConstraint(pair<string, string> p);

    void createOpbFile(string file_name, vector<pair<string, string>> edges);

    void writeOutputSolutionToOutput(string output);

    string getSolution(string outFile);

    void solver();
};


#endif //VERTEXCOVERSOLVER_SATSOLVER_H
