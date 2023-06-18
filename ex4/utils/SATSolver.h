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

class BucketGraph;

class SATSolver {
public:
    int vertexIndex;
    int edgeCount;
    unordered_map<string, int> originalVertexNames;
    vector<string> namesList;
    unordered_map<int, string> indexToNames;
private:

public:
    SATSolver(){};
//    SATSolver(){solver();};

    string solver(vector<pair<string, string>> edges);
//    string solver(vector<pair<string, string>> edges);

    int writeOutputSolutionToOutput(string output);

    void writeSolverSolutionToVC(BucketGraph* G, unordered_map<int, bool>* vc, string output);

    vector<pair<string, string>> readStandardInput();

    void setSATSolverUp(vector<pair<string, string>> edges);

    void setNumEdges(int numEdges){edgeCount = numEdges;};
    void setNumVertices(int numVertices){vertexIndex = numVertices;};
//    vector<pair<string, string>> getPreprocessedInput();
private:

    string eraseLeadingTrailingWhitespacesFromString(string str);

    bool isVertexCharacter(char c);

    string writeOpbMinCond();

    string writeIntObpConstraint(pair<int, int> p);

    string writeStringObpConstraint(pair<string, string> p);

    void createOpbFile(string file_name, vector<pair<string, string>> edges);


    string getSolution(string outFile);


};


#endif //VERTEXCOVERSOLVER_SATSOLVER_H
