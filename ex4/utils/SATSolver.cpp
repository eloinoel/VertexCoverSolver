#include <stdexcept>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <array>
#include <unordered_map>

#include "SATSolver.h"

using namespace std;

string SATSolver::eraseLeadingTrailingWhitespacesFromString(string str)
{
    string whitespace = " \t";
    const auto strBegin = str.find_first_not_of(whitespace);	//filter leading spaces and tabs
    if (strBegin == string::npos) //no content
    {
        return "";
    }

    const auto strEnd = str.find_last_not_of(whitespace);

    return str.substr(strBegin, strEnd - strBegin + 1);
}

bool SATSolver::isVertexCharacter(char c)
{
    if ((c >= '0' && c <= '9') || (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == '_') {
        return true;
    }
    return false;
}

vector<pair<string, string>> SATSolver::readStandardInput()
{
    vertexIndex = 0;
    edgeCount = 0;
    vector<pair<string, string>> edges;  //edges after reading file once, eg. edges[0] = ["a", "b"]

    //iterate through file, map string name to unique index id and degree of vertex
    std::string line;
    while (getline(std::cin, line, '\n'))
    {
        //ignore empty lines
        if (line.empty())
        {
            continue;
        }

        //delete leading and trailing whitespaces
        line = eraseLeadingTrailingWhitespacesFromString(line); //not very fast, better to read byte by byte and ignore spaces

        //extract edge
        string vertex0 = "";
        string vertex1 = "";
        int i;
        bool foundComment = false;
        //first vertex
        for (i = 0; i < (int) line.size(); i++)
        {
            if (isVertexCharacter(line[i]))
            {
                vertex0 += line[i];
            }
                //whitespace
            else if (line[i] == ' ')
            {
                i++;
                break;
            }
            else if (line[i] == '#')
            {
                foundComment = true;
                break;
            }
            else
            {
                std::cerr << "readInput: illegal character read for vertex name 1\n";
                return vector<pair<string, string>>();
            }
        }

        //skip remainder of line if comment found
        if (foundComment)
            continue;

        //second vertex
        for (int j = i; j < (int) line.size(); j++)
        {
            if (isVertexCharacter(line[j]))
            {
                vertex1 += line[j];
            }
                //break if anything else
            else if (line[j] == '#')
            {
                break;
            }
                // fix for OS-side CRLF end of lines
            else if (line[j] == (char) 13)
            {
                continue;
            }
            else
            {
                std::cout << (int) line[j];
                std::cerr << "readInput: illegal character read for vertex name 2\n";
                return vector<pair<string, string>>();
            }
        }

        //add to hash map or update degree
        auto entry = originalVertexNames.find(vertex0);
        if (entry == originalVertexNames.end())
        {
            //vertex not in map
            originalVertexNames.insert({vertex0, vertexIndex});
            indexToNames[vertexIndex] = vertex0;
            namesList.push_back(vertex0);
            vertexIndex++;
        }

        entry = originalVertexNames.find(vertex1);
        if (entry == originalVertexNames.end())
        {
            //vertex not in map
            originalVertexNames.insert({vertex1, vertexIndex});
            namesList.push_back(vertex1);
            indexToNames[vertexIndex] = vertex1;
            vertexIndex++;
        }

        //save edges
        pair<string, string> edge_pair = pair<string, string>({vertex0, vertex1});
        edgeCount++;
        edges.push_back(edge_pair);
    }

    return edges;

}

string SATSolver::writeOpbMinCond()
{
    string res = "min: ";
    string varPre = "+1*x";

//    for (int i = 0; i < vertexIndex; ++i) {
    for (int i = 0; i < vertexIndex; ++i) {
        res += varPre + namesList.at(i);
//        res += varPre + to_string(i);
        if(i != vertexIndex-1)
            res += " ";
    }
    res += ";";

    return res;
}

//string SATSolver::writeIntObpConstraint(pair<int, int> p)
//{
//    string res;
//
//    string varPre = "+1*x";
//    string cond = ">= +1;";
//
//    res = varPre + to_string(p.first) + " ";
//    res += varPre + to_string(p.second) + " ";
//    res += cond;
//
//    return res;
//}

string SATSolver::writeStringObpConstraint(pair<string, string> p)
{
    string res;

    string varPre = "+1*x";
    string cond = ">= +1;";

    res = varPre + p.first + " ";
    res += varPre + p.second + " ";
    res += cond;

    return res;
}

void SATSolver::createOpbFile(string file_name, vector<pair<string, string>> edges)
{
    file_name += ".opb";
    ofstream out(file_name);

    out << writeOpbMinCond() << endl;

    for (int i = 0; i < edgeCount; ++i) {
        out << writeStringObpConstraint(edges.at(i)) << endl;
    }
    out.close();
}

void SATSolver::writeOutputSolutionToOutput(string output)
{
    stringstream ss(output);
    string s;
    vector<string> v;
    while ( getline( ss, s, ' ' ) ) {
        v.push_back(s);
    }
    if((int)v.size() == 1){
        cout << endl;
        return;
    }
    string sol;
    int vertexId;
    for (int i = 0; i < (int)v.size(); ++i) {
        sol = v.at(i);
        if(sol[0] == 'x') {
            sol.erase(0, 1);
            cout << sol << endl;
        }
    }
}

string SATSolver::getSolution(string outFile)
{
    std::ifstream inputFile(outFile); // Open the input file
    string solution;

    if (inputFile.is_open()) { // Check if the file was opened successfully
        std::string line;
        while (std::getline(inputFile, line)) {
            if(line[0] == 'v')
                solution = line;
            // Process each line
//            std::cout << line << std::endl;
        }

        inputFile.close(); // Close the file
    } else {
        std::cerr << "Failed to open the file." << std::endl;
        return "";
    }

    solution.erase(0, 2);
//    cout << solution << endl;
    return solution;
}

void SATSolver::solver() {

    bool printDebug = false;
    vector<pair<string, string>> edges;
    if(printDebug)
        std::cout << "SAT Solver: " << "Reading input" << std::endl;
    edges = readStandardInput();
    if(printDebug)
        std::cout << "SAT Solver: " << "Read input" << std::endl;

    if(edges.empty())
        return;
    if(printDebug)
        std::cout << "SAT Solver: " << "Creating OPB-File" << std::endl;
    string opbFileName = "solvers/input";
    createOpbFile(opbFileName, edges);
    if(printDebug)
        std::cout << "SAT Solver: " << "OPB created" << std::endl;

    // Add ending to file namae for execution
    opbFileName += ".opb ";

    // Open the output file
    string outFileName = "solvers/";
    outFileName += "output.txt";
    std::ofstream outputFile(outFileName);

    // Command as a string
    std::string command = "./solvers/uwrmaxsat ";
    command += opbFileName;
    command += "-of > " + outFileName;
    if(printDebug)
        std::cout << "SAT Solver: " << "Created execution command" << std::endl;


    if(printDebug)
        cout << command << endl;

    // Execute the command and redirect the output to the file
    int result = system(command.c_str());

    if (result == -1) {
        // Handle the error
        std::cerr << "Something went wrong when executing the solver!\n";
    } else {
        // Process the return value as needed
        if(printDebug)
            std::cout << "SAT Solver: " << "Solver was Executed properly" << std::endl;
    }

    if(printDebug)
        std::cout << "SAT Solver: " << "Getting the solution" << std::endl;
    string solverSolution = getSolution(outFileName);

    if(printDebug)
        cout << "The solver solution is: " << solverSolution << endl;

    if(printDebug)
        std::cout << "SAT Solver: " << "Writting the solution out!" << std::endl;

    writeOutputSolutionToOutput(solverSolution);
    // Close the output file
    outputFile.close();


}

