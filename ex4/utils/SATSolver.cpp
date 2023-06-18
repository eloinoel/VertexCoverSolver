#include <stdexcept>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>

#include <csignal>

//#include <unistd.h>
//#include <sys/types.h>
//#include <sys/wait.h>
//#include <csignal>
//#include <thread>
//#include <atomic>

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

    out << writeOpbMinCond() << "\n";

    for (int i = 0; i < edgeCount; ++i) {
        out << writeStringObpConstraint(edges.at(i)) << "\n";
    }
    out.close();
}

string SATSolver::getSolution(string outFile)
{
    std::ifstream inputFile(outFile); // Open the input file
    string solution = "\n";
    string interrupted = "c *** Interrupted ***";

    if (inputFile.is_open()) { // Check if the file was opened successfully
        std::string line;
        while (std::getline(inputFile, line)) {
            if(line.find(interrupted) != string::npos)
                return "-1";

            if(line[0] == 'v')
                solution = line;
        }

        inputFile.close(); // Close the file
    } else {
        std::cerr << "Failed to open the file.\n";
        return "";
    }

    solution.erase(0, 2);
    return solution;
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
//        cout  << "\n";
        return;
    }
    string sol;
    int vertexId;
    for (int i = 0; i < (int)v.size(); ++i) {
        sol = v.at(i);
        if(sol[0] != '-')
        {
            sol.erase(0, 1);
            cout << sol << "\n";
        }
    }
}

string SATSolver::solver() {

    bool printDebug = false;
    vector <pair<string, string>> edges;

    if (printDebug)
        std::cout << "SAT Solver: " << "Reading input\n";
    edges = readStandardInput();
    if (printDebug)
        std::cout << "SAT Solver: " << "Read input\n";

    if (edges.empty())
        return "";
    if (printDebug)
        std::cout << "SAT Solver: " << "Creating OPB-File\n";

    string opbFileName = "solvers/input";
    createOpbFile(opbFileName, edges);

    if (printDebug)
        std::cout << "SAT Solver: " << "OPB created\n";

    // Add ending to file namae for execution
    opbFileName += ".opb ";

    // Open the output file
    string outFileName = "solvers/";
    outFileName += "output.txt";

    // Command as a string
    std::string command = "./solvers/uwrmaxsat ";
    command += opbFileName;
    command += "-of > " + outFileName;

    std::ofstream outputFile(outFileName);

    if (printDebug)
        cout << command << "\n";

    // Execute the command and redirect the output to the file
    int result = system(command.c_str());

    if (result == -1) {
        // Handle the error
        std::cerr << "Something went wrong when executing the solver!\n";
    }

    outputFile.close();

    if (printDebug)
        cout << "SAT Solver: Getting the solution\n";

    string solverSolution = getSolution(outFileName);


    if (solverSolution == "-1") {
        raise(SIGINT);
        return "-1";
    }

    if (printDebug)
        cout << "The solver solution is: " << solverSolution << "\n";

    if (printDebug)
        cout << "SAT Solver: " << "Writting the solution out!\n";


    return solverSolution;

}
