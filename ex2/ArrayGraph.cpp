#include <iostream> //output streams
#include <fstream>  //ifstream file opening

#include "ArrayGraph.h";

/**
 * reads an standard input and creates a graph out of received data
 * assumes that no duplicate edges are present in the read data
*/
ArrayGraph* ArrayGraph::readStandardInput()
{
    //init
	ArrayGraph* G = new ArrayGraph();
    int vertexIndex = 0;
    int edgeCount = 0;
    G->originalVertexNames = unordered_map<string, pair<int, int>>();
    vector<string[2]> edges = vector<string[2]>(); //edges after reading file once, eg. edges[0] = ["a", "b"]
    
    //iterate through file, map string name to unique index id and degree of vertex
    string line;
    while (getline(cin, line, '\n'))
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
                cerr << "readInput: illegal character read for vertex name 1\n";
                return NULL;
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
            else
            {
                cerr << "readInput: illegal character read for vertex name 2\n";
                return NULL;
            }
        }

        //add to hash map or update degree
        auto entry = G->originalVertexNames.find(vertex0);
        if (entry != G->originalVertexNames.end())
        {
            //vertex in map
            entry->second.second++; //update degree
        }
        else
        {
            //vertex not in map
            G->originalVertexNames.insert({vertex0, pair<int, int>(vertexIndex, 0)});
            vertexIndex++;
        }

        entry = G->originalVertexNames.find(vertex1);
        if (entry != G->originalVertexNames.end())
        {
            //vertex in map
            entry->second.second++; //update degree
        }
        else
        {
            //vertex not in map
            G->originalVertexNames.insert({vertex1, pair<int, int>(vertexIndex, 0)});
            vertexIndex++;
        }

        //save edges
        edges[edgeCount][0] = vertex0;
        edges[edgeCount][1] = vertex1;
        edgeCount++;
    }

    //-----------------------------------------------------
    //generate adjacency list of fixed size, and add edges
    //-----------------------------------------------------

    G->adjacencyList = vector<vector<int>*>(G->originalVertexNames.size());
    //clear first
    for (int k = 0; k < G->adjacencyList.size(); k++)
    {
        G->adjacencyList[k] = nullptr;
    }
    //add edges
    for (int k = 0; k < edges.size(); k++)
    {
        //get vertices
        auto firstVertexEntry = G->originalVertexNames.find(edges[k][0]);
        auto secondVertexEntry = G->originalVertexNames.find(edges[k][1]);

        if (firstVertexEntry == G->originalVertexNames.end() || secondVertexEntry == G->originalVertexNames.end())
        {
            throw invalid_argument("readStandardInput: inconsistency: map of original names returned nothing");
        }


        int indexFirst = (*firstVertexEntry).second.first;
        int maxDegFirst = (*firstVertexEntry).second.second;
        int indexSecond = (*secondVertexEntry).second.first;
        int maxDegSecond = (*secondVertexEntry).second.second;

        //first vertex: insert into adjacency list
        if (G->adjacencyList[indexFirst] != nullptr) //edges list initialised
        {
            //find first index to insert adjacent vertex
            for (int insertIndex = 0; insertIndex < G->adjacencyList[indexFirst]->size(); insertIndex++)
            {
                if(G->adjacencyList[indexFirst]->at(insertIndex) != 0)
                {
                    (*G->adjacencyList[indexFirst])[insertIndex] = indexSecond;
                    break;
                }
            }
        }
        else //edges list not yet initialised
        {
            G->adjacencyList[indexFirst] = new vector<int>(maxDegFirst);
            //clear
            for (int l = 0; l < G->adjacencyList[indexFirst]->size(); l++)
            {
                (*G->adjacencyList[indexFirst])[l] = 0;
            }
            //insert
            (*G->adjacencyList[indexFirst])[0] = indexSecond;
        }

        //second vertex: insert into adjacency list
        if (G->adjacencyList[indexSecond] != nullptr) //edges list initialised
        {
            //find first index to insert adjacent vertex
            for (int insertIndex = 0; insertIndex < G->adjacencyList[indexSecond]->size(); insertIndex++)
            {
                if(G->adjacencyList[indexSecond]->at(insertIndex) != 0)
                {
                    (*G->adjacencyList[indexSecond])[insertIndex] = indexFirst;
                    break;
                }
            }
        }
        else //edges list not yet initialised
        {
            G->adjacencyList[indexFirst] = new vector<int>(maxDegSecond);
            //clear
            for (int l = 0; l < G->adjacencyList[indexSecond]->size(); l++)
            {
                (*G->adjacencyList[indexSecond])[l] = 0;
            }
            //insert
            (*G->adjacencyList[indexSecond])[0] = indexFirst;
        }
    }
    return G;
}

string ArrayGraph::eraseLeadingTrailingWhitespacesFromString(string str)
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

bool ArrayGraph::isVertexCharacter(char c)
{
	if ((c >= '0' && c <= '9') || (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == '_') {
		return true;
	}
	return false;
}