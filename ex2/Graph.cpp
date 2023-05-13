#include <string>
#include <iostream> //output streams
#include <fstream>  //ifstream file opening

#include "Graph.h"


using namespace std;

//template <typename T>
//class Vertex {
//private:
//	string id;
//
//	public Vertex(T* data, string id) {
//		this->data = data;
//		this->id = id;
//	}
//};


void Graph::addVertex(string s) 
{
	//add if doesn't already contain
	if (adjacencyMap.find(s) == adjacencyMap.end())
	{
		this->adjacencyMap.insert({ s, {} });
	}
}

/*
* adds an edge to the graph by updating the adjacency map
* adds entries to both vertices' adjacency lists
* adds vertices if not already in graph
*/
void Graph::addEdge(string first, string second) 
{
	auto edgeSet0 = adjacencyMap.find(first);
	auto edgeSet1 = adjacencyMap.find(second);

	//found entry for first vertex
	if (edgeSet0 != adjacencyMap.end()) 
	{
		edgeSet0->second.push_back(second); //doesn't check for duplicates
	}
	//doesn't exist yet
	else 
	{
		adjacencyMap.insert({ first, {second} });
	}

	//found entry for second vertex
	if (edgeSet1 != adjacencyMap.end()) 
	{
		edgeSet1->second.push_back(first); //doesn't check for duplicates
	}
	//doesn't exist yet
	else 
	{
		adjacencyMap.insert({ second, {first} });
	}

}

/*
* deletes a vertex and all its edges from the graph
* also updates adjacency lists of other vertices
* returns a vector* of deleted vertices
* return nullptr if an error occured
*/
vector<string>* Graph::deleteVertex(string v)
{
	vector<string> deletedVertices;
	vector<string> edgesOfV;
	try
	{
		edgesOfV = adjacencyMap.at(v);
	}
	catch (const exception&)
	{
		cerr << "deleteVertex: no such vertex in graph";
		return nullptr;
	}

	//erase v from other vertices' edgelists
	try
	{
		for (string vertex : edgesOfV)
		{
			//remove all occurences
			vector<string> edgesOfVertex = adjacencyMap.at(vertex);
			edgesOfVertex.erase(std::remove_if(
				edgesOfVertex.begin(), edgesOfVertex.end(),
				[v](const string& x) {
					return x == v;
				}), edgesOfVertex.end());

			//delete Vertex if no edges
			vector<string> deletedVerticesTmp;
			if (edgesOfVertex.size() == 0)
			{
				deletedVerticesTmp = *deleteVertex(vertex); //returned vector should only contain one element
			}

			//append deleted vertices to list
			deletedVertices.insert(deletedVertices.end(), deletedVerticesTmp.begin(), deletedVerticesTmp.end());
		}
	}
	catch (const exception&)
	{
		cerr << "deleteVertex: pruning adjacency lists resulted in an error";
		return nullptr;
	}

	adjacencyMap.erase(v);
	deletedVertices.push_back(v);

	return &deletedVertices;
}

/* 
* only deletes one vertex and its edges from the adjacency map
* doesn't update edges of other vertices
*/
void Graph::deleteAdjacencyMapEntry(string v)
{
	adjacencyMap.erase(v);
}

/*
* copy an adjacency map entry from another graph to this graph
*/
void Graph::addAdjacencyMapEntry(Graph* graphToCopyFrom, string vertex)
{
	auto it = graphToCopyFrom->adjacencyMap.find(vertex);
	if (it != graphToCopyFrom->adjacencyMap.end())
	{
		vector<string> edgeCopy = it->second;
		adjacencyMap.insert({ vertex, edgeCopy });

		//print();
	}
	else
		cerr << "addAdjacencyMapEntry: Could not find entry for vertex " + vertex + " in graph to copy from";
}

/*
* returns first edge that has valid vertices that exist in the graph
* returns nullptr if no such edge exists
*/
pair<string, string>* Graph::getFirstValidEdge() 
{
	if (adjacencyMap.size() > 0)
	{
		for (auto it = adjacencyMap.cbegin(); it != adjacencyMap.cend(); ++it)
		{
			//iterate through edges
			for (auto connection : it->second)
			{
				auto search = adjacencyMap.find(connection);
				if (search != adjacencyMap.end())
				{
					return new pair<string, string> (it->first, connection);
				}
			}
		}
	}
	return nullptr;
}

/* 
* returns number of vertices in graph (this is also the minimal number of edges) 
*/
int Graph::getSize()
{
	return adjacencyMap.size();
}

/*
* prints the adjacency list of the graph
*/
void Graph::print() 
{
	if (adjacencyMap.size() > 0)
	{
		cout << "\n";
		for (auto it = adjacencyMap.begin(); it != adjacencyMap.end(); ++it)
		{	
			cout << (*it).first << ": ";
			cout << getVectorContentsString(&(it->second));
			cout << "\n";
		}
		cout << "\n";
	}
}

/*
* generate one string out of vector contents
*/
string Graph::getVectorContentsString(vector<string> *vec)
{
	if (vec->size() > 1)
	{
		string s = "";
		for (int i = 0; i < vec->size() - 1; i++)
		{
			s += vec->at(i) + ", ";
		}
		s += vec->at(vec->size() - 1);
		return s;
	}
	else if (vec->size() == 1)
		return vec->at(0);
	else
		return "<empty vector>";
	
}

Graph* Graph::copy()
{
	Graph* copy = new Graph();
	copy->adjacencyMap = adjacencyMap; //copy-assignment
	//copy->numEdges = numEdges;
	return copy;
}

/*
* create a graph from a text file
* format: 
* - Whitespaces at the beginning and end of every line are ignored.
* - Text after # until the end of the line is ignored (treated as comment).
* - Empty lines are ignored.
* - All other lines define an edge by stating start and end vertex of the edge separated
	by a whitespace. Names of vertices may contain an arbitrary number of letters,
	numbers, and underscores.
*/
Graph* Graph::readInputFromFile(string fileName) 
{
	ifstream file(fileName);
	Graph* G;
		
	//get file contents and create graph
	if (file.is_open())
	{
		G = new Graph();
		string line;
		while (getline(file, line))
		{
			//ignore empty lines
			if (line.empty())
			{
				continue;
			}

			//delete leading and trailing whitespaces
			line = eraseLeadingTrailingWhitespacesFromString(line); //not very fast, better to read byte by byte and ignore spaces
				
			//comments
			if (line[0] == '#')
			{
				continue;
			}

			//extract edge
			string vertex0 = "";
			string vertex1 = "";
			int i;
			//first vertex
			for (i = 0; i < line.size(); i++)
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
				else
				{
					cerr << "readInput: illegal character read for vertex name";
					return nullptr;
				}
			}
			//second vertex
			for (i; i < line.size(); i++)
			{
				if (isVertexCharacter(line[i]))
				{
					vertex1 += line[i];
				}
				//break if anything else
				else
				{
					cerr << "readInput: illegal character read for vertex name";
					break;
				}
			}

			//add edge to graph
			G->addEdge(vertex0, vertex1);
		}
		return G;
	}
	return nullptr;

}

Graph* Graph::readStandardInput()
{
	Graph *G = new Graph();

	string line;
	try
	{
		//get file contents and create graph
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
			for (i = 0; i < line.size(); i++)
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
			for (i; i < line.size(); i++)
			{
				if (isVertexCharacter(line[i]))
				{
					vertex1 += line[i];
				}
				//break if anything else
				else if (line[i] == '#')
				{
					break;
				}
				else
				{
					cerr << "readInput: illegal character read for vertex name 2\n";
					return NULL;
				}
			}

			//add edge to graph
			G->addEdge(vertex0, vertex1);
		}
	}
	catch (const exception& e)
	{
		cout << "exception";
	}
	
	return G;
}

string Graph::eraseLeadingTrailingWhitespacesFromString(string str)
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

bool Graph::isVertexCharacter(char c)
{
	if ((c >= '0' && c <= '9') || (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == '_') {
		return true;
	}
	return false;
}

