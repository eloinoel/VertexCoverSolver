#include <iostream> //output streams
#include <fstream>  //ifstream file opening
#include <cmath> //ceil 

#include "ArrayGraph.h"

// TODO: make more efficent, employ less copies
/**
 * reads an standard input and creates a graph out of received data
 * assumes that no duplicate edges are present in the read data
*/
ArrayGraph* ArrayGraph::readStandardInput()
{
    //init
	ArrayGraph* G = new ArrayGraph();
    int vertexIndex = 0;
    G->originalVertexNames = std::unordered_map<std::string, std::pair<int, int>>();
    std::vector<std::pair<std::string, std::string>> edges = std::vector<std::pair<std::string, std::string>>(); //edges after reading file once, eg. edges[0] = ["a", "b"]
    
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
        std::string vertex0 = "";
        std::string vertex1 = "";
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
            // fix for OS-side CRLF end of lines
            else if (line[j] == (char) 13)
            {
                continue;
            }
            else
            {
                std::cout << (int) line[j];
                std::cerr << "readInput: illegal character read for vertex name 2\n";
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
            G->originalVertexNames.insert({vertex0, std::pair<int, int>(vertexIndex, 1)});
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
            G->originalVertexNames.insert({vertex1, std::pair<int, int>(vertexIndex, 1)});
            vertexIndex++;
        }

        //save edges
        std::pair<std::string, std::string> edge_pair = std::pair<std::string, std::string>({vertex0, vertex1});
        edges.push_back(edge_pair);
    }
    
    //-----------------------------------------------------
    //generate adjacency list of fixed size, and add edges
    //-----------------------------------------------------

    G->adjacencyList = std::vector<std::vector<int>*>(G->originalVertexNames.size());
    //clear first
    for (int k = 0; k < (int) G->adjacencyList.size(); k++)
    {
        G->adjacencyList[k] = nullptr;
    }
    
    //add edges
    for (int k = 0; k < (int) edges.size(); k++)
    {
        //get vertices
        auto firstVertexEntry = G->originalVertexNames.find(edges[k].first);
        auto secondVertexEntry = G->originalVertexNames.find(edges[k].second);

        if (firstVertexEntry == G->originalVertexNames.end() || secondVertexEntry == G->originalVertexNames.end())
        {
            throw std::invalid_argument("readStandardInput: inconsistency: map of original names returned nothing");
        }


        int indexFirst = (*firstVertexEntry).second.first;
        int maxDegFirst = (*firstVertexEntry).second.second;
        int indexSecond = (*secondVertexEntry).second.first;
        int maxDegSecond = (*secondVertexEntry).second.second;

        //first vertex: insert into adjacency list
        if (G->adjacencyList[indexFirst] != nullptr) //edges list initialised
        {
            //find first index to insert adjacent vertex
            for (int insertIndex = 0; insertIndex < (int) G->adjacencyList[indexFirst]->size(); insertIndex++)
            {
                if(G->adjacencyList[indexFirst]->at(insertIndex) == 0)
                {
                    //G->adjacencyList[indexFirst]->insert(G->adjacencyList[indexFirst]->begin() + insertIndex, indexSecond);
                    (*G->adjacencyList[indexFirst])[insertIndex] = indexSecond;
                    break;
                }
            }
        }
        else //edges list not yet initialised
        {
            G->adjacencyList[indexFirst] = new std::vector<int>(maxDegFirst);
            //clear
            for (int l = 0; l < (int) G->adjacencyList[indexFirst]->size(); l++)
            {
                //G->adjacencyList[indexFirst]->insert(G->adjacencyList[indexFirst]->begin() + l, 0);
                (*G->adjacencyList[indexFirst])[l] = 0;
            }
            //insert
            if (maxDegFirst > 0)
            {
                (*G->adjacencyList[indexFirst])[0] = indexSecond;
                //G->adjacencyList[indexFirst]->insert(G->adjacencyList[indexFirst]->begin(), indexSecond);
            }
            else
            {
                throw std::invalid_argument("Trying to insert edge but maxDegFirst is 0");
            }
        }

        //second vertex: insert into adjacency list
        if (G->adjacencyList[indexSecond] != nullptr) //edges list initialised
        {
            //find first index to insert adjacent vertex
            for (int insertIndex = 0; insertIndex < (int) G->adjacencyList[indexSecond]->size(); insertIndex++)
            {
                if(G->adjacencyList[indexSecond]->at(insertIndex) == 0)
                {
                    //G->adjacencyList[indexSecond]->insert(G->adjacencyList[indexSecond]->begin() + insertIndex, indexFirst);
                    (*G->adjacencyList[indexSecond])[insertIndex] = indexFirst;
                    break;
                }
            }
        }
        else //edges list not yet initialised
        {
            G->adjacencyList[indexSecond] = new std::vector<int>(maxDegSecond);
            //clear
            for (int l = 0; l < (int) G->adjacencyList[indexSecond]->size(); l++)
            {
                //G->adjacencyList[indexSecond]->insert(G->adjacencyList[indexSecond]->begin() + l, 0);
                (*G->adjacencyList[indexSecond])[l] = 0;
            }
            //insert
            if (maxDegFirst > 0)
            {
                //G->adjacencyList[indexSecond]->insert(G->adjacencyList[indexSecond]->begin(), indexFirst);
                (*G->adjacencyList[indexSecond])[0] = indexFirst;
            }
            else
            {
                throw std::invalid_argument("Trying to insert edge but maxDegFirst is 0");
            }
        }
    }

    G->initGraphState(G->originalVertexNames.size(), edges.size());
    return G;
}

std::string ArrayGraph::eraseLeadingTrailingWhitespacesFromString(std::string str)
{
	std::string whitespace = " \t";
	const auto strBegin = str.find_first_not_of(whitespace);	//filter leading spaces and tabs
	if (strBegin == std::string::npos) //no content
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

void ArrayGraph::print()
{
    if (adjacencyList.size() > 0)
	{
		std::cout << "\n";
		for (int i = 0; i < (int) adjacencyList.size(); i++)
		{	
			if (adjacencyList[i] != nullptr)
            {
                std::cout << i << ": ";
                for (int j = 0; j < (int) adjacencyList[i]->size() - 1; j++)
                {
                    std::cout << adjacencyList[i]->at(j) << ", ";
                }
                if (adjacencyList[i]->size() > 0)
                {
                    std::cout << adjacencyList[i]->at(adjacencyList[i]->size() - 1);
                }
            }
            else
            {
                std::cout << "Index " << i << " is nullptr";
            }
			std::cout << "\n";
		}
		std::cout << "\n";
	}
}

void ArrayGraph::printOriginalVertexNames()
{
    for (auto entry : originalVertexNames)
    {
        std::cout << entry.first << ": index = " << entry.second.first << ", degree = " << entry.second.second << std::endl;
    }
}

//TODO: implement cycle and clique bound
int ArrayGraph::getLowerBoundVC() {

    //int cycleBound = getCycleBound();
    //return cycleBound;
    return 0;
}

void ArrayGraph::initGraphState(int vertexCount, int edgeCount)
{
    numberOfVertices = vertexCount;
    numberOfEdges = edgeCount;
    graphState = new std::vector<std::pair<bool, int>>(vertexCount);
    for (auto entry : originalVertexNames)
    {
        graphState->at(entry.second.first) = std::make_pair(true, entry.second.second);
    }
}

std::vector<int>* ArrayGraph::getInactiveVertices()
{
    std::vector<int>* inactive = new std::vector<int>();
    for(int i = 0; i < (int) graphState->size(); i++)
    {
        if (!(graphState->at(i).first))
        {
            inactive->push_back(i);
        }
    }
    return inactive;
}

int ArrayGraph::getConnectedVertex()
{
    for (int i = 0; i < (int) graphState->size(); i++)
    {
        if (graphState->at(i).first && graphState->at(i).second >= 1)
        {
            return i;
        }
    }
    return -1;
}

// given an upper bound for the vertex degree (i.e. a pre-calculated max-degree), we can abort as soon as a vertex with a highest possible degree is found
int ArrayGraph::getMaxDegreeVertex()
{
    //int agg = 0;
    int max = -1;
    int maxIndex;

    for (int i = 0; i < (int) adjacencyList.size(); i++)
    {
        //if(2*numberOfEdges - agg + i <= max) { break; }
        if (graphState->at(i).first)
        {
            int degree = getVertexDegree(i);
            //agg += degree;
            if (max < degree)
            {
                max = degree;
                maxIndex = i;
            }
        }
    }
    return maxIndex;
}

int ArrayGraph::getFirstActiveVertex()
{
    for (int i = 0; i < (int) graphState->size(); i++)
    {
        if(graphState->at(i).first)
        {
            return i;
        }
    }
    return -1; // TODO: -1 is dummy return
}

// TODO: maybe provide an array that the neighbours are written into
// Then in cases, where there is no need to allocate a new array, that time can be saved
std::vector<int>* ArrayGraph::getNeighbours(int vertexIndex)
{
    std::vector<int>* neighbours = new std::vector<int>();
    for (int i = 0; i < (int) adjacencyList[vertexIndex]->size(); i++)
    {
        if (graphState->at(adjacencyList[vertexIndex]->at(i)).first)
        {
            neighbours->push_back(adjacencyList[vertexIndex]->at(i));
        }
    }
    return neighbours;
}

void ArrayGraph::setInactive(std::vector<int>* vertexIndices)
{
    for (int i = 0; i < (int) vertexIndices->size(); i++)
    {
        setInactive(vertexIndices->at(i));
        for(int j=0; j< (int) adjacencyList[vertexIndices->at(i)]->size(); j++)
        {
            graphState->at(adjacencyList[vertexIndices->at(i)]->at(j)).second -= 1;
            //if(graphState->at(adjacencyList[vertexIndices->at(i)]->at(j)).first) { numberOfEdges--; }
        }
        //numberOfVertices--;
    }
}

void ArrayGraph::setActive(std::vector<int>* vertexIndices)
{
    for (int i = 0; i < (int) vertexIndices->size(); i++)
    {
        setActive(vertexIndices->at(i));
        for(int j=0; j< (int) adjacencyList[vertexIndices->at(i)]->size(); j++)
        {
            graphState->at(adjacencyList[vertexIndices->at(i)]->at(j)).second += 1;
            //if(graphState->at(adjacencyList[vertexIndices->at(i)]->at(j)).first) { numberOfEdges++; }
        }
        //numberOfVertices++;
    }
}


int ArrayGraph::getCycleBound()
{
    cycleNumber = 0;

    // TODO: Find a better spot to allocate/clear/delete these vectors
    // for now at every bound check => inefficient!
    cycles = new std::vector<std::vector<int>>;

    color = new std::vector<int>(numberOfVertices, -1);
    par = new std::vector<int>(numberOfVertices, -1);

    // TODO: Check with which index to start, maybe there's a more efficient way
    dfs_cycle(0, -1);

    printCycles();

    int lowerBound = 0;
    for (int i = 0; i < cycleNumber; i++) {
        int cycleSize = (*cycles)[i].size();
        if(cycleSize > 2)
            lowerBound += (int) std::ceil(cycleSize/2.f);
    }

    return lowerBound;
}



// Function to mark the vertex with
// different colors for different cycles
void ArrayGraph::dfs_cycle(int u, int p)
{

    // already (completely) visited vertex.
    if ((*color)[u] == 2) {
        return;
    }

    // seen vertex, but was not completely visited -> cycle detected.
    // backtrack based on parents to find the complete cycle.
    if ((*color)[u] == 1) {
        std::vector<int> v;
        cycleNumber++;

        int cur = p;
        v.push_back(cur);

        // backtrack the vertex which are
        // in the current cycle thats found
        while (cur != u) {
            cur = (*par)[cur];
            v.push_back(cur);
        }
        (*cycles).push_back(v);
        return;
    }
    (*par)[u] = p;

    // partially visited.
    (*color)[u] = 1;

    // simple dfs on graph
    std::vector<int> neighbours = *getNeighbours(u);
    for (int v : neighbours) {

        // if it has not been visited previously
        if (v == par->at(u)) {
            continue;
        }
        dfs_cycle(v, u);
    }

    // completely visited.
    (*color)[u] = 2;
}

// Function to print the cycles
void ArrayGraph::printCycles()
{

    // print all the vertex with same cycle
    for (int i = 0; i < cycleNumber; i++) {
        // Print the i-th cycle
        std::cout << "Cycle Number " << i + 1 << ": ";
        for (int x : (*cycles)[i])
            std::cout << x << " ";
        std::cout << std::endl;
    }
}

std::vector<std::string>* ArrayGraph::getStringsFromVertexIndices(std::vector<int>* vertices)
{
    std::vector<std::string>* solution = new std::vector<std::string>();
    for (int i = 0; i < (int) vertices->size(); i++)
    {
        for (auto entry : originalVertexNames)
        {
            if(entry.second.first == vertices->at(i))
            {
                solution->push_back(entry.first);
            }
        }
    }
    return solution;
}

void ArrayGraph::printMappings(std::vector<int>* vertices)
{
    std::cout << "Mapping string --> index:" << std::endl;
    for (int i = 0; i < (int) vertices->size(); i++)
    {
        for (auto entry : originalVertexNames)
        {
            if(entry.second.first == vertices->at(i))
            {
                std::cout << entry.first << " --> " << vertices->at(i) << std::endl;
            }
        }
    }
    std::cout << std::endl;
}