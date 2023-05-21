#include <iostream> //output streams
#include <fstream>  //ifstream file opening
#include <cmath> //ceil 
#include <algorithm> //ceil
#include <stack>          // std::stack

#include "ArrayGraph.h"
#include "BipartiteArrayGraph.h"


/*----------------------------------------------------------*/
/*-----------------   Graph Construction   -----------------*/
/*----------------------------------------------------------*/

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

void ArrayGraph::printOriginalVertexNames()
{
    for (auto entry : originalVertexNames)
    {
        std::cout << entry.first << ": index = " << entry.second.first << ", degree = " << entry.second.second << std::endl;
    }
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

/*----------------------------------------------------------*/
/*-------------------   Graph Utility   --------------------*/
/*----------------------------------------------------------*/

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

void ArrayGraph::printMappings()
{
    std::cout << "Mapping string --> index:" << std::endl;
    for (auto entry : originalVertexNames)
    {
            std::cout << entry.first << " --> " << entry.second.first << std::endl;
    }
    std::cout << std::endl;
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
    int maxIndex = -1;

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

int ArrayGraph::getMaxDegreeVertex(std::vector<int>* candidates)
{
    //int agg = 0;
    int max = -1;
    int maxIndex = -1;

    for (int candidate : *candidates)
    {
        //if(2*numberOfEdges - agg + i <= max) { break; }
        if (graphState->at(candidate).first)
        {
            int degree = getVertexDegree(candidate);
            //agg += degree;
            if (max < degree)
            {
                max = degree;
                maxIndex = candidate;
            }
        }
    }
    return maxIndex;
}

bool ArrayGraph::isVertexCoverFound()
{
    //for each vertex
    for(int i = 0; i < (int) adjacencyList.size(); i++)
    {
        //if active
        if(graphState->at(i).first)
        {
            //check if all neighbours
            for(int j = 0; j < (int) adjacencyList.at(i)->size(); j++)
            {
                //are inactive
                if(graphState->at(j).first)
                {
                    return false;
                }
            }
        }

    }
    return true;
}
/*
 * Iterate through graph and find first still active edge 
*/
std::pair<int, int>* ArrayGraph::getFirstValidEdge()
{
    //iterate through active vertices
    for (int i = 0; i < (int) graphState->size(); i++)
    {
        if(graphState->at(i).first)
        {
            //find valid edge
            for(int j = 0; j < (int) adjacencyList[i]->size(); j++)
            {
                if(graphState->at(adjacencyList[i]->at(j)).first)
                {
                    return new std::pair<int, int> ({i, adjacencyList[i]->at(j)});
                }
            }
        }
    }
    return nullptr;
}

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

/* get unison of neighbours of origins that arent already contained in origins */
std::vector<int>* ArrayGraph::getNeighbours(std::vector<int>* origins)
{
    std::vector<int>* neighbours = new std::vector<int>();
    for (int origin : *origins)
    {
        for (int i = 0; i < (int) adjacencyList[origin]->size(); i++)
        {
            if (!contains(origins, adjacencyList[origin]->at(i)) && graphState->at(adjacencyList[origin]->at(i)).first)
            {
                neighbours->push_back(adjacencyList[origin]->at(i));
            }
        }
    }
    return neighbours;
}

// Find the first component with size > 1 that contains any the origin points
std::vector<int> ArrayGraph::getFirstComponent(std::vector<int>* origins)
{
    std::stack<int> S = std::stack<int>();
    std::vector<int>* visited = new std::vector<int>();
    std::vector<int> component = std::vector<int>();
    int current;
    for (int i = 0; i < (int) origins->size(); i++)
    {
        if(contains(visited, origins->at(i))) continue;
        while(!S.empty()) { S.pop(); }
        while(!visited->empty()) { visited->pop_back(); }
        S.push(origins->at(i));
        visited->push_back(origins->at(i));
        //std::cout << "DFS for " << origins->at(i) << "\n";
        while(!S.empty())
        {
            current = S.top();
            S.pop();
            component.push_back(current);
            for(int neighbour : *getNeighbours(current))
            {
                if(contains(visited, neighbour)) continue;
                visited->push_back(neighbour);
                S.push(neighbour);
            }
        }
        //std::cout << "Pushing\n";
        if(component.size() <= 1) continue;
        return component;
    }
    return component;
}

// Find the components with size > 1 containing the origin points
std::vector<std::vector<int>>* ArrayGraph::getComponents(std::vector<int>* origins)
{
    std::vector<std::vector<int>>* components = new std::vector<std::vector<int>>();
    std::stack<int> S = std::stack<int>();
    std::vector<int>* visited = new std::vector<int>();
    int current;
    for (int i = 0; i < (int) origins->size(); i++)
    {
        if(contains(visited, origins->at(i))) continue;
        std::vector<int> component = std::vector<int>();
        while(!S.empty()) { S.pop(); }
        S.push(origins->at(i));
        visited->push_back(origins->at(i));
        //std::cout << "DFS for " << origins->at(i) << "\n";
        while(!S.empty())
        {
            current = S.top();
            S.pop();
            component.push_back(current);
            for(int neighbour : *getNeighbours(current))
            {
                if(contains(visited, neighbour)) continue;
                visited->push_back(neighbour);
                S.push(neighbour);
            }
        }
        //std::cout << "Pushing\n";
        if(component.size() <= 1) continue;
        components->push_back(component);
    }
    //std::cout << "Returning\n";
    return components;
}


void ArrayGraph::setInactive(int vertexIndex)
{
    graphState->at(vertexIndex).first = false;
    for(int j=0; j< (int) adjacencyList[vertexIndex]->size(); j++)
    {
        graphState->at(adjacencyList[vertexIndex]->at(j)).second -= 1;
        //if(graphState->at(adjacencyList[vertexIndices->at(i)]->at(j)).first) { numberOfEdges--; }
    }
    //numberOfVertices--;
}

void ArrayGraph::setInactive(std::vector<int>* vertexIndices)
{
    for (int i = 0; i < (int) vertexIndices->size(); i++)
    {
        setInactive(vertexIndices->at(i));
    }
}

void ArrayGraph::setActive(int vertexIndex)
{
    graphState->at(vertexIndex).first = true; 
    for(int j=0; j< (int) adjacencyList[vertexIndex]->size(); j++)
    {
        graphState->at(adjacencyList[vertexIndex]->at(j)).second += 1;
        //if(graphState->at(adjacencyList[vertexIndices->at(i)]->at(j)).first) { numberOfEdges++; }
    }
    //numberOfVertices++;
}

void ArrayGraph::setActive(std::vector<int>* vertexIndices)
{
    for (int i = 0; i < (int) vertexIndices->size(); i++)
    {
        setActive(vertexIndices->at(i));
    }
}

/*----------------------------------------------------------*/
/*------------------   Calculate Bounds   ------------------*/
/*----------------------------------------------------------*/

int ArrayGraph::getLowerBoundVC() {

    int cliqueBound = getCliqueBound();
    //int cycleBound = getCycleBound();
    //return cycleBound;
    //return getLPBound();
    //return getLPCycleBound();
    return cliqueBound;
}

std::vector<int> ArrayGraph::getAllLowerBounds() {
    int cliqueBound = getCliqueBound();
    int cycleBound = 0; //getCycleBound(); // TODO: cycle bound throws access error
    int lpBound = getLPBound();
    int lpCycleBound = getLPCycleBound();
    return std::vector<int>({cliqueBound, cycleBound, lpBound, lpCycleBound});
}

std::vector<int>* ArrayGraph::getVerticesSortedByDegree()
{
    std::vector<int>* sorted = new std::vector<int>();
    //accumulate vertices
    for (int i=0; i<numberOfVertices; i++)
    {
        if(graphState->at(i).first) sorted->push_back(i);
    }
    quickSort(sorted, 0, numberOfVertices-1);
    /* for (int i: *sorted)
    	std::cout << sorted->at(i) << ": " << getVertexDegree(i) << "\n"; */
    return sorted;
}

bool ArrayGraph::vertexCanBeAddedToClique(int vertex, std::vector<int>* clique)
{
    //each vertex in clique
    for (int i = 0; i < (int) clique->size(); i++)
    {
        //is a neighbour of vertex
        bool isNeighbour = false;
        for(int j = 0; j < (int) adjacencyList[vertex]->size(); j++) // TODO: add condition whether vertex is active if needed
        {
            if (adjacencyList[vertex]->at(j) == clique->at(i))
            {
                isNeighbour = true;
                break;
            }
        }
        if(!isNeighbour)
        {
            return false;
        }
    }
    return true;
}

/*
* Calculates clique cover with greedy heuristic
*/
int ArrayGraph::getCliqueBound()
{
    std::vector<int>* sorted_vertices = getVerticesSortedByDegree();
    std::vector<std::vector<int>*> cliques = std::vector<std::vector<int>*>();

    //for each vertex
    for(int i = 0; i < (int) sorted_vertices->size(); i++)
    {
        int curVertex = sorted_vertices->at(i);
        //search for clique to add to
        std::pair<int, int> maxClique = std::pair<int, int>({-1, -1}); //clique index, clique size
        for(int j = 0; j < (int) cliques.size(); j++)
        {
            bool canBeAdded = vertexCanBeAddedToClique(curVertex, cliques.at(j));
            if(canBeAdded && (int) cliques.at(j)->size() > maxClique.second)
            {
                maxClique.first = j;
                maxClique.second = cliques.at(j)->size();
            }
        }
        //no clique found
        if(maxClique.first == -1)
        {
            //create clique
            std::vector<int>* newClique = new std::vector<int>();
            newClique->push_back(curVertex);
            cliques.push_back(newClique);
        }
        else
        {
            cliques[maxClique.first]->push_back(curVertex);
        }
    }

    //aggregate clique bound
    int cliqueBound = 0;
    for(int i = 0; i < (int) cliques.size(); i++)
    {
        cliqueBound += cliques.at(i)->size() - 1;
    }
    return cliqueBound;
}

/* Hopcroft Karp Algorithm */
int ArrayGraph::getLPBound()
{
    //generate bipartite graph that splits vertices into left and right
    BipartiteArrayGraph* bipartiteGraph = BipartiteArrayGraph::createBipartiteGraphByVertexSplit(this);

    //execute Hopcroft Karp to the maximum matching size
    int maximumMatchingSize = bipartiteGraph->getMaximumMatching();

    return maximumMatchingSize/2;
};

int ArrayGraph::getCycleBound()
{
    int graphSize = graphState->size();
    if (graphSize == 0)
        return 0;
    cycleNumber = 0;

    // TODO: Find a better spot to allocate/clear/delete these vectors
    // for now at every bound check => inefficient!
    cycles = new std::vector<std::vector<int>>;

    color = new std::vector<int>(numberOfVertices, -1);
    par = new std::vector<int>(numberOfVertices, -1);

    // TODO: Check with which index to start, maybe there's a more efficient way
    // => take the maxDegree!
    int maxDegreeVertex = getMaxDegreeVertex();
    dfs_cycle(maxDegreeVertex, -1);

    std::list<int> disjointCyclesMax;

    disjointCyclesMax = getDisjointCycles();

    int lowerBound = 0;
    for (int i = 0; i < cycleNumber; i++) {
        int cycleSize = (*cycles)[i].size();
        if(cycleSize > 2) {
            int cycleVC = (int) std::ceil(cycleSize / 2.f);
            lowerBound += cycleVC;
        }
    }

    int lowerBoundDisMax = 0;
    int lowerBoundDis = lowerBoundDisMax;
    for (auto i: disjointCyclesMax) {
        int cycleSize = (*cycles)[i].size();
        if(cycleSize > 2) {
            int cycleVC = (int) std::ceil(cycleSize / 2.f);
            lowerBoundDisMax += cycleVC;
        }
    }

    lowerBoundDis = lowerBoundDisMax;

//    std::cout << "#recursive steps: " << lowerBoundDis << std::endl;

    delete cycles;
    delete color;
    delete par;

    return lowerBoundDis;
}

int ArrayGraph::getLPCycleBound()
{
    // generate bipartite graph that splits vertices into left and right
    BipartiteArrayGraph* bipartiteGraph = BipartiteArrayGraph::createBipartiteGraphByVertexSplit(this);

    // execute Hopcroft Karp to the maximum matching size
    int LPCycleBound = bipartiteGraph->getMaximumMatchingCycleBound();
    return LPCycleBound;
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
        int cur = p;
//
//        if ((*color)[cur] == 2)
//            return;
        std::vector<int> v;
        cycleNumber++;


        v.push_back(cur);

        //mark as visited to get disjointed cycles!
        (*color)[cur] = 2;
        // backtrack the vertex which are
        // in the current cycle thats found
        while (cur != u) {
            cur = (*par)[cur];

//            if ((*color)[cur] == 2){
//                cycleNumber--;
//                return;}

            //mark as visited to get disjointed cycles!
            (*color)[cur] = 2;

            v.push_back(cur);
        }
//        (*color)[cur] = 2;
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

std::list<int> ArrayGraph::getDisjointCycles()
{
    using namespace  std;

    list<int> disjointCycles;
    list<int> disjointAllCycles;

    //initialize indexList
    for (int i = 0; i < cycleNumber; ++i) {
        disjointCycles.push_back(i);
        disjointAllCycles.push_back(i);
        int cycleSize = (*cycles)[i].size();
    }

    disjointCycles = getDisjointSet(disjointCycles);

//    printVector(&disjointCycles, "Disjoint Cycles");

    return disjointCycles;

}

std::list<int> ArrayGraph::getDisjointSet(std::list<int> validCycles)
{

    std::unordered_map<int, int> countVertices;
    std::unordered_map<int, std::list<int>> vertexCycle;

//    printVector(&validCycles, "Valid Cycles");

    int max = -1;
    int maxIdx = -1;

    for (int i: validCycles) {
        for (int x : (*cycles)[i]){
            countVertices[x]++;
            vertexCycle[x].push_back(i);

            if(countVertices[x] > max){
                max = countVertices[x];
                maxIdx = x;
            }
        }
    }

    if(max == 1 || max == -1)
        return validCycles;

    std::list<int> fightingCycles = vertexCycle[maxIdx];

//    printVector(&fightingCycles, "Fighting Cycles");

    //remove all but first from validCycles
    int fightingSize = fightingCycles.size();
    std::list<int>::iterator it = fightingCycles.begin();

    int maxSize = -1;

    for (int i = 0; i < fightingSize; ++i) {
        int idx = *it;
        int cycleSize = (*cycles)[idx].size();
        if( cycleSize> maxSize )
            maxSize = idx;
        advance(it, 1);
    }

//    std::cout << "Cycle with biggest size: " << maxSize << std::endl;

    it = fightingCycles.begin();
    for (int i = 0; i < fightingSize; ++i) {
        if(*it != maxSize)
            validCycles.remove(*it);
        advance(it, 1);
    }

    validCycles = getDisjointSet(validCycles);
    return validCycles;
}

// Function to print the cycles
void ArrayGraph::printCycles()
{

    // print all the vertex with same cycle
    for (int i = 0; i < cycleNumber; i++) {
        // Print the i-th cycle
        std::cout << "Cycle Number " << i << ": ";
        for (int x : (*cycles)[i])
            std::cout << x << " ";
        std::cout << std::endl;
    }
}

/*----------------------------------------------------------*/
/*----------------------   Utility   -----------------------*/
/*----------------------------------------------------------*/

int ArrayGraph::partition(std::vector<int>* toSort, int low, int high)
{
    // choose the last element as the pivot
    int pivot = toSort->at(high);

    // temporary pivot index
    int i = low - 1;

    for (int j=low; j<high; j++)
    {
        // if the current element is less than or equal to the pivot
        if (getVertexDegree(toSort->at(j)) <= getVertexDegree(pivot))
        {
            // move the temporary pivot index forward
            i++;
            // swap the current element with the element at the temporary pivot index
            // std::iter_swap(toSort->begin()+i, toSort->begin()+j);
            int tmp = toSort->at(i);
            (*toSort)[i] = toSort->at(j);
            (*toSort)[j] = tmp;
        }
    }

    // move the pivot element to the correct pivot position (between the smaller and larger elements)
    i++;
    // std::iter_swap(toSort->begin()+i, toSort->begin()+high);
    int tmp = toSort->at(i);
    (*toSort)[i] = toSort->at(high);
    (*toSort)[high] = tmp;
    return i;
}

void ArrayGraph::quickSort(std::vector<int>* toSort, int low, int high)
{
    if (low >= 0 && high >= 0 && low < high)
    {
        int p = partition(toSort, low, high);
        quickSort(toSort, low, p - 1);
        quickSort(toSort, p + 1, high);
    }
}

void ArrayGraph::printVector(std::list<int> *vec, std::string name)
{
    std::cout << name << " are: " << std::endl;
    for (auto itr = (*vec).begin();
         itr != (*vec).end(); itr++) {
        std::cout << *itr << " ";
    }
    std::cout << std::endl;
}

/**
 * List of vertices contains a certain vertex
*/
bool ArrayGraph::contains(std::vector<int>* vertexIndices, int vertexIndex)
{
    for(int i=0; i<(int) vertexIndices->size(); i++)
    {
        if(vertexIndices->at(i) == vertexIndex) return true;
    }
    return false;
}





