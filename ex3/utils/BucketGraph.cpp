#include <iostream> //output streams
#include <fstream>  //ifstream file opening
#include <math.h>
#include <chrono>
using namespace std::chrono;

#include "BucketGraph.h"
#include "ColorPrint.h"
#include "BipartiteArrayGraph.h"

typedef ColorPrint cp;

/*----------------------------------------------------------*/
/*-----------------   Graph Construction   -----------------*/
/*----------------------------------------------------------*/

BucketGraph*  BucketGraph::readStandardInput()
{
    //init
	BucketGraph* G = new BucketGraph();
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

    G->initActiveList(edges); //sets vertex references and generates activeList
    G->initAdjMap(); //sets references in adjacency lists of vertices
    G->initBucketQueue(); // initialise bucket queue

    return G;
}

std::string BucketGraph::eraseLeadingTrailingWhitespacesFromString(std::string str)
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

bool BucketGraph::isVertexCharacter(char c)
{
	if ((c >= '0' && c <= '9') || (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == '_') {
		return true;
	}
	return false;
}

void BucketGraph::initActiveList(std::vector<std::pair<std::string, std::string>> edges)
{
    vertexReferences = std::vector<Vertex*>(originalVertexNames.size());

    //clear first
    for (int k = 0; k < (int) vertexReferences.size(); k++)
    {
        vertexReferences[k] = nullptr;
    }

    //add edges
    for (int k = 0; k < (int) edges.size(); k++)
    {
        //get vertices
        auto firstVertexEntry = originalVertexNames.find(edges[k].first);
        auto secondVertexEntry = originalVertexNames.find(edges[k].second);

        if (firstVertexEntry == originalVertexNames.end() || secondVertexEntry == originalVertexNames.end())
        {
            throw std::invalid_argument("readStandardInput: inconsistency: map of original names returned nothing");
        }


        int indexFirst = (*firstVertexEntry).second.first;
        int maxDegFirst = (*firstVertexEntry).second.second;
        int indexSecond = (*secondVertexEntry).second.first;
        int maxDegSecond = (*secondVertexEntry).second.second;

        //first vertex: insert into adjacency list
        if (vertexReferences[indexFirst] != nullptr) //edges list initialised
        {
            if (vertexReferences[indexFirst]->index == indexFirst && vertexReferences[indexFirst]->adj != nullptr)
            {
                vertexReferences[indexFirst]->adj->push_back(indexSecond);
            }
            else
            {
                throw std::invalid_argument("readStandardInput: inconsistency when trying to insert adjacent vertex to v1");
            }
        }
        else //edges list not yet initialised
        {
            Vertex* firstVertex = new Vertex(new std::vector<int>(), edges[k].first, indexFirst, maxDegFirst);
            vertexReferences[indexFirst] = firstVertex;
            //insert
            activeList.push_back(*firstVertex);
            if (maxDegFirst > 0)
            {
                (*firstVertex->adj).push_back(indexSecond);
            }
            else
            {
                throw std::invalid_argument("readStandardInput: trying to insert edge but maxDegFirst is 0");
            }
        }

        //second vertex: insert into adjacency list
        if (vertexReferences[indexSecond] != nullptr) //edges list initialised
        {
            if (vertexReferences[indexSecond]->index == indexSecond && vertexReferences[indexSecond]->adj != nullptr)
            {
                vertexReferences[indexSecond]->adj->push_back(indexFirst);
            }
            else
            {
                throw std::invalid_argument("readStandardInput: inconsistency when trying to insert adjacent vertex to v1");
            }
        }
        else //edges list not yet initialised
        {
            Vertex* secondVertex = new Vertex(new std::vector<int>(), edges[k].second, indexSecond, maxDegSecond);
            vertexReferences[indexSecond] = secondVertex;
            //insert
            activeList.push_back(*secondVertex);
            if (maxDegSecond > 0)
            {
                (*secondVertex->adj).push_back(indexFirst);
            }
            else
            {
                throw std::invalid_argument("readStandardInput: trying to insert edge but maxDegFirst is 0");
            }
        }
    }
}

void BucketGraph::initAdjMap()
{
    for(int i = 0; i < (int) vertexReferences.size(); i++)
    {
        if(vertexReferences[i] != nullptr)
        {
            /* vertexReferences[i]->adj_refs = list<Vertex, member_hook<Vertex, list_member_hook<>, &Vertex::member_hook_>>();
            //copy adj list
            for(int j = 0; j < (int) vertexReferences[i]->adj->size(); j++)
            {
                Vertex* v = vertexReferences[vertexReferences[i]->adj->at(j)];
                std::cout << i << " insert: "<< v->index << std::endl;
                vertexReferences[i]->adj_refs.push_back(*v);
                //vertexReferences[i]->adj_refs.push_back(*vertexReferences[vertexReferences[i]->adj->at(j)]);
            }*/
            vertexReferences[i]->adj_map = new std::unordered_map<int, bool>();
            for(int j = 0; j < (int) vertexReferences[i]->adj->size(); j++)
            {
                vertexReferences[i]->adj_map->insert({vertexReferences[i]->adj->at(j), true});
            }
        }
        else
        {
            throw std::invalid_argument("initAdjRefs: at least one vertex refence is nullptr");
        }
    }
}

void BucketGraph::initBucketQueue()
{
    // init bucketQueue
    bucketQueue = list<Bucket>();
    bool inserted = false;
    int maxDeg = -1;
    // add activeList vertices
    for (auto vertex = activeList.begin(); vertex != activeList.end(); ++vertex)
    {
        inserted = false;
        if(vertex->degree > maxDeg) { maxDeg = vertex->degree; }

        for (auto bucket = bucketQueue.begin(); bucket != bucketQueue.end(); ++bucket)
        {
            if(vertex->degree == bucket->degree)
            {
                bucket->insert({vertex->bucketVertex});
                inserted = true;
                break;
            }
            else if(vertex->degree < bucket->degree)
            {
                bucketQueue.insert(bucket, *(new Bucket(vertex->degree, {vertex->bucketVertex})));
                inserted = true;
                break;
            }
        }
        if(!inserted)
        {
            bucketQueue.insert(bucketQueue.end(), *(new Bucket(vertex->degree, {vertex->bucketVertex})));
        }
    }

    // init bucketReferences
    int deg = 0;
    bucketReferences = std::vector<Bucket*>(maxDeg+1);
    for (auto bucket = bucketQueue.begin(); bucket != bucketQueue.end(); ++bucket)
    {
        while(deg < bucket->degree)
        {
            bucketReferences[deg] = new Bucket(deg, {});
            deg++;
        }
        bucketReferences[bucket->degree] = &*bucket;
        deg++;
    }
}

bool BucketGraph::isAdjMapConsistent()
{
    for(int i = 0; i < (int) vertexReferences.size(); i++)
    {
        //check map map entries
        for(int j = 0; j < (int) vertexReferences[i]->adj->size(); j++)
        {

            auto find = vertexReferences[i]->adj_map->find(vertexReferences[i]->adj->at(j));
            if(find == vertexReferences[i]->adj_map->end())
            {
                return false;
            }
        }
        //check adj vector entries
        for(auto it = vertexReferences[i]->adj_map->begin(); it != vertexReferences[i]->adj_map->end(); ++it)
        {
            bool found = false;
            for(int j = 0; j < (int) vertexReferences[i]->adj->size(); j++)
            {
                if(it->first == vertexReferences[i]->adj->at(j))
                {
                    found = true;
                    break;
                }
            }
            if(!found)
            {
                return false;
            }
        }
    }
    return true;
}

/*----------------------------------------------------------*/
/*-------------------   Graph Utility   --------------------*/
/*----------------------------------------------------------*/

void BucketGraph::print()
{
    if (vertexReferences.size() > 0)
	{
		std::cout << "\n";
        std::cout << "Graph of size: " << vertexReferences.size() << ", active: " << activeList.size() << std::endl;
        std::cout << "BucketQueue of size: " << bucketQueue.size() << ", maxDegree: " << getMaxDegree() << std::endl;
        std::cout << "name <name>, index <index>(<degree>): <neighbours>" << std::endl;
		for (int i = 0; i < (int) vertexReferences.size(); i++)
		{
			if (vertexReferences[i] != nullptr)
            {
                std::cout << "name " << cp::dye(vertexReferences[i]->strName, 'y') << ", index " <<
                cp::dye(std::to_string(i), 'g') << "(" << cp::dye(std::to_string(vertexReferences[i]->degree), 'r') << ")" << ": ";

                //print neighbours
                if(vertexReferences[i]->adj != nullptr)
                {
                    if(vertexReferences[i]->adj->empty())
                    {
                        std::cout << "<empty>";
                    }
                    else
                    {
                        for (int j = 0; j < (int) vertexReferences[i]->adj->size() - 1; j++)
                        {
                            std::cout << vertexReferences[i]->adj->at(j) << ", ";
                        }
                        std::cout << vertexReferences[i]->adj->at(vertexReferences[i]->adj->size() - 1);
                    }
                }
                else
                {
                   std::cout << "<nullptr>";
                }
            }
            else
            {
                std::cout << "<empty reference for index>";
            }
			std::cout << "\n";
		}
		std::cout << "\n";
	}
}

void BucketGraph::printActiveList()
{
    std::cout << "active vertices: ";
    for(auto it = activeList.begin(); it != activeList.end(); ++it)
    {
        std::cout << it->index << ", ";
    }
    if(activeList.size() == 0)
    {
        std::cout << "<empty>";
    }
    std::cout << std::endl;
}

std::string tileStr(std::string toTile, int n) {
	std::string tiling = "";
	for (int i=0; i<n; i++)
	{
		tiling += toTile;
	}
	return tiling;
}

int strRepSize(int n) {
	int spaces = 1;
	while (n > 9)
	{
        n /= 10;
		spaces++;
	}
	return spaces;
}

void BucketGraph::printBucketQueue()
{
    int fieldSize = 4;
    int rowMin = 3;
    int rowMax = 7;
    int row;
    std::cout << "BucketQueue of size " << cp::dye(std::to_string(bucketQueue.size()), 'y') << " with maxDegree: " << getMaxDegree() << std::endl;
    for(auto bucket = bucketQueue.begin(); bucket != bucketQueue.end(); ++bucket)
    {
        row = std::min(rowMax, std::max((int) bucket->vertices.size(), rowMin));
        std::cout << "/-" << cp::dye("Bucket " + std::to_string(bucket->degree), 'w') << tileStr("-", (2+fieldSize)*row-8-strRepSize(bucket->degree)) << "\\" << std::endl;
        for (auto vertex = bucket->vertices.begin(); vertex != bucket->vertices.end();)
        {
            int offs = 0;

            std::cout << "| " << cp::dye(std::to_string(vertex->index), 'w');
            offs = fieldSize-strRepSize(vertex->index);
            ++vertex;
            int i = 1;
            for (; i < row && vertex != bucket->vertices.end(); i++, ++vertex)
            {
                std::cout << ", " << tileStr(" ", offs) << cp::dye(std::to_string(vertex->index), 'w');
                offs = fieldSize-strRepSize(vertex->index);
            }
            std::cout << tileStr(" ", (2+fieldSize)*(row-i));
            std::cout << tileStr(" ", offs) << " |" << std::endl;
        }
        std::cout << "\\" << tileStr("-", (2+fieldSize)*row) << "/" << std::endl << std::endl;
    }
}

std::vector<std::string>* BucketGraph::getStringsFromVertexIndices(std::vector<int>* vertices)
{
    std::vector<std::string>* solution = new std::vector<std::string>();
    for (int i = 0; i < (int) vertices->size(); i++)
    {
        std::string stringcpy = vertexReferences[vertices->at(i)]->strName;
        solution->push_back(stringcpy);
    }
    return solution;
}

bool BucketGraph::vertexHasEdgeTo(int vertex, int secondVertex)
{
    if(vertex >= (int) vertexReferences.size() || secondVertex >= (int) vertexReferences.size())
    {
        throw std::invalid_argument("vertexHasEdgeTo: given vertices do not exist");
    }
    if(vertexReferences[vertex]->adj_map->find(secondVertex) != vertexReferences[vertex]->adj_map->end())
    {
        return true;
    }
    return false;
}

void BucketGraph::removeFromBucketQueue(int degree, std::vector<BucketVertex*> vertices)
{
    bucketReferences[degree]->remove(vertices);
    if(bucketReferences[degree]->vertices.size() == 0) {
        bucketQueue.erase(bucketQueue.iterator_to(*bucketReferences[degree]));
    }
}

void BucketGraph::addToBucketQueue(int degree, std::vector<BucketVertex*> vertices)
{
    // extend bucketReferences
    if(degree > (int) bucketReferences.size()-1)
    {
        for (int i=bucketReferences.size(); i<=degree; i++)
        {
            bucketReferences.push_back(new Bucket(i, {}));
        }
    }
    // insert bucket into queue
    if((int) bucketReferences[degree]->vertices.size() == 0)
    {
        // search larger non-empty bucket in bucketQueue
        auto it = bucketQueue.begin();
        for(; it != bucketQueue.end(); ++it)
        {
            if(it->degree > degree) {
                bucketQueue.insert(it, *bucketReferences[degree]);
                break;
            }
        }

        if(it == bucketQueue.end())
        {
            bucketQueue.insert(bucketQueue.end(), *bucketReferences[degree]);
        }
        // ---------------------------------------------
        /* auto it = *bucketReferences[degree];
        for(; it != bucketQueue.end(); ++it)
        {
            if(it->degree > degree) {
                bucketQueue.insert(it, *bucketReferences[degree]);
                break;
            }
        }

        if(it == bucketQueue.end())
        {
            bucketQueue.insert(bucketQueue.end(), *bucketReferences[degree]);
        } */
    }
    bucketReferences[degree]->insert(vertices);
}

void BucketGraph::setActive(int vertexIndex)
{
    Vertex* v = vertexReferences[vertexIndex];
    if(v == nullptr || v->adj == nullptr)
    {
        throw std::invalid_argument("setActive: vertex " + std::to_string(vertexIndex) + " is nullptr");
    }

    v->isActive = true;
    addToBucketQueue(v->degree, {v->bucketVertex});
    //update degree of all adjacent nodes
    for(int i = 0; i < (int) v->adj->size(); i++)
    {
        if(vertexReferences[i] == nullptr)
        {
            throw std::invalid_argument("setActive: tried updating degree of vertex pointing to nullptr");
        }
        vertexReferences[(*v->adj)[i]]->degree++;
        if(!vertexReferences[(*v->adj)[i]]->isActive) { continue; }
        removeFromBucketQueue(vertexReferences[(*v->adj)[i]]->degree-1, {vertexReferences[(*v->adj)[i]]->bucketVertex});
        addToBucketQueue(vertexReferences[(*v->adj)[i]]->degree, {vertexReferences[(*v->adj)[i]]->bucketVertex});
    }
}

void BucketGraph::setInactive(int vertexIndex)
{
    Vertex* v = vertexReferences[vertexIndex];
    if(v == nullptr || v->adj == nullptr)
    {
        throw std::invalid_argument("setInactive: vertex " + std::to_string(vertexIndex) + " is nullptr");
    }

    v->isActive = false;
    removeFromBucketQueue(v->degree, {v->bucketVertex});
    //update degree of all adjacent nodes
    for(int i = 0; i < (int) v->adj->size(); i++)
    {
        if(vertexReferences[i] == nullptr)
        {
            throw std::invalid_argument("setInactive: tried updating degree of vertex pointing to nullptr");
        }
        vertexReferences[(*v->adj)[i]]->degree--;
        if(!vertexReferences[(*v->adj)[i]]->isActive) { continue; }
        removeFromBucketQueue(vertexReferences[(*v->adj)[i]]->degree+1, {vertexReferences[(*v->adj)[i]]->bucketVertex});
        addToBucketQueue(vertexReferences[(*v->adj)[i]]->degree, {vertexReferences[(*v->adj)[i]]->bucketVertex});
    }
}

void BucketGraph::setActive(std::vector<int>* vertexIndices)
{
    for(int i = 0; i < (int) vertexIndices->size(); i++)
    {
        setActive(vertexIndices->at(i));
    }
}

void BucketGraph::setInactive(std::vector<int>* vertexIndices)
{
    for(int i = 0; i < (int) vertexIndices->size(); i++)
    {
        setInactive(vertexIndices->at(i));
    }
}

std::vector<int>* BucketGraph::getNeighbours(int vertexIndex)
{
    std::vector<int>* neighbours = new std::vector<int>();
    Vertex* v = vertexReferences[vertexIndex];
    if(v == nullptr)
    {
        throw std::invalid_argument("getNeighbours: vertex " + std::to_string(vertexIndex) + " is nullptr");
    }

    for(int i = 0; i < (int) v->adj->size(); i++)
    {
        if(vertexReferences[i] == nullptr)
        {
            throw std::invalid_argument("getNeighbours: tried updating degree of vertex pointing to nullptr");
        }
        if(vertexReferences[(*v->adj)[i]]->isActive)
        {
            neighbours->push_back((*v->adj)[i]);
        }
    }
    return neighbours;
}

int BucketGraph::getMaxDegree()
{
    if(bucketQueue.empty())
        return -1;
    return bucketQueue.back().degree;
}


int BucketGraph::getMaxDegreeVertex()
{
    if(bucketQueue.empty())
        return -1;
    return bucketQueue.back().vertices.front().index;
}

int BucketGraph::getVertexDegree(int vertexIndex)
{
    Vertex* v = vertexReferences[vertexIndex];
    if(v == nullptr)
    {
        throw std::invalid_argument("getVertexDegree: vertex " + std::to_string(vertexIndex) + " is nullptr");
    }
    return v->degree;
}

list<BucketVertex>* BucketGraph::getVerticesOfDegree(int degree)
{
    if(degree < (int) bucketReferences.size())
    {
        return &(bucketReferences[degree]->vertices);
    }
    else
    {
        throw std::invalid_argument("getVerticesOfDegree: degree goes out of bounds");
    }
}

/*----------------------------------------------------------*/
/*--------------   Virtual Matching & Flow   ---------------*/
/*----------------------------------------------------------*/

bool BucketGraph::matchingBFS(std::vector<int>* pairU, std::vector<int>* pairV, std::vector<int>* dist, int NIL)
    {
        std::queue<int> Q = std::queue<int>();
        for (int i=0; i<(int) pairU->size(); i++)
        {
            if ((*pairU)[i] == NIL)
            {
                (*dist)[i] = 0;
                Q.push(i);
            }
            else
            {
                (*dist)[i] = INT32_MAX;
            }
        }
        (*dist)[NIL] = INT32_MAX;
        while(!Q.empty())
        {
            int u = Q.front(); // TODO: check if correct side to peek
            Q.pop();
            if ((*dist)[u] < (*dist)[NIL])
            {
                for (auto v = vertexReferences[u]->adj->begin(); v != vertexReferences[u]->adj->end(); ++v)
                {
                    if ((*dist)[(*pairV)[*v]] == INT32_MAX)
                    {
                        (*dist)[(*pairV)[*v]] = (*dist)[u] + 1;
                        Q.push((*pairV)[*v]);
                    }
                }
            }
        }
        return (*dist)[NIL] != INT32_MAX;
    }

    bool BucketGraph::matchingDFS(std::vector<int>* pairU, std::vector<int>* pairV, std::vector<int>* dist, int u, int NIL)
    {
        if (u != NIL)
        {
            for (auto v = vertexReferences[u]->adj->begin(); v != vertexReferences[u]->adj->end(); ++v)
            {
                if ((*dist)[pairV->at(*v)] == (*dist)[u] + 1)
                {
                    if (matchingDFS(pairU, pairV, dist, pairV->at(*v), NIL))
                    {
                        (*pairV)[*v] = u;
                        (*pairU)[u] = *v;
                        return true;
                    }
                }
            }
            (*dist)[u] = INT32_MAX;
            return false;
        }
        return true;
    }

    int BucketGraph::hopcroftKarpMatchingSize()
    {
        int NIL = vertexReferences.size();
        int matching = 0;
        // initialize pairU/pairV
        std::vector<int> pairU = std::vector<int>(vertexReferences.size()); // TODO: needs to be allocated?
        std::vector<int> pairV = std::vector<int>(vertexReferences.size());
        std::vector<int> dist = std::vector<int>(vertexReferences.size()+1);
        for(int i=0; i<(int) pairU.size(); i++)
        {
            pairU[i] = NIL;
            pairV[i] = NIL;
        }
        while (matchingBFS(&pairU, &pairV, &dist, NIL))
        {
            for (int i=0; i<(int) pairU.size(); i++)
            {
                if (pairU[i] == NIL)
                {
                    if (matchingDFS(&pairU, &pairV, &dist, i, NIL))
                    {
                        matching++;
                    }
                }
            }
        }
        return matching;
    }

/*----------------------------------------------------------*/
/*------------------   Calculate Bounds   ------------------*/
/*----------------------------------------------------------*/

int BucketGraph::getLowerBoundVC() {
    //return getCliqueBound();
    return getLPBound();
    //return getLPCycleBound();
}

int BucketGraph::getLPBound()
{
    return hopcroftKarpMatchingSize()/2;
}

/*
* Calculates clique cover with greedy heuristic
* @param k: stop calculating cliques if bound already surpasses k
*/
int BucketGraph::getCliqueBound(int k)
{
    //TODO:delete time debug
    //auto start = high_resolution_clock::now();

    std::vector<std::vector<int>*> cliques = std::vector<std::vector<int>*>();
    int cliqueBound = 0;

    //int iterations = 0;
    //int cliqueIterations = 0;

    //iterate through vertices in ascending order of degrees (buckets)
    auto bucket = bucketQueue.begin();
    while(bucket != bucketQueue.end())
    {
        //for each vertex in bucket
        for(auto jt = bucket->vertices.begin(); jt != bucket->vertices.end(); ++jt)
        {
            //iterations++;
            int curVertex = jt->index;
            int cliqueIndex = -1;
            std::pair<int, int> maxClique = std::pair<int, int>({-1, -1}); //clique index, clique size
            for(int l = 0; l < (int) cliques.size(); l++)
            {
                //cliqueIterations++;
                if(vertexCanBeAddedToClique(curVertex, cliques.at(l)))
                {
                    cliqueIndex = l;
                    break;
                }
                /* if(vertexCanBeAddedToClique(curVertex, cliques.at(l)) && (int) cliques.at(l)->size() > maxClique.second)
                {
                    maxClique.first = l;
                    maxClique.second = cliques.at(l)->size();
                } */
            }
            //no clique found
            /* if(maxClique.first == -1)
            {
                //create clique
                std::vector<int>* newClique = new std::vector<int>();
                newClique->push_back(curVertex);
                cliques.push_back(newClique);
            }
            else
            {
                cliques[maxClique.first]->push_back(curVertex);
            } */
            //found clique
            if(cliqueIndex >= 0)
            {
                cliques[cliqueIndex]->push_back(curVertex);
                cliqueBound++; //aggregate clique bound
                //early stopping
                if(cliqueBound > k)
                {
                    return cliqueBound;
                }
            }
            else //else create new clique
            {
                std::vector<int>* newClique = new std::vector<int>();
                newClique->push_back(curVertex);
                cliques.push_back(newClique);
            }
        }
        ++bucket;
    }

    /* for(int i = 0; i < (int) cliques.size(); i++)
    {
        cliqueBound += cliques.at(i)->size() - 1;
    } */
    //auto stop = std::chrono::high_resolution_clock::now();
    //auto duration = duration_cast<milliseconds>(stop - start);
    //std::cout << "numIterations: " << iterations << ", numCliqueIterations: " << cliqueIterations << ", duration: " << duration.count() << std::endl;
    return cliqueBound;
}

/* int BucketGraph::getLPBound()
{
    //generate bipartite graph that splits vertices into left and right
    BipartiteArrayGraph* bipartiteGraph = BipartiteArrayGraph::createBipartiteGraphByVertexSplit(this);

    //execute Hopcroft Karp to the maximum matching size
    int maximumMatchingSize = bipartiteGraph->getMaximumMatching();

    return maximumMatchingSize/2;
} */



bool BucketGraph::vertexCanBeAddedToClique(int vertex, std::vector<int>* clique)
{
    //for each vertex in clique
    for (int i = 0; i < (int) clique->size(); i++)
    {
        //is a neighbour of vertex
        if(!vertexHasEdgeTo(vertex, clique->at(i)))
        {
            return false;
        }
    }
    return true;
}
<<<<<<< HEAD


/*----------------------------------------------------------*/
/*------------------   Reduction Rules   -------------------*/
/*----------------------------------------------------------*/


bool ArrayGraph::applyReductionRules(int* k, std::vector<ReductionVertices>* reductionArray)
{
    // if both rules are not applicable
    if(!rule_HighDegree(k, reductionArray) && !rule_DegreeZero(reductionArray))
        //if Buss rule == true => no vertex cover
        if(rule_Buss(k))
            return false;

    // Apply rule 0, 1 & 2 at every deactivation?
//    rule_DegreeOne(k, reductionArray);
//    rule_DegreeTwo(k, reductionArray);

    return true;
}

bool ArrayGraph::rule_HighDegree(int *k, std::vector<ReductionVertices>* reductionArray)
{
    int cnt = 0;

    int maxDegIdx = getMaxDegreeVertex();
    int maxDeg = getVertexDegree(maxDegIdx);

    if (maxDegIdx == -1 || maxDeg == -1)
        return false;

    while(maxDeg > (*k))
    {
        std::vector<int> tooHighToHandle = *getVerticesDegree(maxDeg);
        maxDeg--;
        if(tooHighToHandle.empty())
            continue;

        for (auto i = 0; i < (int)  tooHighToHandle.size(); ++i) {
            cnt++;
            int vertexID = tooHighToHandle.at(i);

            // save deleted vertex
            ReductionVertices delVer;
            delVer.rule = 3;
            delVer.deletedVertices.push_back(vertexID);
            delVer.kDecrement = 1;
            delVer.savedAdjacency = getNeighbours(vertexID);
            reductionArray->push_back(delVer);

            setInactive(vertexID);
            (*k)--;
        }

    }

    if(cnt == 0)
        return false;
    return true;
}


bool ArrayGraph::rule_DegreeZero(std::vector<ReductionVertices>* reductionArray)
{
    std::vector<int> zeroDegree = *getVerticesDegree(0);

    if(zeroDegree.empty())
        return false;

    for (auto i = 0; i < (int)  zeroDegree.size(); ++i) {
        int vertexID = zeroDegree.at(i);

        // save deleted vertex
        ReductionVertices delVer;
        delVer.rule = 0;
        delVer.kDecrement = 0;
        delVer.deletedVertices.push_back(vertexID);
        reductionArray->push_back(delVer);

        setInactive(vertexID);
    }

    return true;
}


bool ArrayGraph::rule_Buss(int* k)
{
    int k_square = std::pow((*k), 2);

    // If |V| > k^2 + k || |E|>k^2 => no vertex cover
    if((numberOfVertices > k_square + (*k)) || (numberOfEdges > k_square))
        return true;
    return false;
}


void ArrayGraph::rule_DegreeOne(int* k, std::vector<ReductionVertices>* reductionArray)
{
    std::vector<int> degreeOne = *getVerticesDegree(1);

    if(degreeOne.empty())
        return;

    for (auto i = 0; i < (int)  degreeOne.size(); ++i) {
        int vertexID = degreeOne.at(i);

        // THIS SHOULD HOLD THROUGH
        int neighbour = getNeighbours(vertexID)->at(0);

        // save deleted vertex
        ReductionVertices delVer;
        delVer.rule = 1;
        delVer.kDecrement = 1;
        delVer.deletedVertices.push_back(neighbour);
        delVer.deletedVertices.push_back(vertexID);
        delVer.savedAdjacency = getNeighbours(neighbour);
//        ->push_back(getNeighbours(neighbour));
        reductionArray->push_back(delVer);

        setInactive(neighbour);
        setInactive(vertexID);
        (*k)--;
    }

}


void ArrayGraph::rule_DegreeTwo(int* k, std::vector<ReductionVertices>* reductionArray)
{
    std::vector<int> degreeTwo = *getVerticesDegree(2);

    if(degreeTwo.empty())
        return;

    for (auto i = 0; i < (int) degreeTwo.size(); ++i) {
        int vertexID = degreeTwo.at(i);

        std::vector<int>* neighbours = getNeighbours(vertexID);
        int neighbourOne = neighbours->at(0);
        int neighbourTwo = neighbours->at(1);

        // Get shortest Neighbourhood to go through less to check for connection
        int shortestNeighbourhoodIdx;
        int otherNeighbourIdx;
        if(getVertexDegree(neighbourOne) > getVertexDegree(neighbourTwo) )
        {
            shortestNeighbourhoodIdx = neighbourTwo;
            otherNeighbourIdx = neighbourOne;
        }
        else
        {
            shortestNeighbourhoodIdx = neighbourOne;
            otherNeighbourIdx = neighbourTwo;
        }

        std::vector<int>* shortestNeighbourhood = getNeighbours(shortestNeighbourhoodIdx);
        std::vector<int>* otherNeighbour = getNeighbours(otherNeighbourIdx);

        // save deleted vertex
        ReductionVertices delVer;
        delVer.rule = 2;
        delVer.deletedVertices.push_back(otherNeighbourIdx);            // at (0)
        delVer.deletedVertices.push_back(shortestNeighbourhoodIdx);     // at (1)
        delVer.deletedVertices.push_back(vertexID);                     // at (2)

        delVer.savedAdjacency = getNeighbours(vertexID);
//        ->push_back(getNeighbours(vertexID));

        // CASE The neighbours know each other
        if(contains(shortestNeighbourhood,otherNeighbourIdx))
        {
            delVer.kDecrement = 2;
//            setInactive(otherNeighbourIdx);
//            setInactive(shortestNeighbourhoodIdx);
            setInactive(vertexID);
            (*k)-=2;
        }
        // CASE Neighbours don't know each other => setInactive or do that in addReducedVertices?
        else
        {
            delVer.kDecrement = 1;
            setVertexAdjacency(vertexID, putAdjacencyTogether(shortestNeighbourhood, otherNeighbour));
//            setInactive(otherNeighbourIdx);
//            setInactive(shortestNeighbourhoodIdx);
            (*k)--;
        }
        setInactive(otherNeighbourIdx);
        setInactive(shortestNeighbourhoodIdx);

        reductionArray->push_back(delVer);
    }
}

void ArrayGraph::rule_Domination(int* k, std::vector<ReductionVertices>* reductionArray)
{
    int maxDegIdx = getMaxDegreeVertex();
    int maxDeg = getVertexDegree(maxDegIdx);

    if (maxDegIdx == -1 || maxDeg == -1)
        return false;

    // for now degree 2
    while(maxDeg > 2) {
        std::vector<int> vertices = *getVerticesDegree(maxDeg);
        maxDeg--;
        if(vertices.empty())
            continue;

        for(auto v: vertices)
        {

        }

    }
}

void ArrayGraph::addReducedVertices(std::vector<int>* S, std::vector<ReductionVertices>* reductionArray)
{
    if(reductionArray->empty())
        return;

    for (auto i = 0; i < (int) reductionArray->size(); ++i) {
        ReductionVertices curReduction = reductionArray->at(i);
        int curAddition = curReduction.deletedVertices.at(0);
        switch (curReduction.rule) {
            // Degree Zero Rule
            case 0:
                break;
            // Degree One Rule
            case 1:
                S->push_back(curAddition);
                break;
            // Degree Two Rule
            case 2:
                if (curReduction.kDecrement == 2)
                {
                   S->push_back(curReduction.deletedVertices.at(1));
                }
                S->push_back(curAddition);
                break;
            // High Degree Rule
            case 3:
                S->push_back(curAddition);
                break;
            // Domination Rule
            case 4:

                break;
            default:
                std::cout<< "There shouldn't be this rule";
                break;
        }
    }

    // TODO: Eventually slow, especially if pointer is getting deleted later
    reductionArray->clear();
}


void ArrayGraph::addBackReducedVertices(int *k, std::vector<ReductionVertices>* reductionArray)
{
    if(reductionArray->empty())
        return;

    for (auto i = 0; i < (int) reductionArray->size(); ++i) {
        ReductionVertices curReduction = reductionArray->at(i);
        int curVertex = curReduction.deletedVertices.at(0);
        switch (curReduction.rule) {
            // Degree Zero
            case 0:
                setActive(curVertex);
                break;
                // Degree One
            case 1:
                setActive(curVertex);
                setActive(curReduction.deletedVertices.at(1));
//                (*k)++;
                break;
            case 2:
                if (curReduction.kDecrement == 2)
                {
                    setActive(curReduction.deletedVertices.at(2));
                }
                setActive(curReduction.deletedVertices.at(1));
                setActive(curVertex);
                break;
            // High Degree Rule
            case 3:
                setVertexAdjacencyBack(curVertex, curReduction.savedAdjacency);
                setActive(curVertex);
//                (*k)++;
                break;
            case 4:

                break;
            default:
                std::cout<< "There shouldn't be this rule";
                break;
        }
        (*k) += curReduction.kDecrement;
    }

    // TODO: Eventually slow, especially if pointer is getting deleted later
    reductionArray->clear();
}
=======
>>>>>>> origin/main
