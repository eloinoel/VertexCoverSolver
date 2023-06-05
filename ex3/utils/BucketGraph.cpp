#include <iostream> //output streams
#include <fstream>  //ifstream file opening
#include <sstream>
#include <math.h>
#include <chrono>
using namespace std::chrono;

#include "BucketGraph.h"
#include "ColorPrint.h"
#include "BipartiteArrayGraph.h"
#include "Reductions.h"

typedef ColorPrint cp;

/*----------------------------------------------------------*/
/*-----------------   Graph Construction   -----------------*/
/*----------------------------------------------------------*/

BucketGraph* BucketGraph::readStandardInput()
{
    //init

	BucketGraph* G = new BucketGraph();
    int vertexIndex = 0;
    G->originalVertexNames = std::unordered_map<std::string, std::pair<int, int>>();
    G->edges = std::vector<std::pair<std::string, std::string>>(); //edges after reading file once, eg. edges[0] = ["a", "b"]
    
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
        G->edges.push_back(edge_pair);
    }
    
    G->initActiveList(); //sets vertex references and generates activeList
    G->initAdjMap(); //sets references in adjacency lists of vertices
    G->initBucketQueue(); // initialise bucket queue
    G->initMatching(); // LP Bound matching fields
    G->reductions = new Reductions();

    G->initDominationHelper(); // Domination
//    G->reductions->initDominationVector(G);
//    G->reductions->printDominationSets();
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

void BucketGraph::initActiveList()
{
    numEdges = edges.size();
    numVertices = originalVertexNames.size();
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
                bucketQueue.insert(bucket, *(new Bucket(vertex->degree, *(new std::vector<BucketVertex*>({vertex->bucketVertex})))));
                inserted = true;
                break;
            }
        }
        if(!inserted)
        {
            bucketQueue.insert(bucketQueue.end(), *(new Bucket(vertex->degree, *(new std::vector<BucketVertex*>({vertex->bucketVertex})))));
        }
    }

    // init bucketReferences
    int deg = 0;
    bucketReferences = std::vector<Bucket*>(maxDeg+1);
    for (auto bucket = bucketQueue.begin(); bucket != bucketQueue.end(); ++bucket)
    {
        while(deg < bucket->degree)
        {
            bucketReferences[deg] = new Bucket(deg, *(new std::vector<BucketVertex*>({})));
            deg++;
        }
        bucketReferences[bucket->degree] = &*bucket;
        deg++;
    }
}

void BucketGraph::initMatching()
{
    pairU = (* new std::vector<int>(vertexReferences.size()));
    pairV = (* new std::vector<int>(vertexReferences.size()));
    dist = (* new std::vector<int>(vertexReferences.size()+1));
    unmatched = new std::vector<int>();
    next_unmatched = new std::vector<int>();
    NIL = vertexReferences.size();
    for(int i=0; i<(int) pairU.size(); i++)
    {
        pairU[i] = NIL;
        pairV[i] = NIL;
        dist[i] = INT32_MAX;
    }
    dist[NIL] = INT32_MAX;
    // run matching algorithm to calculate initial matching
    hopcroftKarpMatchingSize();

    // init flow
    nv = vertexReferences.size();
    s = nv*2;
    t = nv*2+1;
    pred = std::vector<int>(nv*2+2);
    flow = std::vector<std::vector<int>>(nv*2+2);
    for(int i=0; i<nv*2+2; i++)
    {
        flow[i] = std::vector<int>(nv*2+2);
        pred[i] = -1;
        for(int j=0; j<nv*2+2; j++)
        {
            flow[i][j] = 0;
        }
    }
    for (int i=0; i<(int) pairU.size(); i++)
    {
        if(pairU[i] != NIL)
        {
            flow[s][i] = 1;
            flow[i][pairU[i]] = 1;
            flow[pairU[i]][t] = 1;
        }
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

BucketGraph* BucketGraph::resetGraph()
{
    BucketGraph* G = new BucketGraph();
    G->originalVertexNames = originalVertexNames;
    G->edges = edges;

    reductions = new Reductions();
    G->initActiveList();
    G->initAdjMap();
    G->initBucketQueue();
    return G;
}

/*----------------------------------------------------------*/
/*-------------------   Graph Utility   --------------------*/
/*----------------------------------------------------------*/

void BucketGraph::print()
{
    if (vertexReferences.size() > 0)
	{
		std::cout << "\n";
        std::cout << "Graph of size: " << vertexReferences.size() << ", active: " << activeList.size() << ", edges: " << numEdges << std::endl;
        std::cout << "BucketQueue of size: " << bucketQueue.size() << ", maxDegree: " << getMaxDegree() << std::endl;
        std::cout << "name <name>, index <index>(<degree>): <neighbours>" << std::endl;
		for (int i = 0; i < (int) vertexReferences.size(); i++)
		{
			if (vertexReferences[i] != nullptr)
            {
                std::cout << "name " << cp::dye(vertexReferences[i]->strName, 'y') << ", index ";
                if(vertexReferences[i]->isActive)
                {
                    std::cout << cp::dye(std::to_string(i), 'g');
                }
                else 
                {
                    std::cout << cp::dye(std::to_string(i), 'd');
                }
                std::cout << "(" << cp::dye(std::to_string(vertexReferences[i]->degree), 'r') << ")" << ": ";

                //print neighbours
                if(vertexReferences[i]->adj != nullptr)
                {
                    if(vertexReferences[i]->adj->empty())
                    {
                        std::cout << "<empty>";
                    }
                    else
                    {
                        std::vector<int> active = std::vector<int>();
                        std::vector<int> inactive = std::vector<int>();
                        for(int j = 0; j < (int) vertexReferences[i]->adj->size(); j++)
                        {
                            if(vertexReferences[vertexReferences[i]->adj->at(j)]->isActive)
                            {
                                active.push_back(vertexReferences[i]->adj->at(j));
                            }
                            else 
                            {
                                inactive.push_back(vertexReferences[i]->adj->at(j));
                            }
                        }
                        for(int j = 0; j < (int) active.size(); j++)
                        {
                            std::cout << active.at(j) << ", ";
                        }
                        for(int j = 0; j < (int) inactive.size() - 1; j++)
                        {
                            std::cout << cp::dye(std::to_string(inactive.at(j)), 'd') << ", ";
                        }
                        if(inactive.size() > 0)
                            std::cout << cp::dye(std::to_string(inactive.at(inactive.size() - 1)), 'd');

                        //std::cout << vertexReferences[i]->adj->at(vertexReferences[i]->adj->size() - 1);
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

void BucketGraph::printMatching()
{
    std::cout << "Current Matchings in the Graph: " << currentLPBound << std::endl;
    std::cout << "----- Left Matchings: -----" << std::endl;
    for(int i=0; i<(int) pairU.size(); i++) {
        if(pairU[i] != NIL) {
            std::cout << i << " >> " << pairU[i] << std::endl;
        }
    }
    std::cout << "----- Right Matchings: -----" << std::endl;
    for(int i=0; i<(int) pairV.size(); i++) {
        if(pairV[i] != NIL) {
            std::cout << pairV[i] << " << " << i << std::endl;
        }
    }
}

std::vector<std::string>* BucketGraph::getEdgesToConsoleString()
{
    std::vector<std::string>* str = new std::vector<std::string>();
    std::unordered_map<std::pair<int, int>, bool, boost::hash<std::pair<int, int>>>* edges_map = new std::unordered_map<std::pair<int, int>, bool, boost::hash<std::pair<int, int>>>();
    //for each active vertex
    for(int j = 0; j < (int) vertexReferences.size(); j++)
    {
        Vertex* v = vertexReferences.at(j);
        if(v && v->isActive)
        {
            //iterate through edges
            for(int i = 0; i < (int) v->adj->size(); i++)
            {
                if(isActive(v->adj->at(i)))
                {
                    //if edge already in edges map, skip
                    auto map_entry = edges_map->find(std::make_pair(v->index, v->adj->at(i)));
                    if(map_entry != edges_map->end()) continue;
                    map_entry = edges_map->find(std::make_pair(v->adj->at(i), v->index));
                    if(map_entry != edges_map->end()) continue;

                    //else insert edge
                    (*edges_map)[std::make_pair(v->index, v->adj->at(i))] = true;
                    //std::cout << it->strName << " " << vertexReferences[it->adj->at(i)]->strName << std::endl;
                    str->push_back(v->strName + " " + vertexReferences[v->adj->at(i)]->strName + "\n");
                    //str += v->strName + " " + vertexReferences[v->adj->at(i)]->strName + "\n";
                }
            }
        }
    }
    return str;
}

std::vector<std::string>* BucketGraph::getOriginalEdgesToConsoleString()
{
    std::vector<std::string>* str = new std::vector<std::string>();
    for(int i = 0; i < (int) edges.size(); i++)
    {
        str->push_back(edges.at(i).first + " " + edges.at(i).second + "\n");
        //std::cout << edges.at(i).first << " " << edges.at(i).second << std::endl;
    }
    return str;
}

int BucketGraph::getOriginalEdgeCount()
{
    return edges.size();
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

std::vector<std::string>* BucketGraph::getStringsFromVertexIndices(std::unordered_map<int, bool>* vertices)
{
    std::vector<std::string>* solution = new std::vector<std::string>();
    for (auto it = vertices->begin(); it != vertices->end(); ++it)
    {
        std::string stringcpy = vertexReferences[it->first]->strName;
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

void BucketGraph::removeFromBucketQueue(int vertex)
{
    int degree = vertexReferences[vertex]->degree;
    bucketReferences[degree]->remove(vertexReferences[vertex]->bucketVertex);
    if(bucketReferences[degree]->vertices.size() == 0) {
        bucketQueue.erase(bucketQueue.iterator_to(*bucketReferences[degree]));
    }
}

void BucketGraph::addToBucketQueue(int vertex)
{
    int degree = vertexReferences[vertex]->degree;
    // extend bucketReferences if needed
    if(degree > (int) bucketReferences.size()-1)
    {
        for (int i=bucketReferences.size(); i<=degree; i++)
        {
            bucketReferences.push_back(new Bucket(i, *(new std::vector<BucketVertex*>({}))));
        }
    }
    // insert bucket into queue if needed
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
    }
    bucketReferences[degree]->insert(vertexReferences[vertex]->bucketVertex);
}

void BucketGraph::moveToBiggerBucket(int degree, int vertex)
{
    //remove from current bucket
    Bucket* previousBucket = bucketReferences[degree];
    previousBucket->remove(vertexReferences[vertex]->bucketVertex);

    auto previous_bucket_it = bucketQueue.iterator_to(*previousBucket);
    auto next_bucket = previous_bucket_it;
    next_bucket++;
    //list_node<hook_defaults::void_pointer>* next_bucket = nullptr; 
    if(previousBucket->vertices.size() == 0) {
        bucketQueue.erase(previous_bucket_it);
    }

    degree++;
    //extend bucketReferences if bigger bucket doesnt exist yet
    if(degree > (int) bucketReferences.size()-1)
    {
        for (int i=bucketReferences.size(); i<=degree; i++)
        {
            bucketReferences.push_back(new Bucket(i, *(new std::vector<BucketVertex*>({}))));
        }
    }
    // insert bucket into queue
    if((int) bucketReferences[degree]->vertices.size() == 0) //bucket to insert doesn't exist in queue
    {
        //insert bucket into queue
        //next_bucket cannot be the bucket to insert to
        if(next_bucket->vertices.empty())
        {
            throw std::invalid_argument("moveToBiggerBucket: inconsistency, next bucket to insert is empty but still in bucket queue");
        }
        if(next_bucket == bucketQueue.end())
        {
            bucketQueue.insert(bucketQueue.end(), *bucketReferences[degree]);
        }
        else
        {
            bucketQueue.insert(next_bucket, *bucketReferences[degree]);
        }
    }
    bucketReferences[degree]->insert(vertexReferences[vertex]->bucketVertex);
}

void BucketGraph::moveToSmallerBucket(int degree, int vertex)
{
    //std::cout << "moveToSmallerBucket: deg = " << degree << std::endl;
    //printBucketQueue();
    //remove from current bucket
    Bucket* curBucket = bucketReferences[degree];
    curBucket->remove(vertexReferences[vertex]->bucketVertex);
    //std::cout << "after removeFromBucket: deg = " << degree << std::endl;

    //get iterator to previous bucket
    auto cur_bucket_it = bucketQueue.iterator_to(*curBucket);
    auto previous_bucket_it = cur_bucket_it;
    bool curBucketIsBegin = false;
    if(cur_bucket_it != bucketQueue.begin())
    {
        previous_bucket_it--;
    }
    else
    {
        curBucketIsBegin = true;
    }

    //list_node<hook_defaults::void_pointer>* next_bucket = nullptr; 
    if(curBucket->vertices.size() == 0) {
        bucketQueue.erase(cur_bucket_it);
    }

    degree--;
    // insert bucket into queue
    if((int) bucketReferences[degree]->vertices.size() == 0) //bucket to insert doesn't exist in queue
    {
        //previous_bucket cannot be the bucket to insert to
        if(!curBucketIsBegin && previous_bucket_it->vertices.empty()) //if curBucketIsBegin == true, previous_bucket_it hasn't been decremented
        {
            throw std::invalid_argument("moveToSmallerBucket: inconsistency, previous bucket to insert is empty but still in bucket queue");
        }
        //insert bucket into queue
        if(curBucketIsBegin) //beginning of list reached
        {
            bucketQueue.insert(bucketQueue.begin(), *bucketReferences[degree]);
        }
        else
        {
            previous_bucket_it++;
            bucketQueue.insert(previous_bucket_it, *bucketReferences[degree]);
        }
    }
    bucketReferences[degree]->insert(vertexReferences[vertex]->bucketVertex);
}


void BucketGraph::setActive(int vertexIndex)
{
    Vertex* v = vertexReferences[vertexIndex];
    if(v == nullptr || v->adj == nullptr) throw std::invalid_argument("setActive: vertex " + std::to_string(vertexIndex) + " is nullptr");
    if(v->isActive) throw std::invalid_argument("setActive: vertex " + std::to_string(vertexIndex) + "is already active");

    v->isActive = true;
    numVertices++;
    // FIXME: add back if active is needed
    //activeList.push_back(*v);

    addToBucketQueue(v->index); //TODO: optimise this to not loop by storing previous references in recursion
    //update degree of all adjacent nodes
    for(int i = 0; i < (int) v->adj->size(); i++)
    {
        if(vertexReferences[i] == nullptr) throw std::invalid_argument("setActive: tried updating degree of vertex pointing to nullptr");

        vertexReferences[(*v->adj)[i]]->degree++;
        if(!vertexReferences[(*v->adj)[i]]->isActive) { continue; }
        numEdges++;
        moveToBiggerBucket(vertexReferences[(*v->adj)[i]]->degree - 1, (*v->adj)[i]);
    }
    unmatched->push_back(vertexIndex);
    //std::cout << "Restored vertex: " << vertexIndex << std::endl;
}

void BucketGraph::setInactive(int vertexIndex)
{
    Vertex* v = vertexReferences[vertexIndex];
    if(v == nullptr || v->adj == nullptr) throw std::invalid_argument("setInactive: vertex " + std::to_string(vertexIndex) + " is nullptr");
    if(!v->isActive) throw std::invalid_argument("setInactive: vertex " + std::to_string(vertexIndex) + "is already not active");

    v->isActive = false;
    numVertices--;
    // FIXME: add back if active is needed
    //auto iter = activeList.iterator_to(*v);
    //activeList.erase(iter);

    //std::cout << "before removing from queue: " << v->index << std::endl;
    removeFromBucketQueue(v->index);
    //update degree of all adjacent nodes
    for(int i = 0; i < (int) v->adj->size(); i++)
    {
        if(vertexReferences[i] == nullptr) throw std::invalid_argument("setInactive: tried updating degree of vertex pointing to nullptr");
        
        vertexReferences[(*v->adj)[i]]->degree--;
        if(!vertexReferences[(*v->adj)[i]]->isActive) { continue; }
        numEdges--;
        //std::cout << "before moving neighbour " << (*v->adj)[i] << " to smaller bucket" << std::endl;
        moveToSmallerBucket(vertexReferences[(*v->adj)[i]]->degree + 1, (*v->adj)[i]);
        //std::cout << "after moving neighbour " << (*v->adj)[i] << " to smaller bucket" << std::endl;
    }

    //std::cout << "setInactive: after updating neighbours" << std::endl;
    // update matching
    //std::cout << "Checking edge (" << vertexIndex << ", " << pairU[vertexIndex] << "), when deleting " << vertexIndex;
    if(pairU[vertexIndex] != NIL) {
        //std::cout << " --> Decrementing";
        //unmatched.push_back(pairU[vertexIndex]);
        currentLPBound--;
    }
    //std::cout << std::endl;
    //std::cout << "Checking edge (" << pairV[vertexIndex] << ", " << vertexIndex << "), when deleting " << vertexIndex;
    if(pairV[vertexIndex] != NIL) {
        //std::cout << " --> Decrementing";
        unmatched->push_back(pairV[vertexIndex]);
        currentLPBound--;
    }
    //std::cout << std::endl;
    dist[vertexIndex] = INT32_MAX;
    pairV[pairU[vertexIndex]] = NIL;
    pairU[pairV[vertexIndex]] = NIL;
    pairU[vertexIndex] = NIL;
    pairV[vertexIndex] = NIL;
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
    Vertex* v = vertexReferences[vertexIndex];
    std::vector<int>* neighbours = new std::vector<int>();
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

int BucketGraph::getMaxDegreeVertexMinimisingNeighbourEdges()
{
    if(bucketQueue.empty())
        return -1;

    Bucket* bucket = &bucketQueue.back();
    if(bucket->degree > 1)
    {
        int minNeighbourEdgesVertex = -1;
        int minNeighbourEdges = INT_MAX;
        for(auto it = bucket->vertices.begin(); it != bucket->vertices.end(); ++it)
        {
            int numNeighbourEdges = 0;
            for(int i = 0; i < (int) vertexReferences[it->index]->adj->size(); i++)
            {
                if(vertexReferences[vertexReferences[it->index]->adj->at(i)]->isActive)
                {
                    numNeighbourEdges += vertexReferences[vertexReferences[it->index]->adj->at(i)]->degree;
                }
            }
            if(numNeighbourEdges < minNeighbourEdges)
            {
                minNeighbourEdges = numNeighbourEdges;
                minNeighbourEdgesVertex = it->index;
            }
        }
        if(minNeighbourEdgesVertex != -1)
        {
            return minNeighbourEdgesVertex;
        }
    }
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
        return nullptr;
    }
}

int BucketGraph::getFirstVertexOfDegree(int degree)
{
    if(degree < (int) bucketReferences.size())
    {
        if(bucketReferences[degree]->vertices.size() > 0)
        {
            return bucketReferences[degree]->vertices.front().index;
        }
    }
    else
    {
        throw std::invalid_argument("getVerticesOfDegree: degree goes out of bounds");
    }
    return -1;
}

int BucketGraph::getNthActiveNeighbour(int vertex, int n)
{
    if(vertex >= (int) vertexReferences.size())
    {
        throw std::invalid_argument("getFirstActiveNeighbour: vertex");
    }
    Vertex* v = vertexReferences[vertex];
    int counter = 0;
    for(int i = 0; i < (int) v->adj->size(); i++)
    {
        if(vertexReferences[v->adj->at(i)]->isActive)
        {
            if(counter == n)
            {
                return v->adj->at(i);
            }
            counter++;
        }
    }
    return -1;
}

std::pair<int, int>* BucketGraph::getFirstTwoActiveNeighbours(int vertex)
{
    std::pair<int, int>* neighbours = new std::pair<int, int>({-1, -1});
    if(vertex >= (int) vertexReferences.size())
    {
        throw std::invalid_argument("getFirstActiveNeighbour: vertex");
    }
    Vertex* v = vertexReferences[vertex];
    int counter = 0;
    for(int i = 0; i < (int) v->adj->size(); i++)
    {
        if(vertexReferences[v->adj->at(i)]->isActive)
        {
            if(counter == 0)
            {
                neighbours->first = v->adj->at(i);
            }
            else if (counter == 1)
            {
                neighbours->second = v->adj->at(i);
                return neighbours;
            }
            counter++;
        }
    }
    return neighbours;
}

int BucketGraph::getNumVertices()
{
    return numVertices;
}

int BucketGraph::getNumEdges()
{
    return numEdges;
}

int BucketGraph::bruteForceCalculateNumEdges()
{
    int edgeCount = 0;
    for(int i = 0; i < (int) bucketReferences.size(); i++)
    {
        edgeCount += bucketReferences[i]->degree * bucketReferences[i]->vertices.size();
    }
    return edgeCount/2;
}

/*----------------------------------------------------------*/
/*-------------------   Data Reduction   -------------------*/
/*----------------------------------------------------------*/

bool BucketGraph::reduce(int* k)
{
    int c = 0;
    while(true)
    {
        //initialisise
        RULE_APPLICATION_RESULT degreeTwoResult = INAPPLICABLE;
        RULE_APPLICATION_RESULT highDegreeResult = INAPPLICABLE;
        RULE_APPLICATION_RESULT degreeZeroResult = INAPPLICABLE;
        RULE_APPLICATION_RESULT degreeOneResult = INAPPLICABLE;
        RULE_APPLICATION_RESULT dominationResult = INAPPLICABLE;
        RULE_APPLICATION_RESULT LPFlowResult = INAPPLICABLE;

        highDegreeResult = reductions->rule_HighDegree(this, k);
        if(highDegreeResult == INSUFFICIENT_BUDGET) return true; //cut
        /*degreeZeroResult = reductions->rule_DegreeZero(this);
        if(highDegreeResult == INAPPLICABLE && degreeZeroResult == INAPPLICABLE)
        {
            if(reductions->rule_Buss(this, k, getNumVertices(), getNumEdges()) == APPLICABLE)
                return true;
        }
        degreeOneResult = reductions->rule_DegreeOne(this, k);
        if(degreeOneResult == INSUFFICIENT_BUDGET) return true; */ //cut

        //dominationResult = reductions->rule_Domination(this, k);
        dominationResult = reductions->rule_Domination(this, k);
        if(dominationResult == INSUFFICIENT_BUDGET) return true; //cut

        //TODO: debug merge 
        /* degreeTwoResult = reductions->rule_DegreeTwo(this, k);
        if(degreeTwoResult == INSUFFICIENT_BUDGET)
        {
            std::cout << cp::dye("deg2 cut", 'y') << std::endl;
            return true; //cut
        } */

        /* if(c > 0)
        {
            LPFlowResult = reductions->rule_LPFlow(this, k);
            if(LPFlowResult == INSUFFICIENT_BUDGET) return true;
        } */

        if(highDegreeResult == INAPPLICABLE && degreeZeroResult == INAPPLICABLE && degreeOneResult == INAPPLICABLE
         && degreeTwoResult == INAPPLICABLE && dominationResult == INAPPLICABLE && LPFlowResult == INAPPLICABLE) //TODO: add conditions for other rules
        {
            return false;
        }
    }
    return false;
}

void BucketGraph::unreduce(int* k, int previousK, std::unordered_map<int, bool>* vc)
{
    if(reductions->appliedRules->empty())
        return;

    //pop rules
    while(*k < previousK || (!reductions->appliedRules->empty() && reductions->appliedRules->back()->kDecrement == 0))
    {
        Reduction* rule = reductions->appliedRules->back();
        //std::cout << "> unreducing: ";
        switch(rule->rule)
        {
            case DEGREE_ZERO:
                setActive(rule->deletedVertices);
                break;
            case DEGREE_ONE:
                *k = *k + rule->kDecrement;
                setActive(rule->deletedVCVertices);
                if(vc != nullptr)
                {
                    for(int i = 0; i < (int) rule->deletedVCVertices->size(); i++)
                    {
                        vc->insert({rule->deletedVCVertices->at(i), true});
                    }
                    //vc->insert(vc->end(), rule->deletedVCVertices->begin(), rule->deletedVCVertices->end());
                }
                break;
            case DEGREE_TWO:
                //std::cout << "deg2 unreduce" << std::endl;
                *k = *k + rule->kDecrement;
                if(rule->mergeVertexInfo == nullptr)
                {
                    for(int i = 0; i < (int) rule->deletedVCVertices->size(); i++)
                    {
                        setActive(rule->deletedVCVertices->at(i));
                    }
                    setActive(rule->deletedVertices->front());
                    //add to solution
                    if(vc != nullptr)
                    {
                        //didn't merge, take neighbours of deg2 vertex
                        for(int i = 0; i < (int) rule->deletedVCVertices->size(); i++)
                        {
                            vc->insert({rule->deletedVCVertices->at(i), true});
                        }
                    }
                }
                else
                {
                    unmerge(rule); //handles setting vertices back active
                    if(vc != nullptr)
                    {
                        //if vc contains vertex that was merged into, add neighbours of deg2 vertex to the solution
                        auto it = vc->find(std::get<0>(*rule->mergeVertexInfo));
                        if(it != vc->end())
                        {
                            vc->erase(it);
                            for(int i = 0; i < (int) rule->deletedVCVertices->size(); i++)
                            {
                                vc->insert({rule->deletedVCVertices->at(i), true});
                            }
                        }
                        else //else merge vertex wasn't chosen, so its neighbours must be in vc --> add deg2 vertex to solution
                        {
                            vc->insert({rule->deletedVertices->front(), true});
                        }
                    }
                }
                break;
            case HIGH_DEGREE:
                *k = *k + rule->kDecrement;
                setActive(rule->deletedVCVertices);
                if(vc != nullptr)
                {
                    for(int i = 0; i < (int) rule->deletedVCVertices->size(); i++)
                    {
                        vc->insert({rule->deletedVCVertices->at(i), true});
                    }
                    //vc->insert(vc->end(), rule->deletedVCVertices->begin(), rule->deletedVCVertices->end());
                }
                break;
            case DOMINATION:
                //std::cout << cp::dye("Domination unreduce", 'g') << std::endl;
//                if((int) rule->deletedVCVertices->size() == 0)
//                    break;
                //print();
                /* std::cout << "{";
                if(rule->deletedVCVertices != nullptr)
                {
                    for (int j=0; j<(int) rule->deletedVCVertices->size(); j++)
                    {
                        if(vertexReferences[rule->deletedVCVertices->at(j)]->isActive) { std::cout << "!"; }
                        std::cout << rule->deletedVCVertices->at(j) << ", ";
                    }
                
                    std::cout << "}" << std::endl;
                }
                else
                {
                    std::cout << "nullptr" << std::endl;
                } */

                *k = *k + rule->kDecrement;
                setActive(rule->deletedVCVertices);
                //std::cout << "after set active" << std::endl;
                if(vc != nullptr)
                {
                    for(int i = 0; i < (int) rule->deletedVCVertices->size(); i++)
                    {
                        vc->insert({rule->deletedVCVertices->at(i), true});
                    }
                }
                //std::cout << "after unreduce" << std::endl;
                //print();
                break;
            case LPFLOW:
                *k = *k + rule->kDecrement;
                //std::cout << "Restoring component:";
                setActive(rule->deletedVertices);
                setActive(rule->deletedVCVertices);
                /* std::cout << " {";
                for (int j=0; j<(int) rule->deletedVertices->size(); j++)
                {
                    if(vertexReferences[rule->deletedVertices->at(j)]->isActive) { std::cout << "!"; }
                    std::cout << rule->deletedVertices->at(j) << ", ";
                }
                std::cout << "} / ";
                std::cout << "{";
                for (int j=0; j<(int) rule->deletedVCVertices->size(); j++)
                {
                    if(vertexReferences[rule->deletedVCVertices->at(j)]->isActive) { std::cout << "!"; }
                    std::cout << rule->deletedVCVertices->at(j) << ", ";
                }
                std::cout << "}" << std::endl; */
                if(vc != nullptr)
                {
                    for(int i = 0; i < (int) rule->deletedVCVertices->size(); i++)
                    {
                        vc->insert({rule->deletedVCVertices->at(i), true});
                    }
                    //vc->insert(vc->end(), rule->deletedVCVertices->begin(), rule->deletedVCVertices->end());
                }
                break;
            default:
                throw std::invalid_argument("unreduce error: unknown rule");
                break;
        }
        reductions->appliedRules->pop_back();
        //TODO: delete debug at the end
        if(*k > previousK)
        {
            throw std::invalid_argument("unreduce error: " + std::to_string(*k) + " > " + std::to_string(previousK) + ", stop coding garbage");
        }
    }
}

std::tuple<int, std::vector<int>*, std::unordered_map<int, bool>*, std::vector<int>*>* BucketGraph::merge(int v0, int v1, int v2)
{

    //we merge into the vertex with the max degree, for less duplicate checks
    int maxDegVertex = v2;
    int mergev0 = v0;
    int mergev1 = v1;
    if(vertexReferences[mergev0]->degree > vertexReferences[maxDegVertex]->degree) 
    {
        int tmp = maxDegVertex;
        maxDegVertex = mergev0;
        mergev0 = tmp;
    }
    if(vertexReferences[mergev1]->degree > vertexReferences[maxDegVertex]->degree)
    {
        int tmp = maxDegVertex;
        maxDegVertex = mergev1;
        mergev1 = tmp;
    }
    std::cout << cp::dye("merging into vertex " + std::to_string(maxDegVertex) + " <-- " + std::to_string(mergev0) + ", " + std::to_string(mergev1), 'p') << std::endl;

    /* std::cout << "Neighbour " << mergev0 << " (" << getVertexDegree(mergev0) << ") " << " neighbours: ";
    std::vector<int>* neighboursOfFirst = getNeighbours(mergev0);
    for(auto ne = neighboursOfFirst->begin(); ne != neighboursOfFirst->end(); ++ne)
    {
        std::cout << *ne << ", ";
    }
    std::cout << std::endl;
    std::cout << "Neighbour " << mergev1 << " with deg=" << getVertexDegree(mergev1) << " neighbours: ";
    neighboursOfFirst = getNeighbours(mergev1);
    for(auto ne = neighboursOfFirst->begin(); ne != neighboursOfFirst->end(); ++ne)
    {
        std::cout << *ne << ", ";
    }
    std::cout << std::endl; */
    
    /* if(mergev1 == 18)
    {
        print();
    } */
    //std::cout << "merge: maxDegv: " << maxDegVertex << ", mergev0: " << mergev0 << ", mergev1: " << mergev1 << std::endl;

    //std::cout << "merge_debug1.5" << std::endl;

    //save adjacency list of vertex to merge into
    std::vector<int>* adj_copy = new std::vector<int>(*vertexReferences[maxDegVertex]->adj);
    std::unordered_map<int, bool>* adj_map_copy = new std::unordered_map<int, bool>(*vertexReferences[maxDegVertex]->adj_map);
    std::vector<int>* added_vertices = new std::vector<int>();

    //std::cout << "changing adj lists " << std::endl;
    removeFromBucketQueue(maxDegVertex);
    //add edges of merging vertices to adj of vertex to merge into
    for(int i = 0; i < (int) vertexReferences[mergev0]->adj->size(); i++)
    {
        int neighbour = (int) vertexReferences[mergev0]->adj->at(i);
        if(isActive(neighbour) && neighbour != maxDegVertex && neighbour != mergev1 && !vertexHasEdgeTo(maxDegVertex, neighbour))
        {
            //std::cout << "added " << neighbour << std::endl;
            vertexReferences[maxDegVertex]->adj->push_back(neighbour);
            vertexReferences[maxDegVertex]->adj_map->insert({neighbour, true});
            //moveToBiggerBucket(vertexReferences[maxDegVertex]->degree, maxDegVertex);
            vertexReferences[maxDegVertex]->degree++;

            //also change adj lists of neighbour

            vertexReferences[neighbour]->adj->push_back(maxDegVertex);
            vertexReferences[neighbour]->adj_map->insert({maxDegVertex, true});
            vertexReferences[neighbour]->degree++; //because of previously setting inactive, degree is -1 and needs to be incremented, if edge to merge vertex is inserted
            moveToBiggerBucket(vertexReferences[neighbour]->degree-1, neighbour);

            //save for unmerging later
            added_vertices->push_back(neighbour);
        }
    }
    for(int i = 0; i < (int) vertexReferences[mergev1]->adj->size(); i++)
    {
        int neighbour = (int) vertexReferences[mergev1]->adj->at(i);
        if(isActive(neighbour) && neighbour != maxDegVertex && neighbour != mergev0 && !vertexHasEdgeTo(maxDegVertex, neighbour))
        {
            //std::cout << "added " << neighbour << std::endl;
            vertexReferences[maxDegVertex]->adj->push_back(neighbour);
            vertexReferences[maxDegVertex]->adj_map->insert({neighbour, true});
            //moveToBiggerBucket(vertexReferences[maxDegVertex]->degree, maxDegVertex);
            vertexReferences[maxDegVertex]->degree++;

            //also change adj lists of neighbour
            vertexReferences[neighbour]->adj->push_back(maxDegVertex);
            vertexReferences[neighbour]->adj_map->insert({maxDegVertex, true});
            vertexReferences[neighbour]->degree++;
            moveToBiggerBucket(vertexReferences[neighbour]->degree-1, neighbour);

            //save for unmerging later
            added_vertices->push_back(neighbour);
        }
    }
    //std::cout << "merge_debug3" << std::endl;
    //add merge vertex back to appropriate bucket
    addToBucketQueue(maxDegVertex);

    //std::cout << "merge_debug4" << std::endl;
    /* if(maxDegVertex == 18)
    {
        print();
    } */

    //print();
    //printBucketQueue();

    //set vertices that will be merged inactive
    //std::cout << "set inactive: " << mergev0 << std::endl;
    setInactive(mergev0);
    /* if(mergev1 == 18)
    {
        print();
        
    } */

    //std::cout << "set inactive: " << mergev1 << std::endl;
    setInactive(mergev1);

    /* if(maxDegVertex == 18)
    {
        print();
    } */
    /* print();
    printBucketQueue(); */
    std::tuple<int, std::vector<int>*, std::unordered_map<int, bool>*, std::vector<int>*>* mergeInfo = new std::tuple<int, std::vector<int>*,
     std::unordered_map<int, bool>*, std::vector<int>*>(maxDegVertex, adj_copy, adj_map_copy, added_vertices);
    return mergeInfo;
}


void BucketGraph::unmerge(Reduction* mergeRule)
{
    //print();
    //printBucketQueue();
    /* print();
    printBucketQueue(); */
    /* if(mergeRule->mergeVertexInfo == nullptr)
    {
        std::cout << "nullptr" << std::endl;
    } */
    //restore adjacency list and other fields of the merge vertex
    int mergeVertex = std::get<0>(*mergeRule->mergeVertexInfo);
    std::cout << cp::dye("unmerging vertex " + std::to_string(mergeVertex), 'g') << std::endl;
    /* for(int i = 0; i < (int) mergeRule->deletedVertices->size(); i++)
    {
        if(mergeRule->deletedVertices->at(i) != mergeVertex)
        {
            setActive(mergeRule->deletedVertices->at(i));
        }
    } */
    int v0 = mergeRule->deletedVertices->front();
    int v1 = mergeRule->deletedVCVertices->at(0);
    int v2 = mergeRule->deletedVCVertices->at(1);
    
    if(mergeVertex != v0) { setActive(v0); }
    if(mergeVertex != v1) { setActive(v1); }
    if(mergeVertex != v2) { setActive(v2); }

    //std::cout << "mergeVertex: " << mergeVertex << ", N's: " << v0 << ", " << v1 << ", " << v2 << std::endl;



    //std::cout << "before remove from bucket queue" << std::endl;
    removeFromBucketQueue(mergeVertex);
    //std::cout << "after remove from bucket queue" << std::endl;
    //delete merge vertex from adjacency lists of neighbours that he didn't have before merge
    //std::cout << "unmerge_debug1" << std::endl;
    std::vector<int>* added_vertices = std::get<3>(*mergeRule->mergeVertexInfo);
    for(int i = 0; i < (int) added_vertices->size(); i++)
    {
        //std::cout << "decrementing added vertex " << added_vertices->at(i) << std::endl;
        vertexReferences[added_vertices->at(i)]->adj->pop_back(); //should work if all rules are popped in right order
        /* for(auto jt = vertexReferences[added_vertices->at(i)]->adj->begin(); jt != vertexReferences[added_vertices->at(i)]->adj->end(); ++jt) //TODO: more efficient way
        {
            if(*jt == mergeVertex)
            {
                vertexReferences[added_vertices->at(i)]->adj->erase(jt);
                break;
            }
        } */
        //std::cout << "after pop" << std::endl;
        vertexReferences[added_vertices->at(i)]->adj_map->erase(vertexReferences[added_vertices->at(i)]->adj_map->find(mergeVertex));
        //std::cout << "after adj map" << std::endl;
        vertexReferences[added_vertices->at(i)]->degree--;
        //std::cout << "moving vertex " << added_vertices->at(i) << " with deg=" << vertexReferences[added_vertices->at(i)]->degree+1 << " to smaller bucket" << std::endl;
        /* print();
        printBucketQueue(); */
        if(vertexReferences[added_vertices->at(i)]->isActive) { moveToSmallerBucket(vertexReferences[added_vertices->at(i)]->degree+1, added_vertices->at(i)); }
        if(!vertexReferences[added_vertices->at(i)]->isActive) { std::cout << cp::dye("inconsistency, should be active while unmerging: ", 'r') << added_vertices->at(i) << std::endl; }
        //std::cout << "after move to smaller bucket" << std::endl;
        //moving to smaller bucket not needed because still inactive
    }
    //std::cout << "after added vertices edge decrement" << std::endl;
    //std::cout << "mergevertex adj size:" << vertexReferences[mergeVertex]->degree << std::endl;
    vertexReferences[mergeVertex]->adj = std::get<1>(*mergeRule->mergeVertexInfo);
    vertexReferences[mergeVertex]->adj_map = std::get<2>(*mergeRule->mergeVertexInfo);
    vertexReferences[mergeVertex]->degree = vertexReferences[mergeVertex]->adj->size();
    //std::cout << "mergevertex saved adj size:" << vertexReferences[mergeVertex]->degree << std::endl;

    //std::cout << "after mergevertex reset" << std::endl;
    //std::cout << "unmerge_debug2" << std::endl;
    addToBucketQueue(mergeVertex);
    //std::cout << "after add to bucket queue" << std::endl;
    //std::cout << "mergeVertex: " << mergeVertex << ", N's: " << v0 << ", " << v1 << ", " << v2 << std::endl;

    //std::cout << "unmerge_debug3" << std::endl;
    //set merged vertices active again
    //std::cout << "unmerge_debug4" << std::endl;
}




/*----------------------------------------------------------*/
/*--------------   Virtual Matching & Flow   ---------------*/
/*----------------------------------------------------------*/

bool BucketGraph::matchingBFS()
{
    std::queue<int> Q = std::queue<int>();
    if (!didInitialMatchingCalculation)
    {
        for (Vertex active : activeList)
        {
            if (pairU[active.index] == NIL)
            {
                dist[active.index] = 0;
                Q.push(active.index);
            }
            else
            {
                dist[active.index] = INT32_MAX;
            }
        }
    }
    else
    {
        for (auto vertex = unmatched->begin(); vertex != unmatched->end(); ++vertex)
        {
            if(!vertexReferences[*vertex]->isActive) { continue; }
            if (pairU[*vertex] == NIL)
            {
                dist[*vertex] = 0;
                Q.push(*vertex);
            }
            else
            {
                dist[*vertex] = INT32_MAX;
            }
        }
    }
    dist[NIL] = INT32_MAX;
    while(!Q.empty())
    {
        int u = Q.front();
        Q.pop();
        if (dist[u] < dist[NIL])
        {
            for (auto v = vertexReferences[u]->adj->begin(); v != vertexReferences[u]->adj->end(); ++v)
            {
                if(!vertexReferences[*v]->isActive) { continue; }
                if (dist[pairV[*v]] == INT32_MAX)
                {
                    dist[pairV[*v]] = dist[u] + 1;
                    Q.push(pairV[*v]);
                }
            }
        }
    }
    return dist[NIL] != INT32_MAX;
}

bool BucketGraph::matchingDFS(int u)
{
    if (u != NIL)
    {
        for (auto v = vertexReferences[u]->adj->begin(); v != vertexReferences[u]->adj->end(); ++v)
        {
            if(!vertexReferences[*v]->isActive) { continue; }
            if (dist[pairV[*v]] == dist[u] + 1)
            {
                if (matchingDFS(pairV[*v]))
                {
                    /* flow[u][pairU[u]] = 0;
                    flow[*v][pairV[*v]] = 0; */ // TODO: update flow, when updating matching
                    pairV[*v] = u;
                    pairU[u] = *v;
                    /* flow[u][pairU[u]] = 0;
                    flow[pairV[*v]][*v] = 0; */
                    return true;
                }
            }
        }
        dist[u] = INT32_MAX;
        return false;
    }
    return true;
}

int BucketGraph::hopcroftKarpMatchingSize()
{
    int matching = 0;
    if(!didInitialMatchingCalculation) {
        while (matchingBFS())
        {
            for (Vertex active : activeList)
            {
                if (pairU[active.index] == NIL)
                {
                    if (matchingDFS(active.index))
                    {
                        matching++;
                    }
                }
            }
        }
        currentLPBound = matching/* +1 */; //TODO: for some reason matching is always 1 smaller than the actual size of the matching ( Adding 1 for good measure :'P )
        didInitialMatchingCalculation = true;
    }
    else
    {
        /* std::cout << "Matching vertices {";
        for(auto vertex = unmatched->begin(); vertex != unmatched->end(); ++vertex)
        {
            std::cout << *vertex << ", ";
        }
        std::cout << "}" << std::endl; */
        while (matchingBFS())
        {
            for (auto vertex = unmatched->begin(); vertex != unmatched->end(); ++vertex)
            {
                if(!vertexReferences[*vertex]->isActive) { continue; }
                if (pairU[*vertex] == NIL)
                {
                    if (matchingDFS(*vertex))
                    {
                        matching++;
                        //std::cout << "Incrementing for new edge (" << *vertex << ", " << pairU[*vertex] << ")" << std::endl;
                    }
                }
            }
        }
        currentLPBound += matching;
        next_unmatched->clear();
        for(auto vertex = unmatched->begin(); vertex != unmatched->end(); ++vertex)
        {
            if (pairU[*vertex] == NIL)
            {
                next_unmatched->push_back(*vertex);
            }
        }
        unmatched->clear();
        *unmatched = *next_unmatched;
    }
    //std::cout << "post hopcroft karp" << std::endl;
    return currentLPBound;
}

int BucketGraph::edmondsKarpFlow()
{
    //std::cout << ".";
    //std::cout << "Beginning flow" << std::endl;
    // number of vertices helper variable (only for clarity)
    //onst int nv = vertexReferences.size();

    int flow_amt = 0;
    /* int s = nv*2;
    int t = nv*2+1; */
    //std::vector<int> pred = std::vector<int>(nv*2+2);
    std::queue<int> Q = std::queue<int>();
    int current;

    // init flow
    /* flow = std::vector<std::vector<int>>(nv*2+2);
    for(int i=0; i<nv*2+2; i++)
    {
        flow[i] = std::vector<int>(nv*2+2);
        pred[i] = -1;
        for(int j=0; j<nv*2+2; j++)
        {
            flow[i][j] = 0;
        }
    } */
    // set left to right flow according to matching
    /* for (int i=0; i<(int) pairU.size(); i++)
    {
        if(pairU[i] != NIL)
        {
            flow[s][i] = 1;
            flow[i][pairU[i]] = 1;
            flow[pairU[i]][t] = 1;
        }
    } */
    //std::cout << "Initialized flow" << std::endl;
    pred[t] = 1;
    // Until no augmenting path was found last iteration
    while (pred[t] != -1)
    {
        // clear queue & predecessors
        while (!Q.empty()) { Q.pop(); }
        for(int i=0; i<nv*2+2; i++) { pred[i] = -1; }
        // init queue with source vertex
        Q.push(s);
        //std::cout << "Pushed source vertex" << std::endl;
        // BFS to find the shortest s-t path
        while (!Q.empty())
        {
            // retrieve next queued vertex
            current = Q.front();
            Q.pop();

            // if current is s (capacity == 1)
            if(current == s)
            {
                //std::cout << "Evaluating source vertex" << std::endl;
                // s exclusively has edges to left vertices (indices 0 to vertexReferences.size()-1)
                for (int i=0; i<nv; i++)
                {
                    if(!vertexReferences[i]->isActive) { continue; }
                    // if left vertex i wasn't expanded & edge to it isn't already in flow
                    if (pred[i] == -1 && 1 > flow[current][i])
                    {
                        pred[i] = current;
                        Q.push(i);
                    }
                }
                //std::cout << "Evaluated source vertex" << std::endl;
            }
            // if current is a left vertex (capacity == INF for edges to right vertices)
            else if (current < nv)
            {
                //std::cout << "Evaluating left side vertex" << std::endl;
                // current has edges to right vertices (indices vertexReferences.size() to vertexReferences.size()*2-1)
                for (auto v=vertexReferences[current]->adj->begin(); v != vertexReferences[current]->adj->end(); ++v)
                {
                    if(!vertexReferences[*v]->isActive) { continue; }
                    if (pred[nv+*v] == -1)
                    {
                        pred[nv+*v] = current;
                        Q.push(nv+*v);
                    }
                }
                //std::cout << "Evaluated left side vertex" << std::endl;
            }
            // if current is a right vertex (capacity == 1 for edge to t and capacity == INF for the reverse edges to left vertices)
            else if (nv <= current && current < nv*2)
            {
                //std::cout << "Evaluating right side vertex: " << current << std::endl;
                // current has reverse edges to left vertices (indices vertexReferences.size() to vertexReferences.size()*2-1)
                for (auto v=vertexReferences[current-nv]->adj->begin(); v != vertexReferences[current-nv]->adj->end(); ++v)
                {
                    if(!vertexReferences[*v]->isActive) { continue; }
                    //std::cout << "Evaluating right side vertices left side neighbour" << std::endl;
                    if (pred[*v] == -1 && pairU[*v] == current)
                    {
                        pred[*v] = current;
                        Q.push(*v);
                    }
                }
                //std::cout << "Evaluated right side vertices reverse edges" << std::endl;
                if (pred[t] == -1 && 1 > flow[current][t])
                {
                    pred[t] = current;
                    break;
                    //std::cout << "Found s-t path" << std::endl;
                }
            }
        }
        //std::cout << "Finished DFS" << std::endl;

        if (pred[t] != -1)
        {
            /* std::cout << "Found flow-improving path" << std::endl;
            printMatching(); */
            // always: df = 1;
            for (int p = t; p != -1 && p != s; p = pred[p])
            {
                flow[pred[p]][p]++;
                flow[p][pred[p]]--;
                if(pred[p] == s || p == t) { continue; }
                // using reverse edge
                if (p < nv)
                {
                    pairV[pairU[p]] = NIL;
                    pairU[pairV[pred[p]-nv]] = NIL;

                    pairU[p] = NIL;
                    pairV[pred[p]-nv] = NIL;
                    //std::cout << "-" << "(" << p << ", " << pred[p]-nv << ") ";
                }
                // using forward edge
                else
                {
                    pairV[pairU[pred[p]]] = NIL;
                    pairU[pairV[p-nv]] = NIL;

                    pairU[pred[p]] = p-nv;
                    pairV[p-nv] = pred[p];
                    //std::cout << "+" << "(" << pred[p] << ", " << p-nv << ") ";
                }
            }
            /* std::cout << std::endl;
            printMatching(); */
            flow_amt++;
            //std::cout << "Processed flow-improving path" << std::endl;
        }
    }
    //std::cout << "Ending flow" << std::endl;
    return flow_amt;
}

void BucketGraph::getBipartMatchingFlowComponents(std::vector<int>* L, std::vector<int>* R)
{
    std::stack<int> S = std::stack<int>();
    std::vector<int> indices = std::vector<int>(pairU.size());
    std::vector<int> lowlink = std::vector<int>(pairU.size());
    std::vector<bool> onStack = std::vector<bool>(pairU.size());
    for(int i=0; i<(int) indices.size(); i++)
    {
        indices[i] = -1;
        lowlink[i] = -1;
        onStack[i] = false;
    }
    int index = 0;
    for(int i=0; i<(int) pairU.size(); i++)
    {
        if(indices[i] == -1)
        {
            strongconnect(&S, i, index, &indices, &lowlink, &onStack, L, R);
        }
    }
}
void BucketGraph::strongconnect(std::stack<int>* S, int vertex, int index, std::vector<int>* indices, std::vector<int>* lowlink, std::vector<bool>* onStack, std::vector<int>* L, std::vector<int>* R) {
    // Set the depth index for v to the smallest unused index
    S->push(vertex);
    int nv = vertexReferences.size();

    // vertex is left side vertex
    if (vertex < nv)
    {
        (*indices)[vertex] = index;
        (*lowlink)[vertex] = index;
        (*onStack)[vertex] = true;
        index++;
        // vertex has edges to right vertices (indices vertexReferences.size() to vertexReferences.size()*2-1)
        for (auto v=vertexReferences[vertex]->adj->begin(); v != vertexReferences[vertex]->adj->end(); ++v)
        {
            if(!vertexReferences[*v]->isActive) { continue; }
            if (indices->at(*v) == -1)
            {
                // Successor w has not yet been visited; recurse on it
                strongconnect(S, *v+nv, index, indices, lowlink, onStack, L, R);
                (*lowlink)[vertex] = std::min(lowlink->at(vertex), lowlink->at(*v));
            }
            else if (onStack->at(*v))
            {
                // Successor w is in stack S and hence in the current SCC
                // If w is not on stack, then (v, w) is an edge pointing to an SCC already found and must be ignored
                // Note: The next line may look odd - but is correct.
                // It says w.index not w.lowlink; that is deliberate and from the original paper
                (*lowlink)[vertex] = std::min(lowlink->at(vertex), indices->at(*v));
            }
        }
    }
    // vertex is right side vertex
    else if (nv <= vertex && vertex < nv*2)
    {
        (*indices)[vertex-nv] = index;
        (*lowlink)[vertex-nv] = index;
        (*onStack)[vertex-nv] = true;
        index++;
        int n = pairV[vertex-nv];
        if (n != NIL)
        {
            if (indices->at(n) == -1)
            {
                // Successor w has not yet been visited; recurse on it
                strongconnect(S, n, index, indices, lowlink, onStack, L, R);
                (*lowlink)[vertex-nv] = std::min(lowlink->at(vertex-nv), lowlink->at(n));
            }
            else if (onStack->at(n))
            {
                // Successor w is in stack S and hence in the current SCC
                // If w is not on stack, then (v, w) is an edge pointing to an SCC already found and must be ignored
                // Note: The next line may look odd - but is correct.
                // It says w.index not w.lowlink; that is deliberate and from the original paper
                (*lowlink)[vertex-nv] = std::min(lowlink->at(vertex-nv), indices->at(n));
            }
            // If v is a root node, pop the stack and generate an SCC
            if (lowlink->at(vertex-nv) == indices->at(vertex-nv))
            {
                int w = S->top();
                S->pop();
                if(w < nv) { L->push_back(w); }
                else { R->push_back(w-nv); }
                while(w != vertex)
                {
                    w = S->top();
                    S->pop();
                    if(w < nv) { (*onStack)[w] = false; L->push_back(w); }
                    else { (*onStack)[w-nv] = false; R->push_back(w-nv); }
                }
            }
        }
        // If v is a root node, pop the stack and generate an SCC
        if (lowlink->at(vertex-nv) == indices->at(vertex-nv))
        {
            int w = S->top();
            S->pop();
            if(w < nv) { L->push_back(w); }
            else { R->push_back(w-nv); }
            while(w != vertex)
            {
                w = S->top();
                S->pop();
                if(w < nv) { (*onStack)[w] = false; L->push_back(w); }
                else { (*onStack)[w] = false; R->push_back(w-nv); }
            }
        }
    }
}

void BucketGraph::setBipartMatchingFlowComponentsInactive(std::vector<int>* L, std::vector<int>* R, int k, double maxExecTime)
{
    auto start = std::chrono::high_resolution_clock::now();
    //int z = 88;
    //int x = 19;
    std::stack<int> S = std::stack<int>();
    std::vector<int> visited = std::vector<int>(pairU.size()); // 0-unvisited, 1-pending for evaluation, 2-visited
    bool isNotComp = false;
    //std::cout << "Searching for sources: 0-" << pairU.size()-1 << std::endl;
    // for each possible connected component (We exclude degree 0 connected components here)
    //print();
    for (int i=0; i<(int) pairU.size(); i++)
    {
        //if(visited[i] != 0) { continue; }
        //if(i==z) print();
        //if(i==z) std::cout << "Setting up for search from source: " << i << std::endl;
        if (!vertexReferences[i]->isActive || pairU[i] == NIL/* this culls deg=0 rule */) { continue; }
        for(int j=0; j<(int) visited.size(); j++) { visited[j] = 0; }
        while(!S.empty()) { S.pop(); }
        isNotComp = false;
        S.push(i);
        int target = pairU[i];
        std::vector<int> componentL = std::vector<int>();
        std::vector<int> componentR = std::vector<int>();
        //if(i==z) std::cout << cp::dye("Started search from source: ", 'g') << i << std::endl;
        // until component is closed
        while (!S.empty())
        {
            int current = S.top();
            //if(i==z) std::cout << cp::dye("Popped: ", 'p') << current << " from stack" << std::endl;
            if (!vertexReferences[current]->isActive) { isNotComp = true; break; } // TODO:
            // add vertex to component, if visited
            if(visited[current] == 2 && visited[pairU[current]] == 2)
            {
                //if(i==z) std::cout << cp::dye("Traversing back up from: ", 'y') << current << std::endl;
                S.pop();
                continue;
            }

            // expand node
            isNotComp = true;
            visited[current] = 2;
            visited[pairU[current]] = 2;
            //if(i==z) std::cout << cp::dye("Expanding: ", 'y') << current << " with pairU[" << current << "] = " << pairU[current] << std::endl;
            for (auto v=vertexReferences[current]->adj->begin(); v != vertexReferences[current]->adj->end(); ++v)
            {
                // We skip inactive vertices and the matched right vertex and vertices of which we already visited their left matched vertex
                if(!vertexReferences[*v]->isActive) { continue; }
                //if(i==z) std::cout << "Considering neighbour: " << *v << " and matching (" << pairV[*v] << ", " << *v <<")" << std::endl;
                // if we find an unmatched right vertex, abort immediately (This cannot be a component)
                if(pairV[*v] == NIL || !vertexReferences[pairV[*v]]->isActive) {
                    //std::cout << "Current: " << current << " has no left matched vertex." << std::endl;
                    isNotComp = true; break;
                }
                // if found target, this path is valid (we then need to check the other paths to be valid)
                if(pairU[i] == *v) { isNotComp = false; continue; }
                //if(current == pairV[*v]) { componentR.push_back(*v); }
                //if(visited[pairV[*v]] && !visited[*v]) { componentR.push_back(*v); visited[*v] = 2; }    // TODO: debug
                // We skip the matched right vertex and vertices of which we already visited their left matched vertex
                if(current == pairV[*v] || visited[*v] == 2) { continue; }
                if(visited[pairV[*v]] == 2 || visited[*v] == 1 || visited[pairV[*v]] == 1) { isNotComp = true; break; }
                // push next left vertex
                S.push(pairV[*v]);
                isNotComp = false;
                visited[*v] = 1;
                visited[pairV[*v]] = 1;
                //if(i==z) std::cout << cp::dye("Pushing neighbour: ", 'p') << pairV[*v] << " from matching (" << pairV[*v] << ", " << *v <<")" << " to stack" << std::endl;
            }
            if(isNotComp) {
                //if(i==z) std::cout << cp::dye("Current: ", 'r') << current << cp::dye(" has no successor vertex.", 'r') << std::endl;
                break;
            }

            //if(i==z) std::cout << cp::dye("Adding ", 'y') << current << " and: " << pairU[current] << " into the solution" << std::endl;
            componentL.push_back(current);
            componentR.push_back(pairU[current]);
        }
        //std::cout << "Checking if component was found" << std::endl;
        if(isNotComp) { continue; }
        // found component
        // iterate through componentL and calc component right vertices
        /* if((int) componentR.size() > 0) {
            std::cout << "Found component: {";
            for (int j=0; j<(int) componentL.size(); j++)
            {
                std::cout << componentL[j] << ", ";
            }
            std::cout << "} / ";
            std::cout << "{";
            for (int j=0; j<(int) componentR.size(); j++)
            {
                std::cout << componentR[j] << ", ";
            }
            //std::cout << "}" << std::endl;
        } */
        for (int j=0; j<(int) componentL.size(); j++)
        {
            setInactive(componentL[j]);
            L->push_back(componentL[j]);
        }
        for (int j=0; j<(int) componentR.size(); j++)
        {
            setInactive(componentR[j]);
            R->push_back(componentR[j]);
        }
        //std::cout << "}" << std::endl;
        if(componentL.size() != componentR.size()) { throw std::invalid_argument( "component sizes not equal!!!" ); }
        // Early stopping, when k doesnt allow for reduction
        if((int) R->size() > k) { return; }
        auto current = std::chrono::high_resolution_clock::now();
        if(std::chrono::duration_cast<std::chrono::microseconds>(current - start).count() / (double) 1000000 > maxExecTime) { return; }
    }
}

/*----------------------------------------------------------*/
/*------------------   Calculate Bounds   ------------------*/
/*----------------------------------------------------------*/

int BucketGraph::getLowerBoundVC() {
    return getCliqueBound();
    //return getLPBound();
    //return getLPCycleBound();
}

void BucketGraph::resetLPBoundDataStructures()
{
    initMatching();
}

int BucketGraph::getLPBound()
{
    hopcroftKarpMatchingSize();
    return currentLPBound/2;
}

int BucketGraph::getFlow()
{
    return edmondsKarpFlow();
}

// TODO: fix estimation
/* int BucketGraph::getLPCycleBound()
{
    auto matching = hopcroftKarpMatching();
    std::vector<int>* pairU = matching.second.first;
    std::vector<int>* pairV = matching.second.second;
    std::vector<int> covered = std::vector<int>(pairU->size());
    int LPCyclebound = 0;

    std::vector<int> currentCycle = std::vector<int>();
    bool foundCycle;
    int current;
    while(true)
    {
        // find uncovered vertex index
        current = -1;
        foundCycle = false;
        for(int i=0; i<(int) covered.size(); i++)
        {
            if(covered[i] != -1)
            {
                current = i;
                break;
            }
        }
        if(current == -1) break;
        //std::cout << "found uncovered vertex: " << current << "\n";

        bool onLeftSide = true;
        // determine cycle
        while(current != -1 && current != 0)
        {
            if(currentCycle.size() > 2 && currentCycle.front() == current)
            {
                foundCycle = true;
                break;
            }
            //std::cout << "adding " << current << " to cycle\n";
            currentCycle.push_back(current);
            covered[current] = -1;
            if(onLeftSide)
            {
                current = (*pairU)[current];
            }
            else
            {
                current = (*pairV)[current];
            }
            onLeftSide = !onLeftSide;
        }

        if(foundCycle)
        {
            // attempt to subdivide the cycle
            if(currentCycle.size() >= 6)
            {
                for(int c=1; (int) currentCycle.size() >= c + 4; c += 2) {
                    for(int i=0; i<(int) currentCycle.size(); i++)
                    {
                        vertexHasEdgeTo(currentCycle[i+1], currentCycle[i+2+c]);
                        if(vertexHasEdgeTo(currentCycle[i], currentCycle[i+3+c])
                        && vertexHasEdgeTo(currentCycle[i+1], currentCycle[i+2+c]))
                        {
                            LPCyclebound += 1 + c;
                            // TODO: cull i+1 -> i+2+c from currentCycle
                            currentCycle.erase(currentCycle.begin()+i+1, currentCycle.begin()+i+2+c+1);
                            i -= (2 + c); 
                            if((int) currentCycle.size() < c + 4) break;
                        }
                    }
                }
            }
            //std::cout << "found cycle of size: " << currentCycle.size() << "\n";
            // calculate
            LPCyclebound += std::ceil((double) currentCycle.size() / (double) 2);
        }
        else if(currentCycle.size() == 2)
        {
            LPCyclebound++;
        }
        currentCycle.clear();
    }
    return LPCyclebound;
} */

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
