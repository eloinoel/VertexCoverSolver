#include <iostream> //output streams
#include <fstream>  //ifstream file opening
#include <sstream>
#include <math.h>
#include <random>
#include <chrono>
using namespace std::chrono;

#include "BucketGraph.h"
#include "ColorPrint.h"
#include "Reductions.h"

typedef ColorPrint cp;

/*----------------------------------------------------------*/
/*-----------------   Graph Construction   -----------------*/
/*----------------------------------------------------------*/

BucketGraph* BucketGraph::readStandardInput(bool initReductionDataStructures)
{
    //init
    //auto startReadStandardInput = std::chrono::high_resolution_clock::now();
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
    //auto endReadStandardInput = std::chrono::high_resolution_clock::now();

    
    G->initVertexReferences(); // sets vertex references and the adj_maps
    G->initBucketQueue(); // initialise bucket queue

    if(initReductionDataStructures)
    {
        G->initMatching(); // LP Bound matching fields
        G->initUnconfined(); // Unconfined optimisation fields
        G->reductions = new Reductions();

        G->initDominationHelper(); // Domination
    }

//    G->reductions->initDominationVector(G);
//    G->reductions->printDominationSets();
    //auto endInit = std::chrono::high_resolution_clock::now();
    //double readStandardInput = (std::chrono::duration_cast<std::chrono::microseconds>(endReadStandardInput - startReadStandardInput).count() /  1000) / (double) 1000;
    //double initialization = (std::chrono::duration_cast<std::chrono::microseconds>(endInit - endReadStandardInput).count() /  1000) / (double) 1000;
    //std::cout << "Finished readStandardInput in " << readStandardInput + initialization << " seconds (read Input: " << readStandardInput << " + Other initializations: " << initialization << ")" << '\n';
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

void BucketGraph::initVertexReferences()
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
            if (vertexReferences[indexFirst]->index == indexFirst && vertexReferences[indexFirst]->adj_map != nullptr)
            {
                vertexReferences[indexFirst]->adj_map->insert({indexSecond, true});
            }
            else
            {
                throw std::invalid_argument("readStandardInput: inconsistency when trying to insert adjacent vertex to v1");
            }
        }
        else //edges list not yet initialised
        {
            Vertex* firstVertex = new Vertex(new std::unordered_map<int, bool>(), edges[k].first, indexFirst, maxDegFirst);
            vertexReferences[indexFirst] = firstVertex;
            //insert
            if (maxDegFirst > 0)
            {
                //(*firstVertex->adj).push_back(indexSecond);
                firstVertex->adj_map->insert({indexSecond, true});
            }
            else
            {
                throw std::invalid_argument("readStandardInput: trying to insert edge but maxDegFirst is 0");
            }
        }

        //second vertex: insert into adjacency list
        if (vertexReferences[indexSecond] != nullptr) //edges list initialised
        {
            if (vertexReferences[indexSecond]->index == indexSecond && vertexReferences[indexSecond]->adj_map != nullptr)
            {
                //vertexReferences[indexSecond]->adj->push_back(indexFirst);
                vertexReferences[indexSecond]->adj_map->insert({indexFirst, true});
            }
            else
            {
                throw std::invalid_argument("readStandardInput: inconsistency when trying to insert adjacent vertex to v1");
            }
        }
        else //edges list not yet initialised
        {
            Vertex* secondVertex = new Vertex(new std::unordered_map<int, bool>(), edges[k].second, indexSecond, maxDegSecond);
            vertexReferences[indexSecond] = secondVertex;
            //insert
            if (maxDegSecond > 0)
            {
                //(*secondVertex->adj).push_back(indexFirst);
                secondVertex->adj_map->insert({indexFirst, true});
            }
            else
            {
                throw std::invalid_argument("readStandardInput: trying to insert edge but maxDegFirst is 0");
            }
        }
    }
}


void BucketGraph::initBucketQueue()
{
    // init bucketQueue
    bucketQueue = list<Bucket>();
    bool inserted = false;
    int maxDeg = -1;
    // add active vertices to BucketQueue
    for (int i=0; i<(int) vertexReferences.size(); i++)
    {
        Vertex* vertex = vertexReferences[i];
        if(!isActive(vertex->index)) { throw std::invalid_argument("initBucketQueue: found inactive vertex, when no inactive vertices should exist"); }
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
                bucketQueue.insert(bucket, *(new Bucket(vertex->degree, std::vector<BucketVertex*>({vertex->bucketVertex}))));
                inserted = true;
                break;
            }
        }
        if(!inserted)
        {
            bucketQueue.insert(bucketQueue.end(), *(new Bucket(vertex->degree, std::vector<BucketVertex*>({vertex->bucketVertex}))));
        }
    }

    // init bucketReferences
    int deg = 0;
    bucketReferences = std::vector<Bucket*>(maxDeg+1);
    for (auto bucket = bucketQueue.begin(); bucket != bucketQueue.end(); ++bucket)
    {
        while(deg < bucket->degree)
        {
            bucketReferences[deg] = new Bucket(deg, std::vector<BucketVertex*>({}));
            deg++;
        }
        bucketReferences[bucket->degree] = &*bucket;
        deg++;
    }
    // init stable iterator
    stable_bucketQueue_inc_iterator = bucketQueue.iterator_to(*bucketQueue.begin());
    stable_bucketQueue_dec_iterator = bucketQueue.iterator_to(*bucketQueue.end());
}

void BucketGraph::initMatching()
{
    pairU = new std::vector<int>(vertexReferences.size());
    pairV = new std::vector<int>(vertexReferences.size());
    dist = new std::vector<int>(vertexReferences.size()+1);
    unmatched = new std::vector<int>();
    next_unmatched = new std::vector<int>();
    NIL = vertexReferences.size();
    for(int i=0; i<(int) pairU->size(); i++)
    {
        (*pairU)[i] = NIL;
        (*pairV)[i] = NIL;
        (*dist)[i] = INT32_MAX;
    }
    (*dist)[NIL] = INT32_MAX;
    // number of vertices helper variable (only for clarity)
    nv = vertexReferences.size();
    // run matching algorithm to calculate initial matching and setup flow
    hopcroftKarpMatchingSize();
    LP_INITIALISED = true;
}

void BucketGraph::initUnconfined()
{
    /* mayBeUnconfined = new std::vector<bool>(vertexReferences.size());
    for(int i=0; i<(int)mayBeUnconfined->size(); i++)
    {
        (*mayBeUnconfined)[i] = true;
    } */
    UNCONFINED_INITIALISED = true;
}

void BucketGraph::freeMatching()
{
    delete pairU;
    delete pairV;
    delete dist;
    delete unmatched;
    delete next_unmatched;
    NIL = -1;
    // number of vertices helper variable (only for clarity)
    nv = -1;
    LP_INITIALISED = false;
}

void BucketGraph::freeUnconfined()
{
    delete mayBeUnconfined;
    UNCONFINED_INITIALISED = false;
}

BucketGraph* BucketGraph::resetGraph()
{
    BucketGraph* G = new BucketGraph();
    G->originalVertexNames = originalVertexNames;
    G->edges = edges;

    reductions = new Reductions();
    G->initVertexReferences(); // sets vertex references and the adj_maps
    G->initBucketQueue();
    G->initMatching();
    return G;
}

void BucketGraph::freeGraph()
{
    //free bucket queue
    bucketQueue.clear();
    for(int i = 0; i < (int) bucketReferences.size(); ++i)
    {
        if(bucketReferences[i] != nullptr)
        {
            Bucket* bucket = bucketReferences[i];
            //TODO: maybe also free bucketVertex?
            bucket->vertices.clear();
            delete bucket;
            bucket = NULL;
        }
    }

    //free vertex references
    for(int i = 0; i < (int) vertexReferences.size(); ++i)
    {
        if(vertexReferences[i])
        {
            Vertex* vertex = vertexReferences[i];
            if(vertex->adj_map) { delete vertex->adj_map; vertex->adj_map = NULL; }
            if(vertex->bucketVertex) { delete vertex->bucketVertex; vertex->bucketVertex = NULL; }
            delete vertex;
            vertex = NULL;
        }
    }



    //free LP stuff
    if(LP_INITIALISED)
    {
        freeMatching();
    }

    //free Unconfined stuff
    if(UNCONFINED_INITIALISED)
    {
        freeUnconfined();
    }

    //free clique stuff

    //free reductions
    if(dominationHelper) { delete dominationHelper; dominationHelper = NULL; }
    if(reductions)
    {
        reductions->freeReductions();
        delete reductions;
        reductions = NULL;
    }
}

/* copies originalVertexNames and edges (deep copy) to a new graph and inits the graph from there*/
BucketGraph* BucketGraph::copy()
{
    BucketGraph* graphCopy = new BucketGraph();

    //std::cout << "#starting copy" << std::endl;
    graphCopy->originalVertexNames = std::unordered_map<std::string, std::pair<int, int>>();
    for(auto it = originalVertexNames.begin(); it != originalVertexNames.end(); ++it)
    {
        std::string originalName = it->first;
        graphCopy->originalVertexNames.insert({originalName, std::pair<int, int>({it->second.first, it->second.second})});
    }

    //std::cout << "#starting edges copy" << std::endl;

    graphCopy->edges = std::vector<std::pair<std::string, std::string>>((int) edges.size());
    for(int i = 0; i < (int) edges.size(); ++i)
    {
        std::string firstVertex = edges[i].first;
        std::string secondVertex = edges[i].second;
        graphCopy->edges[i] = std::pair<std::string, std::string>({firstVertex, secondVertex});
    }

    //std::cout << "#starting initVertexReferences" << std::endl;
    graphCopy->initVertexReferences(); // sets vertex references and the adj_maps
    //std::cout << "#starting initBucketQueue" << std::endl;
    graphCopy->initBucketQueue(); // initialise bucket queue
    //std::cout << "#starting initMatching" << std::endl;
    graphCopy->initMatching(); // LP Bound matching fields
    //std::cout << "#starting initUnconfined" << std::endl;
    graphCopy->initUnconfined(); // Unconfined optimisation fields
    graphCopy->reductions = new Reductions();
    //std::cout << "#returning" << std::endl;

    return graphCopy;
}

/*----------------------------------------------------------*/
/*-------------------   Graph Utility   --------------------*/
/*----------------------------------------------------------*/

void BucketGraph::print()
{
    if (vertexReferences.size() > 0)
	{
		std::cout << "\n";
        std::cout << "Graph of size: " << vertexReferences.size() << ", active: " << numVertices << ", edges: " << numEdges << '\n';
        std::cout << "BucketQueue of size: " << bucketQueue.size() << ", maxDegree: " << getMaxDegree() << '\n';
        std::cout << "name <name>, index <index>(<degree>): <neighbours>" << '\n';
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
                if(vertexReferences[i]->adj_map != nullptr)
                {
                    if(vertexReferences[i]->adj_map->empty())
                    {
                        std::cout << "<empty>";
                    }
                    else
                    {
                        std::vector<int> active = std::vector<int>();
                        std::vector<int> inactive = std::vector<int>();
                        for(auto it = vertexReferences[i]->adj_map->begin(); it != vertexReferences[i]->adj_map->end(); ++it)
                        {
                            if(vertexReferences[it->first]->isActive)
                            {
                                active.push_back(it->first);
                            }
                            else 
                            {
                                inactive.push_back(it->first);
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

void BucketGraph::printBucketQueue()
{
    int fieldSize = 4;
    int rowMin = 3;
    int rowMax = 7;
    int row;
    std::cout << "BucketQueue of size " << cp::dye(std::to_string(bucketQueue.size()), 'y') << " with maxDegree: " << getMaxDegree() << '\n';
    for(auto bucket = bucketQueue.begin(); bucket != bucketQueue.end(); ++bucket)
    {
        row = std::min(rowMax, std::max((int) bucket->vertices.size(), rowMin));
        std::cout << "/-" << cp::dye("Bucket " + std::to_string(bucket->degree), 'w') << cp::repeatStr("-", (2+fieldSize)*row-8-cp::digitCount(bucket->degree)) << "\\" << '\n';
        for (auto vertex = bucket->vertices.begin(); vertex != bucket->vertices.end();)
        {
            int offs = 0;

            std::cout << "| " << cp::dye(std::to_string(vertex->index), 'w');
            offs = fieldSize-cp::digitCount(vertex->index);
            ++vertex;
            int i = 1;
            for (; i < row && vertex != bucket->vertices.end(); i++, ++vertex)
            {
                std::cout << ", " << cp::repeatStr(" ", offs) << cp::dye(std::to_string(vertex->index), 'w');
                offs = fieldSize-cp::digitCount(vertex->index);
            }
            std::cout << cp::repeatStr(" ", (2+fieldSize)*(row-i));
            std::cout << cp::repeatStr(" ", offs) << " |" << '\n';
        }
        std::cout << "\\" << cp::repeatStr("-", (2+fieldSize)*row) << "/" << '\n' << '\n';
    }
}

void BucketGraph::printBucketSizes()
{
    std::cout << "BucketQueue of size " << cp::dye(std::to_string(bucketQueue.size()), 'y') << " with maxDegree: " << getMaxDegree() << '\n';
    for(auto bucket = bucketQueue.begin(); bucket != bucketQueue.end(); ++bucket)
    {
        std::cout << cp::dye("Bucket " + std::to_string(bucket->degree) + " size: " + std::to_string(bucket->vertices.size()), 'w') << '\n';
    }
    int deg3Reduc = 0;
    for(auto vertex = getBucket(3)->vertices.begin(); vertex != getBucket(3)->vertices.end(); ++vertex)
    {
        int v = vertex->index;
        auto ns = getNeighbours(v);
        int v1 = ns->at(0);
        int v2 = ns->at(1);
        int v3 = ns->at(2);
        if((vertexHasEdgeTo(v1, v2) && !vertexHasEdgeTo(v1, v3) && !vertexHasEdgeTo(v2, v3))
        || (vertexHasEdgeTo(v1, v3) && !vertexHasEdgeTo(v1, v2) && !vertexHasEdgeTo(v2, v3))
        || (vertexHasEdgeTo(v2, v3) && !vertexHasEdgeTo(v1, v2) && !vertexHasEdgeTo(v1, v3))) {
            deg3Reduc++;
        }
    }
    std::cout << "Applicable: " << deg3Reduc << '\n';
}

void BucketGraph::printMatching()
{
    std::cout << "Current Matchings in the Graph: " << currentLPBound << '\n';
    std::cout << "----- Left Matchings: -----" << '\n';
    for(int i=0; i<(int) pairU->size(); i++) {
        if((*pairU)[i] != NIL) {
            std::cout << i << " >> " << (*pairU)[i] << '\n';
        }
    }
    std::cout << "----- Right Matchings: -----" << '\n';
    for(int i=0; i<(int) pairV->size(); i++) {
        if((*pairV)[i] != NIL) {
            std::cout << (*pairV)[i] << " << " << i << '\n';
        }
    }
}

void BucketGraph::printVC(std::unordered_map<int, bool>* vc) {
    if(vc->size() == 0) {
        std::cout << "{}" << std::endl;
        return;
    }
    std::cout << "{" << cp::dye(std::to_string(vc->begin()->first), 'y');
    for(auto vertex=vc->begin(); vertex!=vc->end(); ++vertex)
    {
        if(vertex == vc->begin()) { continue; }
        std::cout << ", " + cp::dye(std::to_string(vertex->first), 'y');
    }
    std::cout << "}" << std::endl;
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
            for(auto it = v->adj_map->begin(); it != v->adj_map->end(); ++it) 
            {
                if(isActive(it->first))
                {
                    auto map_entry = edges_map->find(std::make_pair(v->index, it->first));
                    if(map_entry != edges_map->end()) continue;
                    map_entry = edges_map->find(std::make_pair(it->first, v->index));
                    if(map_entry != edges_map->end()) continue;

                    //else insert edge
                    (*edges_map)[std::make_pair(v->index, it->first)] = true;
                    //std::cout << it->strName << " " << vertexReferences[it->adj->at(i)]->strName << '\n';
                    str->push_back(v->strName + " " + vertexReferences[it->first]->strName + "\n");
                    //str += v->strName + " " + vertexReferences[v->adj->at(i)]->strName + "\n";
                }
            }
            /* for(int i = 0; i < (int) v->adj->size(); i++)
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
                    //std::cout << it->strName << " " << vertexReferences[it->adj->at(i)]->strName << '\n';
                    str->push_back(v->strName + " " + vertexReferences[v->adj->at(i)]->strName + "\n");
                    //str += v->strName + " " + vertexReferences[v->adj->at(i)]->strName + "\n";
                }
            } */
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
        //std::cout << edges.at(i).first << " " << edges.at(i).second << '\n';
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

void BucketGraph::printVertices(std::vector<int>* vertices)
{
    if(vertices == nullptr)
        throw std::invalid_argument("printVertices: nullptr");
        
    for (int i = 0; i < (int) vertices->size(); ++i)
    {
        std::cout << vertexReferences[vertices->at(i)]->strName << "\n";
    }
}

void BucketGraph::printVertices(std::unordered_map<int, bool>* vertices)
{
    if(vertices == nullptr)
        throw std::invalid_argument("printVertices: nullptr");

    for (auto it = vertices->begin(); it != vertices->end(); ++it)
    {
        std::cout << vertexReferences[it->first]->strName << "\n";
    }
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
    //std::cout << "removeFromBucketQueue: before bucket removal" << std::endl;
    bucketReferences[degree]->remove(vertexReferences[vertex]->bucketVertex);
    //std::cout << "removeFromBucketQueue: after bucket removal" << std::endl;
    if(bucketReferences[degree]->vertices.size() == 0) {
        //std::cout << "removeFromBucketQueue: before iterator updates" << std::endl;
        if(degree == stable_bucketQueue_inc_iterator->degree) { ++stable_bucketQueue_inc_iterator; }
        if(degree == stable_bucketQueue_dec_iterator->degree) { --stable_bucketQueue_dec_iterator; }
        //std::cout << "removeFromBucketQueue: before erase" << std::endl;
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
        if(previous_bucket_it->degree == stable_bucketQueue_inc_iterator->degree) { ++stable_bucketQueue_inc_iterator; }
        if(previous_bucket_it->degree == stable_bucketQueue_dec_iterator->degree) { --stable_bucketQueue_dec_iterator; }
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
    //std::cout << "moveToSmallerBucket: deg = " << degree << '\n';
    //printBucketQueue();
    //remove from current bucket
    Bucket* curBucket = bucketReferences[degree];
    curBucket->remove(vertexReferences[vertex]->bucketVertex);
    //std::cout << "after removeFromBucket: deg = " << degree << '\n';

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
        if(cur_bucket_it->degree == stable_bucketQueue_inc_iterator->degree) { ++stable_bucketQueue_inc_iterator; }
        if(cur_bucket_it->degree == stable_bucketQueue_dec_iterator->degree) { --stable_bucketQueue_dec_iterator; }
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
    if(v == nullptr || v->adj_map == nullptr) throw std::invalid_argument("setActive: vertex " + std::to_string(vertexIndex) + " is nullptr");
    if(v->isActive) throw std::invalid_argument("setActive: vertex " + std::to_string(vertexIndex) + " is already active");

    v->isActive = true;
    numVertices++;

    addToBucketQueue(v->index); //TODO: optimise this to not loop by storing previous references in recursion
    //update degree of all adjacent nodes
    for(auto it = v->adj_map->begin(); it != v->adj_map->end(); ++it)
    {
        if(vertexReferences[it->first] == nullptr) throw std::invalid_argument("setActive: tried updating degree of vertex pointing to nullptr");

        vertexReferences[it->first]->degree++;
        if(!vertexReferences[it->first]->isActive) { continue; }
        numEdges++;
        moveToBiggerBucket(vertexReferences[it->first]->degree - 1, it->first);
    }
    if (!LP_INITIALISED) { return; }
    queueForMatching(vertexIndex);
    //std::cout << "Restored vertex: " << vertexIndex << '\n';
}

void BucketGraph::setInactive(int vertexIndex)
{
    Vertex* v = vertexReferences[vertexIndex];
    if(v == nullptr || v->adj_map == nullptr) throw std::invalid_argument("setInactive: vertex " + std::to_string(vertexIndex) + " is nullptr");
    if(!v->isActive) throw std::invalid_argument("setInactive: vertex " + std::to_string(vertexIndex) + " is already not active");

    v->isActive = false;
    numVertices--;

    //std::cout << "before removing from queue: " << v->index << '\n';
    removeFromBucketQueue(v->index);
    //update degree of all adjacent nodes
    for(auto it = v->adj_map->begin(); it != v->adj_map->end(); ++it)
    {
        if(vertexReferences[it->first] == nullptr) throw std::invalid_argument("setInactive: tried updating degree of vertex pointing to nullptr");

        vertexReferences[it->first]->degree--;
        if(!vertexReferences[it->first]->isActive) { continue; }
        numEdges--;
        //std::cout << "before moving neighbour " << (*v->adj)[i] << " to smaller bucket" << '\n';
        moveToSmallerBucket(vertexReferences[it->first]->degree + 1, it->first);
        //std::cout << "after moving neighbour " << (*v->adj)[i] << " to smaller bucket" << '\n';
    }

    // update matching
    if(LP_INITIALISED)
    {
        removeFromMatching(vertexIndex);
    }

    if(UNCONFINED_INITIALISED)
    {
        // This is a O(m+n) runtime risk // TODO: so remove probably
        //scheduleComponentForUnconfined(v->index);
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
    Vertex* v = vertexReferences[vertexIndex];
    std::vector<int>* neighbours = new std::vector<int>();
    if(v == nullptr)
    {
        throw std::invalid_argument("getNeighbours: vertex " + std::to_string(vertexIndex) + " is nullptr");
    }

    for(auto it = v->adj_map->begin(); it != v->adj_map->end(); ++it)
    {
        if(vertexReferences[it->first] == nullptr) { throw std::invalid_argument("getNeighbours: vertex pointing to nullptr"); }
        if(vertexReferences[it->first]->isActive)
        {
            neighbours->push_back(it->first);
        }
    }
/*     for(int i = 0; i < (int) v->adj->size(); i++)
    {
        if(vertexReferences[i] == nullptr)
        {
            throw std::invalid_argument("getNeighbours: tried updating degree of vertex pointing to nullptr");
        }
        if(vertexReferences[(*v->adj)[i]]->isActive)
        {
            neighbours->push_back((*v->adj)[i]);
        }
    } */
    return neighbours;
}

bool BucketGraph::containsConnectedVertex()
{
    return (bucketQueue.back().degree != 0);
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


int BucketGraph::getRandomConnectedVertex(int randomRangeCap)
{
    //std::cout << "bucketQueue size: " << bucketQueue.size() << std::endl;
    if(bucketQueue.empty())
        return -1;

    int minBucket = 0;
    int numBuckets = bucketQueue.size();
    int smallestBucketDegree = bucketQueue.front().degree;
    if (smallestBucketDegree == 0) { minBucket = 1; }
    if (bucketQueue.size() == 1 && smallestBucketDegree == 0) { return -1; }

    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> bucketDistr(minBucket, numBuckets-1); // define the range

    int randomBucketIndex = (int) bucketDistr(gen);
    auto bucketIt = bucketQueue.begin();
    for(int i = 0; i < randomBucketIndex ; ++i)
    {
        ++bucketIt;
    }
    //std::cout << "chose bucket of degree: " << bucketIt->degree << std::endl;
    list<BucketVertex>* bucket = &(bucketIt->vertices);
    if(bucket == nullptr || bucket->size() == 0)
        throw std::invalid_argument("getRandomConnectedVertex: random bucket doesn't exist or is empty");

    auto startRandomNrCalc = std::chrono::high_resolution_clock::now();
    //choose random vertex
    int maxRange = bucket->size() - 1;
    if(randomRangeCap > 0 && randomRangeCap < maxRange)
        maxRange = randomRangeCap;

    std::uniform_int_distribution<> distr(0, maxRange); // define the range
    int randomIndex = (int) distr(gen);
    auto endRandomNrCalc = std::chrono::high_resolution_clock::now();

    auto startIteratorCalc = std::chrono::high_resolution_clock::now();
    auto it = bucket->begin();
    for(int i = 0; i < randomIndex; ++i)
    {
        ++it;
    }
    auto endIteratorCalc = std::chrono::high_resolution_clock::now();

    double randomNrCalc = (std::chrono::duration_cast<std::chrono::microseconds>(endRandomNrCalc - startRandomNrCalc).count());
    double iteratorCalc = (std::chrono::duration_cast<std::chrono::microseconds>(endIteratorCalc - startIteratorCalc).count());
    //std::cout << "randomNrCalc: " << randomNrCalc << ", iteratorCalc: " << iteratorCalc << '\n';

    return it->index;

}

int BucketGraph::getRandomMaxDegreeVertex(int randomRangeCap)
{
    if(bucketQueue.empty())
        return -1;

    int maxDegree = bucketQueue.back().degree;
    list<BucketVertex>* bucket = getVerticesOfDegree(maxDegree);
    
    if(bucket == nullptr || bucket->size() == 0)
        throw std::invalid_argument("getRandomMaxDegreeVertex: max degree bucket doesn't exist or is empty");

    if(maxDegree == 0) //no need to randomise here
        return bucket->begin()->index;

    auto startRandomNrCalc = std::chrono::high_resolution_clock::now();
    //choose random vertex
    int maxRange = bucket->size() - 1;
    if(randomRangeCap > 0 && randomRangeCap < maxRange)
        maxRange = randomRangeCap;
 
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, maxRange); // define the range
    int randomIndex = (int) distr(gen);
    auto endRandomNrCalc = std::chrono::high_resolution_clock::now();

    auto startIteratorCalc = std::chrono::high_resolution_clock::now();
    auto it = bucket->begin();
    for(int i = 0; i < randomIndex; ++i)
    {
        ++it;
    }
    auto endIteratorCalc = std::chrono::high_resolution_clock::now();

    double randomNrCalc = (std::chrono::duration_cast<std::chrono::microseconds>(endRandomNrCalc - startRandomNrCalc).count());
    double iteratorCalc = (std::chrono::duration_cast<std::chrono::microseconds>(endIteratorCalc - startIteratorCalc).count());
    //std::cout << "randomNrCalc: " << randomNrCalc << ", iteratorCalc: " << iteratorCalc << '\n';

    return it->index;

}

/* returns min degree vertex of degree > 0 and -1 if doesn't exist */
int BucketGraph::getMinDegreeVertex()
{
    if(bucketQueue.empty())
        return -1;
    
    //return a vertex of degree > 0
    if(!containsConnectedVertex())
        return -1;

    auto smallestDegreeBucket = bucketQueue.begin();
    if(bucketQueue.size() < 2 && smallestDegreeBucket->degree == 0)
        throw std::invalid_argument("getMinDegreeVertex: Inconsistency error because bucketQueue does not include enough buckets");

    if(smallestDegreeBucket->degree == 0)
    {
        ++smallestDegreeBucket;
    }

    if(smallestDegreeBucket->vertices.size() == 0)
        throw std::invalid_argument("getMinDegreeVertex: Inconsistency error because min degree bucket is empty");

    return smallestDegreeBucket->vertices.back().index;
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
    for(auto it = v->adj_map->begin(); it != v->adj_map->end(); ++it)
    {
        if(vertexReferences[it->first]->isActive)
        {
            if(counter == n)
            {
                return it->first;
            }
            counter++;
        }
    }
    return -1;
}

std::pair<int, int> BucketGraph::getFirstTwoActiveNeighbours(int vertex)
{
    std::pair<int, int> neighbours = std::pair<int, int>({-1, -1});

    if(vertex >= (int) vertexReferences.size())
    {
        throw std::invalid_argument("getFirstActiveNeighbour: vertex");
    }
    Vertex* v = vertexReferences[vertex];
    int counter = 0;
    for(auto it = v->adj_map->begin(); it != v->adj_map->end(); ++it)
    {
        if(vertexReferences[it->first]->isActive)
        {
            if(counter == 0)
            {
                neighbours.first = it->first;
            }
            else if (counter == 1)
            {
                neighbours.second = it->first;
                return neighbours;
            }
            counter++;
        }
    }
    return neighbours;
}

int BucketGraph::getNumConnectedVertices()
{
    if(getNumVertices()>0)
    {
        return numVertices-getVerticesOfDegree(0)->size();
    }
    return numVertices;
}

int BucketGraph::getTotalNumVertices()
{
    return vertexReferences.size();
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
void BucketGraph::preprocess(int* k, bool printDebug)
{
    while(true)
    {
        if(reductions->rule_DegreeOne(this, k, -1, false, printDebug) == APPLICABLE) continue;
        if(reductions->rule_DegreeTwo(this, k, -1, false, printDebug) == APPLICABLE) continue;

        //TODO: add depth to rules
        if(reductions->rule_DegreeThree_Domination(this, k, -1, deg3domDeg2, false, printDebug) == APPLICABLE) continue;
        if(reductions->rule_DegreeThree_Clique(this, k, -1, false, printDebug) == APPLICABLE){ continue; }
        if(reductions->rule_DegreeThree_Independent(this, -1, printDebug) == APPLICABLE) continue;
        //TODO:
        //if(reductions->rule_Domination(this, k, -1, false) == APPLICABLE) continue;
        if(reductions->rule_Unconfined(this, k, -1, false, printDebug) == APPLICABLE) continue;
        if(reductions->rule_LPFlow(this, k, -1, false, printDebug) == APPLICABLE) continue;

        return;
    }
}

/*
 * 0: deg1,
 * 1: deg2,
 * 2: domination,
 * 3: unconfined,
 * 4: LP,
 * 5: highDeg (not in preprocessing),
 * 6: Buss (not in preprocessing),
 * 7: Deg3IndepentSet,
 * 8: Deg3Clique,
 * 9: Deg3Domination
*/
void BucketGraph::preprocess(int* k, std::vector<bool>& rulesToApply, bool printDebug)
{
    while(true)
    {
        if(rulesToApply[0] && reductions->rule_DegreeOne(this, k, -1, false, printDebug) == APPLICABLE) continue;
//        if(rulesToApply[1] && reductions->rule_DegreeTwo_Simple_Case(this, k, -1, false, printDebug) == APPLICABLE) continue;
        if(rulesToApply[1] && reductions->rule_DegreeTwo(this, k, -1, false, printDebug) == APPLICABLE) continue;
        //TODO:
        if(rulesToApply[9] && reductions->rule_DegreeThree_Domination(this, k, -1, deg3domDeg2, false, printDebug) == APPLICABLE) continue;
        if(rulesToApply[8] && reductions->rule_DegreeThree_Clique(this, k, -1, false, printDebug) == APPLICABLE){ continue; }
        if(rulesToApply[7] && reductions->rule_DegreeThree_Independent(this, -1, printDebug) == APPLICABLE) continue;

        if(rulesToApply[10] && reductions->rule_DegreeFour_Clique(this, k, -1, false, printDebug) == APPLICABLE) continue;
        if(rulesToApply[11] && reductions->rule_DegreeFour_Domination(this, k, -1, false, printDebug) == APPLICABLE) continue;
        //TODO:
        if(rulesToApply[2] && reductions->rule_Domination(this, k, -1, false) == APPLICABLE) continue;
        if(rulesToApply[3] && reductions->rule_Unconfined(this, k, -1, false, printDebug) == APPLICABLE) continue;
        if(rulesToApply[4] && reductions->rule_LPFlow(this, k, -1, false, printDebug) == APPLICABLE) continue;
        return;
    }
}

void BucketGraph::printReductionStack() { reductions->printReductionStack(); }
/*
 * 0: deg1, 1: deg2, 2: domination, 3: unconfined, 4: LP, 5: highDeg, 6: buss
*/
bool BucketGraph::dynamicReduce(int* k, int depth, bool printDebug)
{
                                                    //  0     1     2     3     4     5     6      7    8     9
    std::vector<bool> reductions = std::vector<bool>{true, true, false, false, false, true, true, false, false, true, false, true};

    if(depth % period_deg3 == 0){
        reductions.at(7) = true; // Deg4 Ind
        reductions.at(8) = true; // Deg4 2-Clique
        if(period_deg3 > 1) period_deg3--;
    }
    if(depth % period_unc == 0)
    {
////         + unconfined
        reductions.at(4) = UNCONFINED_INITIALISED;
//        reductions.at(10) = true; // Deg4 2-Clique
        if(period_unc > 1) period_unc--;

    }
    if(depth % period_lp == 0)
    {
        // + LP
//        reductions.at(7) = true; // Deg4 Ind
//        reductions.at(8) = true; // Deg4 2-Clique
        reductions.at(4) = LP_INITIALISED;
        reductions.at(10) = true; // Deg4 2-Clique
        if(period_lp > 1) period_lp--;
    }

    return reduce(k, depth, &reductions, printDebug);
}

/*
 * 0: deg1,
 * 1: deg2,
 * 2: domination,
 * 3: unconfined,
 * 4: LP,
 * 5: highDeg,
 * 6: Buss
 * 7: Deg3IndepentSet,
 * 8: Deg3Clique,
 * 9: Deg3Domination
*/
bool BucketGraph::reduce(int* k, int depth, std::vector<bool>* rulesToApply, bool printDebug)
{
    if(rulesToApply == nullptr) {
                                    //  0     1       2       3    4      5      6   7       8   9
        rulesToApply = new std::vector{true, true, false, true, true, true, true, true, false, false};
//        rulesToApply = new std::vector{true, true, true, true, true, true, true, true};
    }
    //initialisise
    RULE_APPLICATION_RESULT degreeTwoResult = INAPPLICABLE;
    RULE_APPLICATION_RESULT highDegreeResult = INAPPLICABLE;
    RULE_APPLICATION_RESULT degreeZeroResult = INAPPLICABLE;
    RULE_APPLICATION_RESULT degreeOneResult = INAPPLICABLE;
    RULE_APPLICATION_RESULT dominationResult = INAPPLICABLE;
    RULE_APPLICATION_RESULT unconfinedResult = INAPPLICABLE;
    RULE_APPLICATION_RESULT LPFlowResult = INAPPLICABLE;
    RULE_APPLICATION_RESULT degreeThreeIndResult = INAPPLICABLE;
    RULE_APPLICATION_RESULT degreeThreeCliqResult = INAPPLICABLE;
    RULE_APPLICATION_RESULT degreeThreeDomResult = INAPPLICABLE;
    RULE_APPLICATION_RESULT degreeFourCliqResult = INAPPLICABLE;
    RULE_APPLICATION_RESULT degreeFourDomResult = INAPPLICABLE;

    while(true)
    {
//        std::cout << "highdeg " << '\n';
        if(rulesToApply->at(5)) { highDegreeResult = reductions->rule_HighDegree(this, k, depth); }
        if(highDegreeResult == APPLICABLE) continue;
        if(highDegreeResult == INSUFFICIENT_BUDGET) return true;

        if(rulesToApply->at(6) && reductions->rule_Buss(this, k, getNumConnectedVertices(), getNumEdges()) == APPLICABLE) return true;

        if(rulesToApply->at(0)) { degreeOneResult = reductions->rule_DegreeOne(this, k, depth, true, printDebug); }
        if(degreeOneResult == APPLICABLE) continue;
        if(degreeOneResult == INSUFFICIENT_BUDGET) return true;

        //std::cout << "deg2 " << '\n';
//        if(rulesToApply->at(1)) { degreeTwoResult = reductions->rule_DegreeTwo_Simple_Case(this, k, depth, true, printDebug); }
        if(rulesToApply->at(1)) { degreeTwoResult = reductions->rule_DegreeTwo(this, k, depth, true, printDebug ); }
        if(degreeTwoResult == APPLICABLE) continue;
        if(degreeTwoResult == INSUFFICIENT_BUDGET) return true;

        if(rulesToApply->at(2)) { dominationResult = reductions->rule_Domination(this, k, depth, true); }

        if (rulesToApply->at(9)) { degreeThreeDomResult = reductions->rule_DegreeThree_Domination(this, k, depth, false, true, deg3dom);}
        if (degreeThreeDomResult == INSUFFICIENT_BUDGET) return true;
        if (degreeThreeDomResult == APPLICABLE) continue;
        if (rulesToApply->at(8)) { degreeThreeCliqResult = reductions->rule_DegreeThree_Clique(this, k, depth, true, printDebug); }
        if (degreeThreeCliqResult == INSUFFICIENT_BUDGET) return true;
        if (degreeThreeCliqResult == APPLICABLE) continue;
        if (rulesToApply->at(7)) { degreeThreeIndResult = reductions->rule_DegreeThree_Independent(this, depth, printDebug); }
        if (degreeThreeIndResult == INSUFFICIENT_BUDGET) return true;
        if (degreeThreeIndResult == APPLICABLE) continue;

        if (rulesToApply->at(10)) { degreeFourCliqResult = reductions->rule_DegreeFour_Clique(this, k, depth, true, printDebug); }
        if (degreeFourCliqResult == INSUFFICIENT_BUDGET) return true;
        if (degreeFourCliqResult == APPLICABLE) continue;
        if (rulesToApply->at(11)) { degreeFourDomResult = reductions->rule_DegreeFour_Domination(this, k, depth, true, printDebug); }
        if (degreeFourDomResult == INSUFFICIENT_BUDGET) return true;
        if (degreeFourDomResult == APPLICABLE) continue;

//        if(rulesToApply->at(7) && rulesToApply->at(8) && rulesToApply->at(9) && getBucketSize(3) != 0)
//            throw std::invalid_argument("reduce error: There shouldn't be any degree 3 vertices left");

        if(rulesToApply->at(2)) { dominationResult = reductions->rule_Domination(this, k, depth, true); }
        if(dominationResult == APPLICABLE) continue;
        if(dominationResult == INSUFFICIENT_BUDGET) return true;

        if(rulesToApply->at(3)) { unconfinedResult = reductions->rule_Unconfined(this, k, depth, true, printDebug); }
        if(unconfinedResult == APPLICABLE) continue;
        if(unconfinedResult == INSUFFICIENT_BUDGET) return true;

        if(rulesToApply->at(4)) { LPFlowResult = reductions->rule_LPFlow(this, k, depth, true, printDebug); }
        if(LPFlowResult == APPLICABLE) continue;
        if(LPFlowResult == INSUFFICIENT_BUDGET) return true;

        return false;
    }
    delete rulesToApply;
    return false;
}

void BucketGraph::unreduce(int* k, int previousK, int depth, std::unordered_map<int, bool>* vc)
{
    if(reductions->appliedRules == nullptr)
        throw std::invalid_argument("unreduce: appliedRules is nullptr");
    if(reductions->appliedRules->empty())
        return;


    //pop rules
    while(!reductions->appliedRules->empty() && (*k < previousK || reductions->appliedRules->back()->kDecrement == 0))
    {
        Reduction* rule = reductions->appliedRules->back();
        //std::cout << "unreduce: rule: " << rule->rule << ", depth: " << depth << ", rule->depth: " << rule->rDepth << std::endl;
        if(rule->rDepth != depth)
        {
            if(depth < rule->rDepth)
            {
                std::cout << "unreduce: rule: " << rule->rule << ", depth: " << depth << ", rule->depth: " << rule->rDepth << "\n";
                throw std::invalid_argument("unreduce: attempted to unreduce rules of incorrect depth, unreducing is inconsistent with reducing");
            }
            break;
        }
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
                if(rule->mergeVertexInfo == nullptr) //didn't merge
                {
                    //std::cout << "deg2 non-merged" << std::endl;
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
                    //std::cout << "before unmerge" << std::endl;
                    int mergeVertex = std::get<0>(*rule->mergeVertexInfo);
                    //std::cout << "deg2 merged" << std::endl;
                    unmerge(rule); //handles setting vertices back active
                    //std::cout << "deg2 post unmerge" << std::endl;
                    if(vc != nullptr)
                    {
                        //if vc contains vertex that was merged into, add neighbours of deg2 vertex to the solution
                        auto it = vc->find(mergeVertex);
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
                //std::cout << cp::dye("Domination unreduce", 'g') << '\n';
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
                
                    std::cout << "}" << '\n';
                }
                else
                {
                    std::cout << "nullptr" << '\n';
                } */
                *k = *k + rule->kDecrement;
                setActive(rule->deletedVCVertices);
                //std::cout << "after set active" << '\n';
                if(vc != nullptr)
                {
                    for(int i = 0; i < (int) rule->deletedVCVertices->size(); i++)
                    {
                        vc->insert({rule->deletedVCVertices->at(i), true});
                    }
                }
                //std::cout << "after unreduce" << '\n';
                //print();
                break;
            case UNCONFINED:
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
                std::cout << "}" << '\n'; */
                if(vc != nullptr)
                {
                    for(int i = 0; i < (int) rule->deletedVCVertices->size(); i++)
                    {
                        vc->insert({rule->deletedVCVertices->at(i), true});
                    }
                    //vc->insert(vc->end(), rule->deletedVCVertices->begin(), rule->deletedVCVertices->end());
                }
                break;
            case DEGREE_THREE_IND: {
                if (depth != rule->rDepth)
                    return;

                if(deg3ind) {
                    std::cout << "\n";
                    std::cout << "Recursion depth coincide => unreduce rule: Degree 3: Independent Set\n";
                    std::cout << depth << " == " << rule->rDepth << '\n';
                }

                unreduceDeg3Ind(rule, vc);

                break;
            }
            case DEGREE_THREE_CLIQ:
            {
                if (depth != rule->rDepth)
                    return;

                if(deg3clique) {
                    std::cout << "\n";
                    std::cout << "Recursion depth coincide => unreduce rule: Degree 3: 2-Clique-Neighbourhood\n";
                    std::cout << depth << " == " << rule->rDepth << '\n';
                }
                *k = *k + rule->kDecrement;

                unreduceDeg3Clique(rule, vc);

                break;
            }
            case DEGREE_THREE_DOM:
            {
                if (depth != rule->rDepth)
                    return;

                if(deg3dom) {
                    std::cout << "\n";
                    std::cout << "Recursion depth coincide => unreduce rule: Degree 3: Domination\n";
                    std::cout << depth << " == " << rule->rDepth << '\n';
                }
                *k = *k + rule->kDecrement;

                unreduceDeg3Dom(rule, vc);

                break;
            }
            case DEGREE_FOUR_CLIQUE:
            {
                if (depth != rule->rDepth)
                    return;

                if(deg4clique) {
                    std::cout << "\n";
                    std::cout << "Recursion depth coincide => unreduce rule: Degree 4: Clique\n";
                    std::cout << depth << " == " << rule->rDepth << '\n';
                }
                if(deg4clique)
                    std::cout << "UNREDUCE: Degree 4: Clique\n";
                unreduceDeg4Clique(rule, k, vc);

                break;
            }
            case DEGREE_FOUR_DOM:
            {
                if (depth != rule->rDepth)
                    return;

                if(deg4dom) {
                    std::cout << "\n";
                    std::cout << "Recursion depth coincide => unreduce rule: Degree 4: Domination\n";
                    std::cout << depth << " == " << rule->rDepth << '\n';
                }
                if(deg4dom)
                    std::cout << "UNREDUCE: Degree 4: Domination\n";
                unreduceDeg4Dom(rule, k, vc);

                break;
            }
            default:
                throw std::invalid_argument("unreduce error: unknown rule");
                break;
        }
        //std::cout << "unreduce: before delete" << std::endl;
        if (rule->deletedVCVertices != nullptr) delete rule->deletedVCVertices;
        if (rule->deletedVertices != nullptr) delete rule->deletedVertices;
        if(rule->addedEdges != nullptr) delete rule->addedEdges;
        reductions->appliedRules->pop_back();
        delete rule;
        if(*k > previousK)
        {
            throw std::invalid_argument("unreduce error: " + std::to_string(*k) + " > " + std::to_string(previousK) + ", stop coding garbage");
        }
    }
}

std::tuple<int, std::unordered_map<int, bool>*, std::vector<int>*>* BucketGraph::merge(int v0, int v1, int v2)
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
    // TODO:
    //std::cout << cp::dye("merging into vertex " + std::to_string(maxDegVertex) + " <-- " + std::to_string(mergev0) + ", " + std::to_string(mergev1), 'p') << '\n';
    //if(maxDegVertex == 21) { print(); }
    /* std::cout << "Neighbour " << mergev0 << " (" << getVertexDegree(mergev0) << ") " << " neighbours: ";
    std::vector<int>* neighboursOfFirst = getNeighbours(mergev0);
    for(auto ne = neighboursOfFirst->begin(); ne != neighboursOfFirst->end(); ++ne)
    {
        std::cout << *ne << ", ";
    }
    std::cout << '\n';
    std::cout << "Neighbour " << mergev1 << " with deg=" << getVertexDegree(mergev1) << " neighbours: ";
    neighboursOfFirst = getNeighbours(mergev1);
    for(auto ne = neighboursOfFirst->begin(); ne != neighboursOfFirst->end(); ++ne)
    {
        std::cout << *ne << ", ";
    }
    std::cout << '\n'; */
    
    /* if(mergev1 == 18)
    {
        print();
    } */
    //std::cout << "merge: maxDegv: " << maxDegVertex << ", mergev0: " << mergev0 << ", mergev1: " << mergev1 << '\n';

    //std::cout << "merge_debug1.5" << '\n';

    //save adjacency list of vertex to merge into
    std::unordered_map<int, bool>* adj_map_copy = new std::unordered_map<int, bool>(*vertexReferences[maxDegVertex]->adj_map);
    std::vector<int>* added_vertices = new std::vector<int>();

    //std::cout << "changing adj lists " << '\n';
    removeFromBucketQueue(maxDegVertex);
    //add edges of merging vertices to adj of vertex to merge into
    for(auto it = vertexReferences[mergev0]->adj_map->begin(); it != vertexReferences[mergev0]->adj_map->end(); ++it)
    {
        int neighbour = it->first;
        if(isActive(neighbour) && neighbour != maxDegVertex && neighbour != mergev1 && !vertexHasEdgeTo(maxDegVertex, neighbour))
        {
            numEdges++;
            //std::cout << "added " << neighbour << '\n';
            vertexReferences[maxDegVertex]->adj_map->insert({neighbour, true});
            //moveToBiggerBucket(vertexReferences[maxDegVertex]->degree, maxDegVertex);
            vertexReferences[maxDegVertex]->degree++;

            //also change adj lists of neighbour
            vertexReferences[neighbour]->adj_map->insert({maxDegVertex, true});
            vertexReferences[neighbour]->degree++; //because of previously setting inactive, degree is -1 and needs to be incremented, if edge to merge vertex is inserted
            moveToBiggerBucket(vertexReferences[neighbour]->degree-1, neighbour);

            //save for unmerging later
            added_vertices->push_back(neighbour);
        }
    }
    for(auto it = vertexReferences[mergev1]->adj_map->begin(); it != vertexReferences[mergev1]->adj_map->end(); ++it)
    {
        int neighbour = it->first;
        if(isActive(neighbour) && neighbour != maxDegVertex && neighbour != mergev0 && !vertexHasEdgeTo(maxDegVertex, neighbour))
        {
            numEdges++;
            //std::cout << "added " << neighbour << '\n';
            vertexReferences[maxDegVertex]->adj_map->insert({neighbour, true});
            //moveToBiggerBucket(vertexReferences[maxDegVertex]->degree, maxDegVertex);
            vertexReferences[maxDegVertex]->degree++;

            //also change adj lists of neighbour
            vertexReferences[neighbour]->adj_map->insert({maxDegVertex, true});
            vertexReferences[neighbour]->degree++;
            moveToBiggerBucket(vertexReferences[neighbour]->degree-1, neighbour);

            //save for unmerging later
            added_vertices->push_back(neighbour);
        }
    }
    //std::cout << "merge_debug3" << '\n';
    //add merge vertex back to appropriate bucket
    addToBucketQueue(maxDegVertex);

    //std::cout << "merge_debug4" << '\n';

    //print();
    //printBucketQueue();

    //set vertices that will be merged inactive
    //std::cout << "set inactive: " << mergev0 << '\n';
    setInactive(mergev0);


    //std::cout << "set inactive: " << mergev1 << '\n';
    setInactive(mergev1);

    if(LP_INITIALISED)
    {
        removeFromMatching(maxDegVertex);
        queueForMatching(maxDegVertex);
    }

    std::tuple<int, std::unordered_map<int, bool>*, std::vector<int>*>* mergeInfo = new std::tuple<int, 
        std::unordered_map<int, bool>*, std::vector<int>*>(maxDegVertex, adj_map_copy, added_vertices);
    return mergeInfo;
}


void BucketGraph::unmerge(Reduction* mergeRule)
{
    //print();
    //printBucketQueue();
    /* print();
    printBucketQueue(); */
    if(mergeRule->mergeVertexInfo == nullptr)
    {
        throw std::invalid_argument("unmerge: mergeVertexInfo is nullptr");
    }
    //restore adjacency list and other fields of the merge vertex
    int mergeVertex = std::get<0>(*mergeRule->mergeVertexInfo);
    int v0 = mergeRule->deletedVertices->front();
    int v1 = mergeRule->deletedVCVertices->at(0);
    int v2 = mergeRule->deletedVCVertices->at(1);
    
    //std::cout << cp::dye("unmerging vertex " + std::to_string(mergeVertex) + " --> " + std::to_string(v0) + ", " + std::to_string(v1) + ", " + std::to_string(v2), 'g') << std::endl;

    if(mergeVertex != v0) { setActive(v0); }
    if(mergeVertex != v1) { setActive(v1); }
    if(mergeVertex != v2) { setActive(v2); }

    //std::cout << "mergeVertex: " << mergeVertex << ", N's: " << v0 << ", " << v1 << ", " << v2 << std::endl;

    //std::cout << "before remove from bucket queue" << std::endl;
    removeFromBucketQueue(mergeVertex);
    //std::cout << "after remove from bucket queue" << std::endl;
    //delete merge vertex from adjacency lists of neighbours that he didn't have before merge
    //std::cout << "unmerge_debug1" << '\n';
    std::vector<int>* added_vertices = std::get<2>(*mergeRule->mergeVertexInfo);
    for(int i = 0; i < (int) added_vertices->size(); i++)
    {
        numEdges--;
        //std::cout << "deg2 post adj pop" << std::endl;
        vertexReferences[mergeVertex]->degree--;
        //std::cout << "after pop" << std::endl;
        vertexReferences[added_vertices->at(i)]->adj_map->erase(mergeVertex);
        //vertexReferences[added_vertices->at(i)]->adj_map->erase(vertexReferences[added_vertices->at(i)]->adj_map->find(mergeVertex));
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
    //std::cout << "after added vertices edge decrement" << '\n';
    //std::cout << "mergevertex adj size:" << vertexReferences[mergeVertex]->degree << '\n';
    delete vertexReferences[mergeVertex]->adj_map;
    vertexReferences[mergeVertex]->adj_map = std::get<1>(*mergeRule->mergeVertexInfo);
    //std::cout << "mergevertex saved adj size:" << vertexReferences[mergeVertex]->degree << '\n';

    //std::cout << "after mergevertex reset" << std::endl;
    addToBucketQueue(mergeVertex);
    //std::cout << "after add to bucket queue" << std::endl;
    //std::cout << "mergeVertex: " << mergeVertex << ", N's: " << v0 << ", " << v1 << ", " << v2 << '\n';
    if(LP_INITIALISED)
    {
        removeFromMatching(mergeVertex);
        queueForMatching(mergeVertex);
    }

    // TODO: free mergeVertexInfo members
    delete std::get<2>(*mergeRule->mergeVertexInfo);
    delete mergeRule->mergeVertexInfo;
}

int BucketGraph::getReductionStackSize()
{
    if(reductions != nullptr && reductions->appliedRules != nullptr)
        return reductions->appliedRules->size();
    else
        return -1;
}

bool BucketGraph::addEdgeToVertex(int vertex, int edge)
{
    if(vertex >= (int) vertexReferences.size() || edge >= (int) vertexReferences.size() || vertex < 0 || edge < 0 || vertex == edge)
        throw std::invalid_argument("addEdgeToVertex: vertex " + std::to_string(vertex) + " or edge " + std::to_string(edge) +  " invalid\n");
    if(!isActive(vertex) || !isActive(edge))
        throw std::invalid_argument("addEdgeToVertex: vertex " + std::to_string(vertex) + " or edge " + std::to_string(edge) +  " inactive\n");
    if(vertexHasEdgeTo(vertex, edge)) { return false; }

    //add edge to vertex adj list
    vertexReferences[vertex]->adj_map->insert({edge, true});
    vertexReferences[vertex]->degree++;
    moveToBiggerBucket(vertexReferences[vertex]->degree-1, vertex);

    //add vertex to adj list of @edge
    vertexReferences[edge]->adj_map->insert({vertex, true});
    vertexReferences[edge]->degree++;
    moveToBiggerBucket(vertexReferences[edge]->degree-1, edge);

    numEdges++;
    return true;
}

void BucketGraph::removeEdgeFromVertex(int vertex, int edge)
{
    if(vertex >= (int) vertexReferences.size() || edge >= (int) vertexReferences.size() || vertex < 0 || edge < 0 || vertex == edge)
        throw std::invalid_argument("deletEdgeFromVertex: vertex " + std::to_string(vertex) + " or edge " + std::to_string(edge) +  " invalid.\n");

    if(!isActive(vertex) || !isActive(edge))
        throw std::invalid_argument("deletEdgeFromVertex: vertices " + std::to_string(vertex) + " and " + std::to_string(edge) +  " already inactive.\n");


    if(!vertexHasEdgeTo(vertex, edge))
        throw std::invalid_argument("deleteEdgeFromVertex: vertex " + std::to_string(vertex) + " doesn't have edge " + std::to_string(edge) +  " to delete.\n");

    //delete edge from vertex adj list
    //vertexReferences[vertex]->adj_map->erase(vertexReferences[vertex]->adj_map->find(edge));
    vertexReferences[vertex]->adj_map->erase(edge);
    vertexReferences[vertex]->degree--;
    moveToSmallerBucket(vertexReferences[vertex]->degree+1, vertex);

    //delete vertex from adj list of @edge
    vertexReferences[edge]->adj_map->erase(vertex);
   // vertexReferences[edge]->adj_map->erase(vertexReferences[edge]->adj_map->find(vertex));
    vertexReferences[edge]->degree--;
    moveToSmallerBucket(vertexReferences[edge]->degree+1, edge);

    numEdges--;
    return;
}




/*----------------------------------------------------------*/
/*--------------   Virtual Matching & Flow   ---------------*/
/*----------------------------------------------------------*/

bool BucketGraph::matchingBFS()
{
    std::queue<int> Q = std::queue<int>();
    if (!didInitialMatchingCalculation)
    {
        for (int i=0; i<(int) vertexReferences.size(); i++)
        {
            Vertex* active = vertexReferences[i];
            if(!isActive(active->index)) { continue; }
            if ((*pairU)[active->index] == NIL)
            {
                (*dist)[active->index] = 0;
                Q.push(active->index);
            }
            else
            {
                (*dist)[active->index] = INT32_MAX;
            }
        }
    }
    else
    {
        for (auto vertex = unmatched->begin(); vertex != unmatched->end(); ++vertex)
        {
            if(!isActive(*vertex)) { continue; }
            if ((*pairU)[*vertex] == NIL)
            {
                (*dist)[*vertex] = 0;
                Q.push(*vertex);
            }
            else
            {
                (*dist)[*vertex] = INT32_MAX;
            }
        }
    }
    (*dist)[NIL] = INT32_MAX;
    while(!Q.empty())
    {
        int u = Q.front();
        Q.pop();
        if ((*dist)[u] < (*dist)[NIL])
        {
            for (auto v = vertexReferences[u]->adj_map->begin(); v != vertexReferences[u]->adj_map->end(); ++v)
            {
                int vertex = v->first;
                if(!isActive(vertex)) { continue; }
                if ((*dist)[(*pairV)[vertex]] == INT32_MAX)
                {
                    (*dist)[(*pairV)[vertex]] = (*dist)[u] + 1;
                    Q.push((*pairV)[vertex]);
                }
            }
        }
    }
    return (*dist)[NIL] != INT32_MAX;
}

bool BucketGraph::matchingDFS(int u)
{
    if (u != NIL)
    {
        for (auto v = vertexReferences[u]->adj_map->begin(); v != vertexReferences[u]->adj_map->end(); ++v)
        {
            int vertex = v->first;
            if(!isActive(vertex)) { continue; }
            if ((*dist)[(*pairV)[vertex]] == (*dist)[u] + 1)
            {
                if (matchingDFS((*pairV)[vertex]))
                {
                    (*pairV)[vertex] = u;
                    (*pairU)[u] = vertex;
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
    int matching = 0;
    if(!didInitialMatchingCalculation) {
        while (matchingBFS())
        {
            for (int i=0; i<(int) vertexReferences.size(); i++)
            {
                Vertex* active = vertexReferences[i];
                if(!isActive(active->index)) { continue; }
                if ((*pairU)[active->index] == NIL)
                {
                    if (matchingDFS(active->index))
                    {
                        matching++;
                    }
                }
            }
        }
        currentLPBound = matching;
        didInitialMatchingCalculation = true;
    }
    else
    {
        /* std::cout << "Matching vertices {";
        for(auto vertex = unmatched->begin(); vertex != unmatched->end(); ++vertex)
        {
            std::cout << *vertex << ", ";
        }
        std::cout << "}" << '\n'; */
        while (matchingBFS())
        {
            for (auto vertex = unmatched->begin(); vertex != unmatched->end(); ++vertex)
            {
                if(!isActive(*vertex)) { continue; }
                if ((*pairU)[*vertex] == NIL)
                {
                    if (matchingDFS(*vertex))
                    {
                        matching++;
                        //std::cout << "Incrementing for new edge (" << *vertex << ", " << (*pairU)[*vertex] << ")" << '\n';
                    }
                }
            }
        }
        currentLPBound += matching;
        next_unmatched->clear();
        for(auto vertex = unmatched->begin(); vertex != unmatched->end(); ++vertex)
        {
            if(!isActive(*vertex)) { continue; }
            if ((*pairU)[*vertex] == NIL)
            {
                next_unmatched->push_back(*vertex);
            }
        }
        unmatched->clear();
        *unmatched = *next_unmatched;
    }
    //std::cout << "post hopcroft karp" << '\n';
    return currentLPBound;
}

void BucketGraph::removeFromMatching(int vertexIndex)
{
    //std::cout << "Checking edge (" << vertexIndex << ", " << (*pairU)[vertexIndex] << "), when deleting " << vertexIndex;
    if((*pairU)[vertexIndex] != NIL) {
        //std::cout << " --> Decrementing";
        //unmatched.push_back((*pairU)[vertexIndex]);
        currentLPBound--;
    }
    //std::cout << '\n';
    //std::cout << "Checking edge (" << pairV[vertexIndex] << ", " << vertexIndex << "), when deleting " << vertexIndex;
    if((*pairV)[vertexIndex] != NIL) {
        //std::cout << " --> Decrementing";
        unmatched->push_back((*pairV)[vertexIndex]);
        currentLPBound--;
    }
    //std::cout << '\n';
    (*dist)[vertexIndex] = INT32_MAX;
    if((*pairU)[vertexIndex] != NIL) { (*pairV)[(*pairU)[vertexIndex]] = NIL; }
    if((*pairV)[vertexIndex] != NIL) { (*pairU)[(*pairV)[vertexIndex]] = NIL; }
    (*pairU)[vertexIndex] = NIL;
    (*pairV)[vertexIndex] = NIL;
}

void BucketGraph::queueForMatching(int vertexIndex)
{
    unmatched->push_back(vertexIndex);
}

void BucketGraph::setBipartMatchingFlowComponentsInactive(std::vector<int>* L, std::vector<int>* R, int k, double maxExecTime)
{
    auto start = std::chrono::high_resolution_clock::now();
    //int z = 88;
    //int x = 19;
    std::stack<int> S = std::stack<int>();
    std::vector<int> visited = std::vector<int>(pairU->size()); // 0-unvisited, 1-pending for evaluation, 2-visited
    bool isNotComp = false;
    //std::cout << "Searching for sources: 0-" << pairU->size()-1 << '\n';
    // for each possible connected component (We exclude degree 0 connected components here)
    //print();
    for (int i=0; i<(int) pairU->size(); i++)
    {
        //if(visited[i] != 0) { continue; }
        //if(i==z) print();
        //if(i==z) std::cout << "Setting up for search from source: " << i << '\n';
        if (!vertexReferences[i]->isActive || (*pairU)[i] == NIL/* this culls deg=0 rule */) { continue; }
        for(int j=0; j<(int) visited.size(); j++) { visited[j] = 0; }
        while(!S.empty()) { S.pop(); }
        isNotComp = false;
        S.push(i);
        int target = (*pairU)[i];
        std::vector<int> componentL = std::vector<int>();
        std::vector<int> componentR = std::vector<int>();
        //if(i==z) std::cout << cp::dye("Started search from source: ", 'g') << i << '\n';
        // until component is closed
        while (!S.empty())
        {
            int current = S.top();
            //if(i==z) std::cout << cp::dye("Popped: ", 'p') << current << " from stack" << '\n';
            if (!vertexReferences[current]->isActive) { isNotComp = true; break; } // TODO:
            // add vertex to component, if visited
            if(visited[current] == 2 && visited[(*pairU)[current]] == 2)
            {
                //if(i==z) std::cout << cp::dye("Traversing back up from: ", 'y') << current << '\n';
                S.pop();
                continue;
            }

            // expand node
            isNotComp = true;
            visited[current] = 2;
            visited[(*pairU)[current]] = 2;
            //if(i==z) std::cout << cp::dye("Expanding: ", 'y') << current << " with (*pairU)[" << current << "] = " << (*pairU)[current] << '\n';
            for (auto v=vertexReferences[current]->adj_map->begin(); v != vertexReferences[current]->adj_map->end(); ++v)
            {
                // We skip inactive vertices and the matched right vertex and vertices of which we already visited their left matched vertex
                int vertex = v->first;
                if(!isActive(vertex)) { continue; }
                //if(i==z) std::cout << "Considering neighbour: " << *v << " and matching (" << (*pairV)[*v] << ", " << *v <<")" << '\n';
                // if we find an unmatched right vertex, abort immediately (This cannot be a component)
                if((*pairV)[vertex] == NIL || !vertexReferences[(*pairV)[vertex]]->isActive) {
                    //std::cout << "Current: " << current << " has no left matched vertex." << '\n';
                    isNotComp = true; break;
                }
                // if found target, this path is valid (we then need to check the other paths to be valid)
                if((*pairU)[i] == vertex) { isNotComp = false; continue; }
                //if(current == (*pairV)[*v]) { componentR.push_back(*v); }
                //if(visited[(*pairV)[*v]] && !visited[*v]) { componentR.push_back(*v); visited[*v] = 2; }    // TODO: debug
                // We skip the matched right vertex and vertices of which we already visited their left matched vertex
                if(current == (*pairV)[vertex] || visited[vertex] == 2) { continue; }
                if(visited[(*pairV)[vertex]] == 2 || visited[vertex] == 1 || visited[(*pairV)[vertex]] == 1) { isNotComp = true; break; }
                // push next left vertex
                S.push((*pairV)[vertex]);
                isNotComp = false;
                visited[vertex] = 1;
                visited[(*pairV)[vertex]] = 1;
                //if(i==z) std::cout << cp::dye("Pushing neighbour: ", 'p') << (*pairV)[*v] << " from matching (" << (*pairV)[*v] << ", " << *v <<")" << " to stack" << '\n';
            }
            if(isNotComp) {
                //if(i==z) std::cout << cp::dye("Current: ", 'r') << current << cp::dye(" has no successor vertex.", 'r') << '\n';
                break;
            }

            //if(i==z) std::cout << cp::dye("Adding ", 'y') << current << " and: " << (*pairU)[current] << " into the solution" << '\n';
            componentL.push_back(current);
            componentR.push_back((*pairU)[current]);
        }
        //std::cout << "Checking if component was found" << '\n';
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
            //std::cout << "}" << '\n';
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
        //std::cout << "}" << '\n';
        if(componentL.size() != componentR.size()) { throw std::invalid_argument( "component sizes not equal!!!" ); }
        // Early stopping, when k doesnt allow for reduction
        if((int) R->size() > k) { return; }
        auto current = std::chrono::high_resolution_clock::now();
        if(std::chrono::duration_cast<std::chrono::microseconds>(current - start).count() / (double) 1000000 > maxExecTime) { return; }
    }
}

void BucketGraph::scheduleComponentForUnconfined(int vertexIndex)
{
    std::queue<int> Q = std::queue<int>();
    // init Queue with neighbours of deleted vertex
    for(auto v = vertexReferences[vertexIndex]->adj_map->begin(); v != vertexReferences[vertexIndex]->adj_map->end(); ++v)
    {
        int vertex = v->first;
        if(!isActive(vertex)) { continue; }
        if(isVertexScheduledForUnconfined(vertex)) { continue; }
        Q.push(vertex);
    }
    while(!Q.empty())
    {
        int u = Q.front();
        Q.pop();
        if (!isVertexScheduledForUnconfined(u))
        {
            scheduleForUnconfined(u);
            for (auto v = vertexReferences[u]->adj_map->begin(); v != vertexReferences[u]->adj_map->end(); ++v)
            {
                int vertex = v->first;
                if(!isActive(vertex)) { continue; }
                if(isVertexScheduledForUnconfined(vertex)) { continue; }
                Q.push(vertex);
            }
        }
    }
    // TODO:
    unscheduleForUnconfined(vertexIndex);
}

/*----------------------------------------------------------*/
/*------------------   Calculate Bounds   ------------------*/
/*----------------------------------------------------------*/

int BucketGraph::getLowerBoundVC() {
    return getCliqueBound();
    //return getLPBound();
    //return getLPCycleBound();
}

int BucketGraph::getLPBound()
{
    hopcroftKarpMatchingSize();
    return currentLPBound/2;
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

    for(int i = 0; i < (int) cliques.size(); i++)
    {
        delete cliques[i];
    }
    //auto stop = std::chrono::high_resolution_clock::now();
    //auto duration = duration_cast<milliseconds>(stop - start);
    //std::cout << "numIterations: " << iterations << ", numCliqueIterations: " << cliqueIterations << ", duration: " << duration.count() << '\n';
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
        if(!vertexHasEdgeTo(vertex, (*clique)[i]))
        {
            return false;
        }
    }
    return true;
}


std::vector<std::pair<std::string,std::string>> BucketGraph::getPreprocessedEdges()
{
    using namespace std;

    vector<pair<string,string>> activeEdges;

    for (int i = 0; i < (int)edges.size(); ++i) {
        pair<string,string> edge = edges.at(i);

        // Get vertex index
        int vertex1 = originalVertexNames[edge.first].first;
        int vertex2 = originalVertexNames[edge.second].first;

        if(isActive(vertex1) && isActive(vertex2))
            activeEdges.push_back(edge);
    }

    return activeEdges;
}

void BucketGraph::unreduceDeg3Ind(Reduction *rule, std::unordered_map<int, bool> *vc)
{

    int vDeg3 = rule->deletedVertices->at(0);

    int a = rule->deletedVCVertices->at(0);
    int b = rule->deletedVCVertices->at(1);
    int c = rule->deletedVCVertices->at(2);

    setActive(vDeg3);

    if(deg3ind) {
        std::cout << "Unreduce: Degree 3 Independent Set" << '\n';
        std::cout << "Degree 3 vertex: " << vDeg3 << '\n';
        std::cout << a << ", " << b << ", " << c << '\n';
    }

    // Solution S'
    if (vc != nullptr)
    {
        if(deg3ind) {
            std::cout << "VC: ";
            if (vc->empty()) {
                std::cout << "Empty...";
            } else {
                if((int)vc->size() < 50)
                    for (const auto &pair: *vc) {
                        std::cout << pair.first << " ";
                    }
            }
            std::cout << "\nVC size: " << (int)vc->size() << '\n';
            std::cout << '\n';
        }

        int a = rule->deletedVCVertices->at(0);
        int b = rule->deletedVCVertices->at(1);
        int c = rule->deletedVCVertices->at(2);

        int commonSolution = 0;
        int inSolution[3] = {0, 0, 0};

        auto ita = vc->find(a);
        if (ita != vc->end()) {
            commonSolution++;
            inSolution[0] = 1;
        }
        auto itb = vc->find(b);
        if (itb != vc->end()) {
            commonSolution++;
            inSolution[1] = 1;
        }
        auto itc = vc->find(c);
        if (itc != vc->end()) {
            commonSolution++;
            inSolution[2] = 1;
        }

        if(deg3ind)
            std::cout << commonSolution <<" Neighbours are in S'\n";

        if(commonSolution == 1)
        {
            if(deg3ind)
            {
                std::string indcase1 = "CASE 1\n";
                std::cout << ColorPrint::dye(indcase1, 'g');
            }
            if (inSolution[0] == 1) {
                if(deg3ind)
                    std::cout << "Erasing a" << '\n';
                vc->erase(ita);
            }
            else if (inSolution[1] == 1) {
                if(deg3ind)
                    std::cout << "Erasing b" << '\n';
                vc->erase(itb);
            }
            else if (inSolution[2] == 1) {
                if(deg3ind)
                    std::cout << "Erasing c" << '\n';
                vc->erase(itc);
            }
            else
                throw std::invalid_argument("unreduce error: Ind Deg-3: unknown in case 1");
            if(deg3ind)
                std::cout << "Insert v: "<< vDeg3 << " into VC!\n" << '\n';
            vc->insert({vDeg3, true});
        }
        else if (commonSolution == 2)
        {
            if(deg3ind)
            {
                std::string indcase2 = "CASE 2\n";
                std::cout << ColorPrint::dye(indcase2, 'b');
            }
            if (inSolution[0] == 1 && inSolution[1] == 1) {
                if(deg3ind) {
                    std::cout << "a and b in Solution \n";
                    std::cout << " a: " << a << '\n';
                    std::cout << " b: " << b << '\n';
                    std::cout << "Erasing a: " << a << '\n'; }
                vc->erase(ita);
            }
            else if (inSolution[1] == 1 && inSolution[2] == 1) {
                if(deg3ind) {
                    std::cout << "b and c in Solution \n";
                    std::cout << " b: " << b << '\n';
                    std::cout << " c: " << c << '\n';
                    std::cout << "Erasing b: " << b << '\n'; }
                vc->erase(itb);
            }
            else if (inSolution[0] == 1 && inSolution[2] == 1) {
                if(deg3ind) {
                    std::cout << "c and a in Solution \n";
                    std::cout << " c: " << c << '\n';
                    std::cout << " a: " << a << '\n';
                    std::cout << "Erasing c: " << c << '\n'; }
                vc->erase(itc);
            }
            vc->insert({vDeg3, true});
            if(deg3ind) {
                std::cout << "Adding v to S: " << vDeg3 << '\n';
                std::cout << "VC size: " << (int)vc->size() << '\n';
            }
        }
        else if (commonSolution == 3)
        {
            if(deg3ind)
            {
                std::string indcase3 = "CASE 3\n";
                std::cout << ColorPrint::dye(indcase3, 'r');
            }
            if(deg3ind)
                std::cout << "All 3 Neighbours are part of the solution!" << '\n';
        }
        else{
            std::cout << "Common solution: "<< commonSolution << '\n';
            std::string errorMsg = "unreduce error: Independent Deg-3: unknown case: ";
            errorMsg += "Common solution: " + std::to_string(commonSolution) + '\n';
            throw std::invalid_argument(errorMsg);
        }
    }
        // No vc
    else {
        // Removing Edges!
        std::vector<int> addedEdgesToA = rule->addedEdges->at(0);
        std::vector<int> addedEdgesToB = rule->addedEdges->at(1);
        std::vector<int> addedEdgesToC = rule->addedEdges->at(2);

        if(deg3ind)
            std::cout << "Removing Edge:" << a <<  " with:" << '\n';
        for (int j = 0; j < (int)addedEdgesToA.size(); ++j) {
            if(deg3ind)
                std::cout << addedEdgesToA.at(j) << '\n';
            if(vertexHasEdgeTo(a, addedEdgesToA.at(j)))
                removeEdgeFromVertex(a, addedEdgesToA.at(j));
        }

        if(deg3ind)
            std::cout << '\n';

        if(deg3ind)
            std::cout << "Removing Edge:" << b <<  " with:" << '\n';
        for (int j = 0; j < (int)addedEdgesToB.size(); ++j) {
            if(deg3ind)
                std::cout << addedEdgesToB.at(j) << '\n';
            if(vertexHasEdgeTo(b, addedEdgesToB.at(j)))
                removeEdgeFromVertex(b, addedEdgesToB.at(j));
        }
        if(deg3ind)
            std::cout << '\n';

        if(deg3ind)
            std::cout << "Removing Edge:" << c <<  " with:" << '\n';
        for (int j = 0; j < (int)addedEdgesToC.size(); ++j) {
            if(deg3ind)
                std::cout << addedEdgesToC.at(j) << '\n';
            if(vertexHasEdgeTo(c, addedEdgesToC.at(j)))
                removeEdgeFromVertex(c, addedEdgesToC.at(j));
        }
        if(deg3ind)
            std::cout << '\n';
    }
}

void BucketGraph::unreduceDeg3Clique(Reduction *rule, std::unordered_map<int, bool> *vc)
{
    int vDeg3 = rule->deletedVertices->at(0);

    // 3 Case
    int c11 = rule->deletedVCVertices->at(0);
    int c12 = rule->deletedVCVertices->at(1);
    int c2 = rule->deletedVCVertices->at(2);

    if(deg3clique) {
        std::cout << "Unreduce: Deg3 Clique" << '\n';
        std::cout << "Degree 3 vertex: " << vDeg3 << '\n';
        std::cout << c11 << ", " << c12 << ", " << c2 << '\n';
    }

    setActive(vDeg3);
    setActive(c2);

    // Solution S'
    if (vc != nullptr)
    {
        if(deg3clique) {
            std::cout << "VC: ";
            if (vc->empty()) {
                std::cout << "Empty...";
            } else {
                for (const auto &pair: *vc) {
                    if(pair.first == c11 || pair.first == c12)
                        std::cout << pair.first << " ";
                }
            }
            std::cout << '\n';
        }

        int commonSolution = 0;

        auto itc11 = vc->find(c11);
        if (itc11 != vc->end()) {
            commonSolution++;
        }
        auto itc12 = vc->find(c12);
        if (itc12 != vc->end()) {
            commonSolution++;
        }

        if(deg3clique)
            std::cout << commonSolution <<" Neighbours are in S'\n";

        if(commonSolution == 1)
        {
            if(deg3clique)
                std::cout << "Inserting v to S: "<< vDeg3 << '\n';
            vc->insert({vDeg3, true});
        }
        else if (commonSolution == 2)
        {
            if(deg3clique)
                std::cout << "Adding c2 to S: " << c2 << '\n';
            vc->insert({c2, true});
        }
        else if (commonSolution == 0)
            throw std::invalid_argument("unreduce error: 2-Clique: There should be at least 1 in the common solution");
        else
            throw std::invalid_argument("unreduce error: 2-Clique: unknown case");
    }
        // No vc
    else {

        // Removing Edges!
        std::vector<int> addedEdgesToC11 = rule->addedEdges->at(0);
        std::vector<int> addedEdgesToC12 = rule->addedEdges->at(1);

        if(deg3clique)
            std::cout << "Removing Edge:" << c11 <<  " with:" << '\n';
        for (int j = 0; j < (int)addedEdgesToC11.size(); ++j) {
            if(deg3clique)
                std::cout << addedEdgesToC11.at(j) << '\n';
            if(vertexHasEdgeTo(c11, addedEdgesToC11.at(j)))
                removeEdgeFromVertex(c11, addedEdgesToC11.at(j));
        }
        if(deg3clique)
            std::cout << '\n';

        if(deg3clique)
            std::cout << "Removing Edge:" << c12 <<  " with:" << '\n';
        for (int j = 0; j < (int)addedEdgesToC12.size(); ++j) {
            if(deg3clique)
                std::cout << addedEdgesToC12.at(j) << '\n';
            if(vertexHasEdgeTo(c12, addedEdgesToC12.at(j)))
                removeEdgeFromVertex(c12, addedEdgesToC12.at(j));
        }
        if(deg3clique)
            std::cout << "\n";
    }
}

void BucketGraph::unreduceDeg3Dom(Reduction *rule, std::unordered_map<int, bool> *vc)
{
    int vDeg3 = rule->deletedVertices->at(0);

    // 3 Case
    int dom = rule->deletedVCVertices->at(0);
    int c1 = rule->deletedVCVertices->at(1);
    int c2 = rule->deletedVCVertices->at(2);

    if(deg3dom) {
        std::cout << "Unreduce: Deg3: Domination" << '\n';
        std::cout << "Degree 3 vertex: " << vDeg3 << '\n';
        std::cout << dom << ", " << c1 << ", " << c2 << '\n';
    }

    setActive(dom);

    bool clique = false;
    if (rule->kDecrement == 3) {
        setActive(c1);
        setActive(c2);
        clique = true;
    }

    // Solution S'
    if (vc != nullptr)
    {
        if(deg3dom) {
            std::cout << "VC: ";
            if (vc->empty()) {
                std::cout << "Empty...";
            } else {
                for (const auto &pair: *vc) {
                    std::cout << pair.first << " ";
                }
            }
            std::cout << '\n';
        }

        if(deg3dom) {
            if(clique){
                std::cout << "Clique!" << '\n';
                std::cout << "Inserting to VC: " << dom << ", " << c1 << ", " << c2 << '\n';
            }
            else{
                std::cout << "No Clique!" << '\n';
                std::cout << "Inserting dominator to VC: " << dom << '\n';
            }
        }

        if(clique)
        {
            vc->insert({c1, true});
            vc->insert({c2, true});
        }
        else
        {
//            if(deg3domDeg2){
//                int c1 = rule->deletedVCVertices->at(1);
//                int c2 = rule->deletedVCVertices->at(2);
//
//                int commonSolution = 0;
////                        int inSolution[3] = {0, 0, 0};
//
//                auto itc1 = vc->find(c1);
//                if (itc1 != vc->end()) {
//                    commonSolution++;
////                            inSolution[0] = 1;
//                }
//                auto itc2 = vc->find(c2);
//                if (itc2 != vc->end()) {
//                    commonSolution++;
////                            inSolution[1] = 1;
//                }
//
//                if (commonSolution == 2) {
//                    vc->insert({c1, true});
//                    vc->insert({c2, true});
//                } else if (commonSolution == 0) {
//                    vc->insert({vDeg3, true});
//                } else {
//                    std::string errorMsg = "unreduce error: Deg-3 Domination: Deg-2 Replacement: unknown case: ";
//                    errorMsg += "Common solution: " + std::to_string(commonSolution) + '\n';
//                    throw std::invalid_argument(errorMsg);
//                }
//            }
        }

        vc->insert({dom, true});

        if(deg3dom) {
            std::cout << "VC: ";
            if (vc->empty()) {
                std::cout << "Empty...";
            } else {
                for (const auto &pair: *vc) {
                    std::cout << pair.first << " ";
                }
            }
            std::cout << '\n';
            std::cout << '\n';
        }
    }
        // No vc
    else {

        if(deg3dom) {
            std::cout << "Not a Solution!" << '\n';
            if(clique){
                std::cout << "Clique!" << '\n';
                std::cout << "Setting active again: " << dom << ", " << c1 << ", " << c2 << '\n';
            }
            else{
                std::cout << "No Clique!" << '\n';
                std::cout << "Activating dominator to VC: " << dom << '\n';
            }
        }

        {
//                    if(clique){
//                        setActive(c1);
//                        setActive(c2);
//                    }
//                    else
//                    {
//                        if(deg3domDeg2)
//                        {// Removing Edges!
//                            std::vector<int> addedEdgesToC1 = rule->addedEdges->at(0);
//                            std::vector<int> addedEdgesToC2 = rule->addedEdges->at(1);
//
//                            if (deg3dom)
//                                std::cout << "Removing Edge:" << c1 << " with:" << '\n';
//                            for (int j = 0; j < (int) addedEdgesToC1.size(); ++j) {
//                                if (!isActive(c1))
//                                    break;
//                                if (deg3dom)
//                                    std::cout << addedEdgesToC1.at(j) << '\n';
//                                if (vertexHasEdgeTo(c1, addedEdgesToC1.at(j)) && isActive(addedEdgesToC1.at(j)))
//                                    removeEdgeFromVertex(c1, addedEdgesToC1.at(j));
//                            }
//                            for (int j = 0; j < (int) addedEdgesToC2.size(); ++j) {
//                                if (!isActive(c2))
//                                    break;
//                                if (deg3dom)
//                                    std::cout << addedEdgesToC2.at(j) << '\n';
//                                if (vertexHasEdgeTo(c2, addedEdgesToC2.at(j)) && isActive(addedEdgesToC2.at(j)))
//                                    removeEdgeFromVertex(c2, addedEdgesToC2.at(j));
//                            }
//                        }
//                    }
        }
    }
}

void BucketGraph::unreduceDeg4Clique(Reduction* rule, int* k, std::unordered_map<int, bool>* vc)
{
    int v = rule->deletedVertices->at(0);

    int kDecrement = rule->kDecrement;

    (*k) = (*k) + kDecrement;

    // Case 3 Connections
    int a = rule->deletedVCVertices->at(0);
    int b = rule->deletedVCVertices->at(1);
    int c = rule->deletedVCVertices->at(2);
    int d = rule->deletedVCVertices->at(3);

    if(deg4clique) {
        std::cout << "Vertex: " << v << '\n';
        std::cout << a << ", " << b << ", " << c << ", " << d << '\n';
    }


    // TODO: set the correct vertices active!
    if(kDecrement == 1)
    {
        setActive(d);
    }else
    {
        setActive(c);
        setActive(d);
    }
    setActive(v);

    if(vc != nullptr){
        if(deg4clique) {
            std::cout << "VC: ";
            if (vc->empty())
                std::cout << "empty...\n";
            else
                std::cout << "VC size = " << (int) vc->size() << '\n';
        }

        int commonSolution = 0;
        int solutions[3] = {0, 0, 0};

        auto itc1_1 = vc->find(a);
        if (itc1_1 != vc->end()) {
            commonSolution++;
            solutions[0] = 1;
        }
        auto itc1_2 = vc->find(b);
        if (itc1_2 != vc->end()) {
            commonSolution++;
            solutions[1] = 1;
        }

        if (kDecrement == 1) {

            auto itc1_3 = vc->find(c);
            if (itc1_3 != vc->end()) {
                commonSolution++;
                solutions[2] = 1;
            }

            if (commonSolution == 2) {
                vc->insert({v, true});
            } else if (commonSolution == 3) {
                vc->insert({d, true});
            } else if (commonSolution < 2)
                throw std::invalid_argument(
                        "unreduce error: deg4: 2-Clique: at least 2 Vertices in the common solution");
            else
                throw std::invalid_argument("unreduce error: deg4: 2-Clique: unknown case");
        } else {
            if (commonSolution == 1) {
                vc->insert({v, true});
                // c1_1(reduce) = a(unreduce) opposite vertex will be c(unreduce) = c2_1(reduce)
                if (solutions[0] == 1) {
                    vc->insert({c, true});
                } else if (solutions[1] == 1)
                    vc->insert({d, true});
                else
                    throw std::invalid_argument("unreduce error: deg4: 2-Clique: It should be a or b!");
            } else if (commonSolution == 2) {
                vc->insert({c, true});
                vc->insert({d, true});
            }
        }
    }
    else {
        if(deg4clique) {
            std::cout << "Not a VC, deleting edges\n";
            if(kDecrement == 1)
                std::cout << "Case 3\n";
            else
                std::cout << "Case 4\n";
        }
        // Removing Edges!
        std::vector<int> addedEdgesToC11 = rule->addedEdges->at(0);
        std::vector<int> addedEdgesToC12 = rule->addedEdges->at(1);

        if(deg4clique)
            std::cout << "Removing Edge:" << a <<  " with:" << '\n';
        for (int j = 0; j < (int)addedEdgesToC11.size(); ++j) {
            if(deg4clique)
                std::cout << addedEdgesToC11.at(j) << '\n';
            if(vertexHasEdgeTo(a, addedEdgesToC11.at(j)))
                removeEdgeFromVertex(a, addedEdgesToC11.at(j));
        }
        if(deg4clique)
            std::cout << '\n';

        if(deg4clique)
            std::cout << "Removing Edge:" << b <<  " with:" << '\n';
        for (int j = 0; j < (int)addedEdgesToC12.size(); ++j) {
            if(deg4clique)
                std::cout << addedEdgesToC12.at(j) << '\n';
            if(vertexHasEdgeTo(b, addedEdgesToC12.at(j)))
                removeEdgeFromVertex(b, addedEdgesToC12.at(j));
        }
        if(deg4clique)
            std::cout << "\n";

        if(kDecrement == 1)
        {
            std::vector<int> addedEdgesToC13 = rule->addedEdges->at(2);

            if(deg4clique)
                std::cout << "Removing Edge:" << c <<  " with:" << '\n';
            for (int j = 0; j < (int)addedEdgesToC13.size(); ++j) {
                if(deg4clique)
                    std::cout << addedEdgesToC13.at(j) << '\n';
                if(vertexHasEdgeTo(c, addedEdgesToC13.at(j)))
                    removeEdgeFromVertex(c, addedEdgesToC13.at(j));
            }
            if(deg4clique)
                std::cout << "\n";
        }
    }
}

void BucketGraph::unreduceDeg4Dom(Reduction *rule, int* k, std::unordered_map<int, bool> *vc)
{
    int v = rule->deletedVertices->at(0);

//    if(deg3dom) {
//        std::cout << "Unreduce: Deg3: Domination" << '\n';
//        std::cout << "Degree 3 vertex: " << vDeg3 << '\n';
//        std::cout << dom << ", " << c1 << ", " << c2 << '\n';
//    }
    int kDecrement = rule->kDecrement;

    (*k) = (*k) + kDecrement;

    if(kDecrement == 4)
        setActive(v);

    for (int i = 0; i < kDecrement; ++i) {
        setActive(rule->deletedVCVertices->at(i));
    }

    // Solution S'
    if (vc != nullptr)
    {
        if(deg4dom) {
            std::cout << "VC: ";
            if (vc->empty())
                std::cout << "Empty...";
            else
                std::cout << "VC size = " << (int) vc->size() << '\n';
        }

        for (int i = 0; i < kDecrement; ++i) {
            vc->insert({rule->deletedVCVertices->at(i), true});
        }

    }
        // No vc
    else {
        if(deg4dom) {
            std::cout << "Not a Solution!" << '\n';
        }
    }
}
