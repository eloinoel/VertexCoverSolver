#include <iostream> //output streams
#include <fstream>  //ifstream file opening
#include <math.h>

#include "BucketGraph.h"
#include "ColorPrint.h"

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

void BucketGraph::initBucketQueue()
{
    bucketQueue = list<Bucket>();
    bool foundDeg;
    int maxDeg = -1;
    for (auto elem : activeList)
    {
        if(elem.degree > maxDeg)
        {
            maxDeg = elem.degree;
        }

        foundDeg = false;
        for (auto bucket = bucketQueue.begin(); bucket != bucketQueue.end(); ++bucket)
        {
            if(elem.degree == bucket->degree)
            {
                bucket->insert({elem.bucketVertex});
                foundDeg = true;
            }
        }
        if(!foundDeg)
        {
            addBucket(elem.degree, {elem.bucketVertex});
        }
    }

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

void BucketGraph::removeFromBucketQueue(int degree, std::vector<BucketVertex*> vertices)
{
    bucketReferences[degree]->remove(vertices);
    if(bucketReferences[degree]->vertices.size() == 0) {
        bucketQueue.erase(bucketQueue.iterator_to(*bucketReferences[degree]));
    }
    /* for(Bucket bucket : bucketQueue)
    {
        if(bucket.degree == degree)
        {
            bucket.remove(vertices);
            if(bucket.vertices.size() == 0) {
                bucketQueue.erase(bucketQueue.iterator_to(bucket));
            }
            break;
        }
    } */
}

void BucketGraph::addToBucketQueue(int degree, std::vector<BucketVertex*> vertices)
{
    if(bucketReferences[degree]->vertices.size() == 0)
    {
        bucketQueue.insert(bucketQueue.iterator_to(*bucketReferences[degree+1]), *bucketReferences[degree]);
    }
    bucketReferences[degree]->insert(vertices);
    /* for(Bucket bucket : bucketQueue)
    {
        if(bucket.degree == degree)
        {
            bucket.insert(vertices);
            break;
        }
        else if(bucket.degree > degree)
        {
            bucketQueue.insert(bucketQueue.iterator_to(bucket), Bucket(degree, vertices));
            break;
        }
    } */
}

void BucketGraph::moveInBucketQueue(int degree, std::vector<BucketVertex*> vertices, int newDegree)
{
    /* for(Bucket bucket : bucketQueue)
    {
        if(bucket.degree == degree)
        {
            bucket.remove(vertices);
            if(bucket.vertices.size() == 0) {
                bucketQueue.erase(bucketQueue.iterator_to(bucket));
            }
            break;
        }

        if(bucket.degree == newDegree)
        {
            bucket.insert(vertices);
            break;
        }
        else if(bucket.degree > newDegree)
        {
            bucketQueue.insert(bucketQueue.iterator_to(bucket), Bucket(newDegree, vertices));
            break;
        }
    } */
}

/*----------------------------------------------------------*/
/*-------------------   Graph Utility   --------------------*/
/*----------------------------------------------------------*/

void BucketGraph::print()
{
    if (vertexReferences.size() > 0)
	{
		std::cout << "\n";
		for (int i = 0; i < (int) vertexReferences.size(); i++)
		{
			if (vertexReferences[i] != nullptr)
            {
                std::cout << "name " << dye(vertexReferences[i]->strName, 'y') << ", index " <<
                dye(std::to_string(i), 'g') << "(" << dye(std::to_string(vertexReferences[i]->degree), 'r') << ")" << ": ";

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
        moveInBucketQueue(vertexReferences[(*v->adj)[i]]->degree-1, {vertexReferences[(*v->adj)[i]]->bucketVertex}, vertexReferences[(*v->adj)[i]]->degree);
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
        moveInBucketQueue(vertexReferences[(*v->adj)[i]]->degree+1, {vertexReferences[(*v->adj)[i]]->bucketVertex}, vertexReferences[(*v->adj)[i]]->degree);
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

int BucketGraph::getMaxDegreeVertex()
{
    return bucketQueue.back().vertices[0];
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

int BucketGraph::getVerticesOfDegree(int degree)
{
    //TODO: get from bucket queue
    return 0;
}

void BucketGraph::removeFromBucketQueue(int degree, std::vector<BucketVertex*> vertices)
{
    for(auto bucket = bucketQueue.begin(); bucket != bucketQueue.end(); ++bucket)
    {
        if(bucket->degree == degree)
        {
            bucket->remove(vertices);
            if(bucket->vertices.size() == 0) {
                bucketQueue.erase(bucketQueue.iterator_to(*bucket));
            }
            break;
        }
    }
}

void BucketGraph::addToBucketQueue(int degree, std::vector<BucketVertex*> vertices)
{
    for(auto bucket = bucketQueue.begin(); bucket != bucketQueue.end(); ++bucket)
    {
        if(bucket->degree == degree)
        {
            bucket->insert(vertices);
            break;
        }
        else if(bucket->degree > degree)
        {
            bucketQueue.insert(bucketQueue.iterator_to(bucket), Bucket(degree, vertices));
            break;
        }
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

/*
* Calculates clique cover with greedy heuristic
* @param k: stop calculating cliques if bound already surpasses k
*/
int BucketGraph::getCliqueBound(int k)
{
    std::vector<std::vector<int>*> cliques = std::vector<std::vector<int>*>();
    int cliqueBound = 0;

    //iterate through vertices in ascending order of degrees (buckets)
    for(auto it = bucketQueue.begin(); it != bucketQueue.end(); ++it)
    {
        Bucket bucket = *it;
        //for each vertex in bucket
        for(auto jt = bucket.vertices.begin(); jt != bucket.vertices.end(); ++jt)
        {
            int curVertex = jt->index;
            int cliqueIndex = -1;
            for(int k = 0; k < (int) cliques.size(); k++)
            {
                if(vertexCanBeAddedToClique(curVertex, cliques.at(k)))
                {
                    cliqueIndex = k;
                    break;
                }

            }
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
    }
    return cliqueBound;
}

bool BucketGraph::vertexCanBeAddedToClique(int vertex, std::vector<int>* clique)
{
    //for each vertex in clique
    for (int i = 0; i < (int) clique->size(); i++)
    {
        //is a neighbour of vertex
        bool isNeighbour = false;
        std::vector<int>* neighbours = getNeighbours(vertex);
        for(int j = 0; j < (int) neighbours->size(); j++)
        {
            if (neighbours->at(j) == clique->at(i))
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