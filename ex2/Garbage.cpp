
//Component Subdivision Iterative Branching Algorithm Bruno, Stand: Samstag 20.05.2023 12 Uhr

/* 

// Choose max degree vertex from whole graph if there are no stacked subcomponents and
// only from component otherwise
//
int chooseVertex(ArrayGraph* G, std::stack<std::vector<std::vector<int>>>* subComponents)
{	// TODO: this doesnt work
	if(subComponents->empty())
	{
		return G->getMaxDegreeVertex();
	}
	else
	{
		// TODO: could always choose largest or smallest component, as desired, as a heuristic
		return G->getMaxDegreeVertex(&(subComponents->top().back()));
	}
}

vector<int>* VCVertexBranchingIterative(ArrayGraph* G, int k, std::vector<int>* vc, bool useDegLEQ2Alg, VCDebugMode debug)
{
	// stack storing the differentials of partial solutions, currently under evaluation and whether a partial solution was already expanded
    std::stack<std::pair<std::vector<int>, bool>> S = std::stack<std::pair<std::vector<int>, bool>>();
	std::stack<std::vector<std::vector<int>>> subComponents = std::stack<std::vector<std::vector<int>>>();
	std::pair<std::vector<int>, bool> current;
	int branchVertex;
	int branchVertexDegree;
	// number of currently active vertices
	int partialVCSize = 0;
	
	// initialize stack
	std::vector<int> bv = {};
	std::pair<std::vector<int>, bool> childDifferential = {bv, false};
	S.push(childDifferential);

	std::vector<std::vector<int>> init = std::vector<std::vector<int>>();
	init.push_back(std::vector<int>());
	for(int i=0; i<G->getNumberOfVertices(); i++)
	{
		init.at(0).push_back(i);
	}
	std::vector<int> origins = std::vector<int>();
	for(int i=0; i<G->getNumberOfVertices(); i++)
	{
		origins.push_back(i);
	}
	std::cout << "Init: ";
	for (auto i: init) {
		for (auto j: i)
			std::cout << j << ' ';
		std::cout << '\n';
	}

	subComponents.push(*(G->getComponents(&origins)));
	if(debug == ExecutionTranscript)
	{
		std::cout << dye("started", 'g') << " BnB with stack size: " << S.size() << " and k=" << k << "\n";
	}

	int it = 0;
	while (!S.empty())
	{
		if(debug == IterationsTaken) it++;
		// retrieve the differential of the current partial vertex cover solution to its parent solution (+ expanded tag)
		current = S.top();

		if (current.second)
		{
			if(debug == ExecutionTranscript)
			{
				std::cout << tileStr("--", partialVCSize) << "- " << dye("Traversing", 'y') << " back up a previously expanded search tree node\n";
				std::cout << tileStr("--", partialVCSize) << "- " << dye("Restoring", 'p') << " vertices: {";
				if (current.first.size() > 0) std::cout << current.first.at(0);
				for (int i=1; i < (int) current.first.size(); i++)
				{
					std::cout << ", " << current.first.at(i);
				}
				std::cout << "}\n";
			}
			G->setActive(&current.first);
			partialVCSize -= current.first.size();
			subComponents.pop();
			S.pop();
			continue;
		}

		// otherwise delete the vertices (defined by the differential calculated in the expansion of the parent search tree node)
		G->setInactive(&current.first);
		partialVCSize += current.first.size();

		// calculate subcomponents and add them to the component partition stack
		auto nextComponents = subComponents.top();
		if(!current.first.empty()) nextComponents.pop_back();
		for(auto component : *(G->getComponents(G->getNeighbours(&current.first))))
		{
			nextComponents.push_back(component);
		}
		subComponents.push(nextComponents);

		std::cout << "next: ";
		for (auto i: nextComponents) {
			for (auto j: i)
				std::cout << j << ' ';
			std::cout << '\n';
		}
		// unify (top() / component with current.first) with *(G->getComponents(&current.first))

		if(debug == ExecutionTranscript)
		{
			auto SP = S;
			std::cout << tileStr("--", partialVCSize) << "- " << dye("peeking", 'y') << " stack: {";
			while (!SP.empty())
			{
				auto ccurrent = SP.top();
				SP.pop();
				std::cout << "{";
				if (ccurrent.first.size() > 0) std::cout << ccurrent.first.at(0);
				for (int i=1; i< (int) ccurrent.first.size(); i++)
				{
					std::cout << ", " << ccurrent.first.at(i);
				}
				std::cout << "}";
				if(!SP.empty()) std::cout << ", ";
			}
			std::cout << "} of size: " << S.size() << "\n";

			std::cout << tileStr("--", partialVCSize) << "- " << dye("Deleting", 'p') << " vertices: {";
			if (current.first.size() > 0) std::cout << current.first.at(0);
			for (int i=1; i< (int) current.first.size(); i++)
			{
				std::cout << ", " << current.first.at(i);
			}
			std::cout << "}\n";
			std::cout << tileStr("--", partialVCSize) << "- " << dye("checking", 'y') << " partial solution:\n";
			std::cout << tileStr("--", partialVCSize) << "- " << "k=" << k << ", partialVCSize: " << partialVCSize << "\n";
		}

		if(debug == StackAccessTranscript)
		{
			std::cout << tileStr("-", partialVCSize) << dye("{", 'r');
			if (current.first.size() > 0) std::cout << current.first.at(0);
			for (int i=1; i< (int) current.first.size(); i++)
			{
				std::cout << dye(", ", 'r') << current.first.at(i);
			}
			std::cout << dye("}", 'r') << "\n";
		}

		if (subComponents.top().size() <= 0)
		{
			// if current solution is actually better than the current best solution: update k & vc
			vc = G->getInactiveVertices();
			k = partialVCSize;

			if(debug == ExecutionTranscript)
			{
				std::cout << tileStr("--", partialVCSize) << "> " << dye("found", 'g') << " VC: {";
				if (vc->size() > 0) std::cout << vc->at(0);
				for (int i=1; i < (int) vc->size(); i++)
				{
					std::cout << ", " << vc->at(i);
				}
				std::cout << "} of size: " << partialVCSize << "\n";
			}
			if(debug == ExecutionTranscript)
			{
				std::cout << tileStr("--", partialVCSize) << "- " << dye("Restoring", 'p') << " vertices: {";
				if (current.first.size() > 0) std::cout << current.first.at(0);
				for (int i=1; i < (int) current.first.size(); i++)
				{
					std::cout << ", " << current.first.at(i);
				}
				std::cout << "}\n";
			}
			
			if(debug == IterationsTaken) std::cout << "(k=" << k << ") solved in " << it << " iterations.\n";
			// return vertex cover
			return vc;
		}

		// get maxDegVert of remaining active vertices (the way the algorithm branches is determined by the choice of this vertex)
		branchVertex = G->getMaxDegreeVertex(&(subComponents.top().back())); // chooseVertex(G, &subComponents);
		branchVertexDegree = G->getVertexDegree(branchVertex);

		if(debug == ExecutionTranscript)
		{
			std::cout << tileStr("--", partialVCSize) << "> " << dye("selected", 'b') << " branchVertex: " << branchVertex << " with degree " << branchVertexDegree << "\n";
		}

		// if no active vertex left in graph or no vertex with degree >= 0: (We found a solution)
		if (branchVertex == -1 || branchVertexDegree == 0)
		{
			//if(subComponents.top().empty()) {
				// if current solution is actually better than the current best solution: update k & vc
				vc = G->getInactiveVertices();
				k = partialVCSize;

				if(debug == ExecutionTranscript)
				{
					std::cout << tileStr("--", partialVCSize) << "> " << dye("found", 'g') << " VC: {";
					if (vc->size() > 0) std::cout << vc->at(0);
					for (int i=1; i < (int) vc->size(); i++)
					{
						std::cout << ", " << vc->at(i);
					}
					std::cout << "} of size: " << partialVCSize << "\n";
				}
				if(debug == ExecutionTranscript)
				{
					std::cout << tileStr("--", partialVCSize) << "- " << dye("Restoring", 'p') << " vertices: {";
					if (current.first.size() > 0) std::cout << current.first.at(0);
					for (int i=1; i < (int) current.first.size(); i++)
					{
						std::cout << ", " << current.first.at(i);
					}
					std::cout << "}\n";
				}
				
				if(debug == IterationsTaken) std::cout << "(k=" << k << ") solved in " << it << " iterations.\n";
				// return vertex cover
				return vc;
			}
			else
			{
				//TODO: traverse up
			}
		}

		// if maximum search depth k is reached
		// or if the current solution has already been expanded, revert its deletion of vertices and traverse back up to the next partial solution
		if (k <= partialVCSize)
		{
			if(debug == ExecutionTranscript)
			{
				std::cout << tileStr("--", partialVCSize) << "> " << dye("reached", 'c') << " search tree depth k=" << k << "\n";
				std::cout << tileStr("--", partialVCSize) << "- " << dye("Restoring", 'p') << " vertices: {";
				if (current.first.size() > 0) std::cout << current.first.at(0);
				for (int i=1; i < (int) current.first.size(); i++)
				{
					std::cout << ", " << current.first.at(i);
				}
				std::cout << "}\n";
			}
			if(debug == StackAccessTranscript) {
				std::cout << tileStr("-", partialVCSize) << dye("{", 'g');
				if (current.first.size() > 0) std::cout << current.first.at(0);
				for (int i=1; i< (int) current.first.size(); i++)
				{
					std::cout << dye(", ", 'g') << current.first.at(i);
				}
				std::cout << dye("}", 'g') << "\n";
			}

			G->setActive(&current.first);
			partialVCSize -= current.first.size();
			subComponents.pop();
			S.pop();
			continue;
		}

		// mark this partial solution as expanded and expand subsequently
		current.second = true;
		S.pop();
		S.push(current);

		// solve graph with maxVertDegree <= 2 in linear time // TODO: Linear time alg doesnt work with subComponents yet (needs to only calc it for currennt subcomponent)
		if (useDegLEQ2Alg && branchVertexDegree <= 2)
		{
			// determine partial 2 VC for each connected component
			std::vector<int> deleted;
			std::vector<int> A;
			std::vector<int> CoA;
			bool maxDepthReached = false;
			while(true)
			{

				int origin = G->getConnectedVertex();
				// if vertex cover is found, return it
				if(origin == -1) {
					vc = G->getInactiveVertices();
					k = partialVCSize;
					return vc;
				}
				auto neighbours = G->getNeighbours(origin);
				int current1 = -1;
				int current2 = -1;

				if(G->getVertexDegree(origin) < 2)
				{
					current1 = neighbours->at(0);
				}
				else
				{
					current1 = neighbours->at(0);
					current2 = neighbours->at(1);
				}
				int lastCurrent1 = origin;
				int lastCurrent2 = origin;
				bool addToA = false;
				A.push_back(origin);

				while(current1 != -1 || current2 != -1)
				{
					// check if maximum search depth k is reached
					if((current1 != -1 && current2 != -1 && k - partialVCSize - std::min(A.size(), CoA.size()) < 2) || k - partialVCSize - std::min(A.size(), CoA.size()) < 1)
					{
						G->setActive(&deleted);
						partialVCSize -= deleted.size();
						G->setActive(&(current.first));
						partialVCSize -= current.first.size();
						S.pop();
						maxDepthReached = true;
						break;
					}

					if(current1 == current2)
					{
						if(addToA)
						{
							A.push_back(current1);
						}
						else
						{
							CoA.push_back(current1);
						}
						break;
					}

					if(current1 != -1)
					{
						if(addToA)
						{
							A.push_back(current1);
						}
						else
						{
							CoA.push_back(current1);
						}
						
						int nextCurrent1 = -1;
						for(auto neighbour : *(G->getNeighbours(current1)))
						{
							if(lastCurrent1 != neighbour) { nextCurrent1 = neighbour; }
						}
						lastCurrent1 = current1;
						current1 = nextCurrent1;
					}

					if(current2 != -1)
					{
						if(addToA)
						{
							A.push_back(current2);
						}
						else
						{
							CoA.push_back(current2);
						}
						
						int nextCurrent2 = -1;
						for(auto neighbour : *(G->getNeighbours(current2)))
						{
							if(lastCurrent2 != neighbour) { nextCurrent2 = neighbour; }
						}
						lastCurrent2 = current2;
						current2 = nextCurrent2;
					}
					addToA = !addToA;
				}
				if(maxDepthReached) break;

				// choose smaller partial vc, delete the selected vertices from the graph and save their indices in deleted
				if(A.size() <= CoA.size())
				{
					deleted.insert(
						deleted.end(),
						std::make_move_iterator(A.begin()),
						std::make_move_iterator(A.end())
					);
					G->setInactive(&A);
					partialVCSize += A.size();
				}
				else
				{
					deleted.insert(
						deleted.end(),
						std::make_move_iterator(CoA.begin()),
						std::make_move_iterator(CoA.end())
					);
					G->setInactive(&CoA);
					partialVCSize += CoA.size();
				}
				A.clear();
				CoA.clear();
			}
		}

		// refined search tree branching for maxVertDegree >= 3
		else if ((!useDegLEQ2Alg && branchVertexDegree >= 1) || (useDegLEQ2Alg && branchVertexDegree >= 3))
		{
			// if k and current partial VC size permit adding the neighbours
			if (k - partialVCSize >= branchVertexDegree)
			{
				if(debug == ExecutionTranscript)
				{
					std::cout << tileStr("--", partialVCSize) << "> " << dye("pushing", 'p') << " maxDegreeVertex neighbourhood deletion: {";
					if (G->getNeighbours(branchVertex)->size() > 0) std::cout << G->getNeighbours(branchVertex)->at(0);
					for (int i=1; i<(int) G->getNeighbours(branchVertex)->size(); i++)
					{
						std::cout << ", " << G->getNeighbours(branchVertex)->at(i);
					}
					std::cout << "}\n";
				}

				// add neighbours of the current vertex to the child differential for evaluating the partial vertex cover where all of the branchVertex's neighbours where taken into the vertex cover
				// TODO: is there no way to create this inline?
				std::pair<std::vector<int>, bool> childDifferentialN({*(G->getNeighbours(branchVertex)), false});
				S.push(childDifferentialN);
			}
			// if k and current partial VC size permit adding the current vertex
			if (k - partialVCSize >= 1)
			{
				if(debug == ExecutionTranscript)
				{
					std::cout << tileStr("--", partialVCSize) << "> " << dye("pushing", 'p') << " maxDegreeVertex deletion: {" << branchVertex << "}\n";
				}

				// TODO: is there no way to create this inline?
				std::vector<int> bv = std::vector<int>();
				bv.push_back(branchVertex);
				std::pair<std::vector<int>, bool> childDifferentialV({bv, false});
				S.push(childDifferentialV);
			}

		}
	}

	if(debug == IterationsTaken && vc == nullptr) { std::cout << "(k=" << k << ") determined to have no solution in " << it << " iterations.\n"; }
	return vc;
} */