from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from builtins import open
from builtins import range
from builtins import str
from builtins import int
from future import standard_library
standard_library.install_aliases()
# coding: utf-8

# In[ ]:

from networkit import *
from betterRepair import repairPartition
import math, sys, subprocess, functools, operator, time, random, os


# In[ ]:

def getCutWeight(G, part, v, block):
    n = G.numberOfNodes()
    z = G.upperNodeIdBound()
    assert(G.hasNode(v))
    assert(len(part) == n)
    assert(block in part)
    
    return sum([G.weight(v, u) for u in G.nodes() if G.hasEdge(v,u) and part[u] == block])


# In[ ]:

def dpPartition(G, k, imbalance, isCharged=[], useLowerBounds=False):
    """
    Partition G into subsets of size at most math.ceil(n/k)*(1+imbalance) and with consecutive node ids.
    Charged nodes are not grouped into the same subset.
    """

    # validate input
    n = G.numberOfNodes()
    if len(isCharged) > 0:
        assert(len(isCharged)==n)
        assert(sum(isCharged) <= k)
    else:
        isCharged = [False for i in range(n)]
    assert(k > 1)
    assert(k <= n)
    assert(imbalance >= 0)
    maxBlockSize = int(math.ceil(n / k)*(1+imbalance))
    minBlockSize = max(math.ceil(math.floor(n / k)*(1-imbalance)), 1) if useLowerBounds else 1
    
    # allocate cut and predecessor table
    table = [[float("inf") for j in range(k)] for i in range(n)]
    pred = [[-1 for j in range(k)] for i in range(n)]
    
    # fill values for the first fragment
    chargeEncountered = False
    weightSum = 0
    for i in range(min(maxBlockSize, n)):    
        # a fragment may only contain one charge, stop when a second one is encountered
        if isCharged[i]:
            if chargeEncountered:
                break
            else:
                chargeEncountered = True
                
        # update current weight sum
        for neighbor in G.neighbors(i):
            if neighbor > i:
                weightSum += G.weight(i,neighbor)
            elif neighbor < i:
                weightSum -= G.weight(i,neighbor)
                
        table[i][0] = weightSum if i >= minBlockSize -1 else float("inf")
        
    # fill remaining values 
    for i in range(n):      
        windowStart = max(i-maxBlockSize,0)
              
        # make sure that no two charged nodes are in the same partition
        chargeEncountered = False
        for l in reversed(range(windowStart, i+1)):
            assert(l >= windowStart)
            if isCharged[l]:
                if chargeEncountered:
                    windowStart = l
                    break
                else:
                    chargeEncountered = True
                    
        # fill cost array
        costArray = []
        cutWeight = 0
        for l in reversed(range(windowStart, i+1)):
            for neighbor in G.neighbors(l):
                # we count only edges to nodes with higher ids, to avoid double counting
                if neighbor > i:
                    cutWeight += G.weight(l,neighbor)
            costArray.append(cutWeight)
        
        # calculate optimal next fragment
        for j in range(1,k):
            predList = [table[l][j-1] + costArray[i-l-1] for l in range(windowStart, i-minBlockSize+1)]
            if (len(predList) > 0):
                minPred = min(predList)
                table[i][j] = minPred
                pred[i][j] = predList.index(minPred) + windowStart
                                      
    # get result from table
    bestCutValue = table[n-1][k-1]
        
    if (bestCutValue == float("inf")):
        raise ValueError("Combination n="+str(n)+", k="+str(k)+", epsilon="+str(imbalance)+" and chargedNodes="+str([i for i in G.nodes() if isCharged[i]])+" allows no partition!")
    result = partitioning.Partition(n)
    result.setUpperBound(k)
    
    # search best path backwards
    j = k-1
    i = n-1
    
    while (j > 0):
        nextI = pred[i][j]
        assert(nextI >= 0)
        # assign partitions to nodes
        for l in range(nextI+1, i+1):
            result[l] = j
        j -= 1
        i = nextI
        
    # assign partitions to first nodes not covered by previous loop
    for l in range(0, nextI+1):
        result[l] = 0
        
    # check results:
    for i in range(n):
        assert(result[i] >= 0)
        assert(result[i] < k)
        
    #if table[n-1][k-1] != partitioning.computeEdgeCut(result, G):
    #    print(table[n-1][k-1], 'vs', partitioning.computeEdgeCut(result, G))
        
    for size in result.subsetSizes():
        if (size > maxBlockSize):
            print("For n=", n, ", k=", k, "imbalance=", maxImbalance , ", ", size, " is wrong.")
        assert(size <= maxBlockSize)
    
    return result


# In[ ]:

def naivePartition(G, k):
    """
    Chop a new fragment off G every n/k nodes
    """
    n = G.numberOfNodes()
    naivePart = partitioning.Partition(n)
    naivePart.allToSingletons()
    for i in range(n):
        naivePart.moveToSubset(int(i/math.ceil(n/k)), i)
    naivePart.compact()
    return naivePart


# In[ ]:

def greedyPartition(G, k, imbalance, isCharged=[]):
    """
    Starting with singleton clusters, greedily merge the heaviest edge as long as it is smaller than sizelimit.
    """
    n = G.numberOfNodes()
    if len(isCharged) > 0:
        assert(len(isCharged)==n)
    else:
        isCharged = [False for i in range(n)]
    n = G.numberOfNodes()
    part = partitioning.Partition(n)
    part.allToSingletons()
    chargedPartitions = set([part.subsetOf(i) for i in range(n) if isCharged[i]])
    sizelimit = int(math.ceil(n / k)*(1+imbalance))
    remainingFragments = n

    
    def getWeight(edge):
        return G.weight(edge[0], edge[1])
    
    sortedEdges = sorted(G.edges(), key=getWeight)
    
    # merge heaviest edge, as long as allowed
    while len(sortedEdges) > 0 and remainingFragments > k:
        allowed = True
        heaviestEdge = sortedEdges.pop()
        firstPart = part.subsetOf(heaviestEdge[0])
        secondPart = part.subsetOf(heaviestEdge[1])
        if firstPart in chargedPartitions and secondPart in chargedPartitions:
            allowed = False
        sizeMap = part.subsetSizeMap()
        if sizeMap[firstPart] + sizeMap[secondPart] > sizelimit:
            allowed = False
        partSet = {firstPart, secondPart}
        for i in range(n-2):
            if part[i] in partSet and part[i+2] in partSet and not part[i+1] in partSet:
                allowed = False #otherwise, would create single embedded node
        if allowed:
            part.mergeSubsets(firstPart, secondPart)
            remainingFragments -= 1
            if firstPart in chargedPartitions or secondPart in chargedPartitions:
                chargedPartitions.add(part.subsetOf(heaviestEdge[0]))
    
    part.compact()
    return part


# In[ ]:

def mlPartition(G, k, imbalance, isCharged=[], bisectRecursively = False, avoidGaps = False):
    """
    Use a multi-level approach with Fiduccia-Matheyses to partition G.

    Subsets have size at most (1+imbalance)*ceil(n/k)
    """
    n = G.numberOfNodes()
    if len(isCharged) > 0:
        assert(len(isCharged)==n)
        if k > 0:
            assert(sum(isCharged) <= k)
    else:
        isCharged = [False for i in range(n)]
    
    listOfChargedNodes = [i for i in range(n) if isCharged[i]]
    greedy = greedyPartition(G, k, imbalance, isCharged)
    try:
        dynamic = dpPartition(G, k, imbalance, isCharged)
        if partitioning.computeEdgeCut(greedy, G) < partitioning.computeEdgeCut(dynamic, G):
            initial = greedy
        else:
            initial = dynamic
    except ValueError:
        initial = greedy
    # Problem: The single node repair in C++ may create invalid partitions.
    # Better to not use it and repair in Python.
    mlp = partitioning.MultiLevelPartitioner(G, k, imbalance, bisectRecursively, listOfChargedNodes, avoidGaps, initial)
    mlp.run()
    return repairPartition(G, mlp.getPartition(), imbalance, isCharged)


# In[ ]:

def kaHiPWrapper(G, k, imbalance = 0.2, pathToKaHiP = '/home/moritzl/Gadgets/KaHIP/deploy/kaffpa', multiple=False):
    """
    Calls KaHiP, an external partitioner.
    """
    
    tempFileName = 'tempForKaHiP.graph'
    outputFileName = 'tmppartition'+str(k)
    n = G.numberOfNodes()
    
    maxWeight = max([G.weight(u,v) for (u,v) in G.edges()])
    
    """
    KaHiP only accepts integer weights, thus we scale and round them.
    Weights must be under 1 million, otherwise the METIS graph writer switches to scientific notation,
    which confuses KaHiP
    """
    scalingFactor = int((10**6-1)/maxWeight)
    
    
    #copy and scale graph
    Gscaled = G.copyNodes()
    for (u,v) in G.edges():
        Gscaled.addEdge(u,v,int(G.weight(u,v)*scalingFactor))
    
    # write out temporary file
    writeGraph(Gscaled, tempFileName, Format.METIS)
    
    # call KaHIP
    callList = [pathToKaHiP, '--k='+str(k), '--imbalance='+str(int(imbalance*100)), '--preconfiguration=strong']
    if multiple:
        callList.append('--time_limit=1')
    callList.append(tempFileName)
    subprocess.call(callList)
    
    # read in partition
    part = community.PartitionReader().read(outputFileName)
    
    # remove temporary files
    subprocess.call(['rm', tempFileName])
    subprocess.call(['rm', outputFileName])
    
    return part


# In[ ]:

def getBestCut(G, k, imbalance, isCharged = []):
    """
    Executes the multilevel, greedy and dynamic programming algorithm, also calls KaHiP if available.
    Returns the result yielding the best cut weight.
    """
    n = G.numberOfNodes()
    if len(isCharged) == 0:
        isCharged = [False for v in range(G.numberOfNodes())]
    sizelimit = int(math.ceil(n / k)*(1+imbalance))
    
    ml = mlPartition(G, k, imbalance, isCharged)
    if not partitionValid(G, ml, sizelimit, isCharged):
        ml = repairPartition(G, ml, imbalance, isCharged)
    result = ml
    resultWeight = partitioning.computeEdgeCut(result, G)
        
    greedy = greedyPartition(G, k, imbalance, isCharged)
    if not partitionValid(G, greedy, sizelimit, isCharged):
        greedy = repairPartition(G, greedy, imbalance, isCharged)
        assert(partitionValid(G, greedy, sizelimit, isCharged))
    cutWeight = partitioning.computeEdgeCut(greedy, G)
    if cutWeight < resultWeight:
        result = greedy
        resultWeight = cutWeight
    
    try:
        cont = dpPartition(G, k, imbalance, isCharged)
        assert(partitionValid(G, cont, sizelimit, isCharged))
        cutWeight = partitioning.computeEdgeCut(cont, G)
        if cutWeight < resultWeight:
            result = cont
            resultWeight = cutWeight
    except ValueError as e:
        print(e)
        print("Continuing with other partitioners.")

    try:
        ka = kaHiPWrapper(G, k, imbalance)
        if not partitionValid(G, ka, sizelimit, isCharged):
            ka = repairPartition(G, ka, imbalance, isCharged)
        cutWeight = partitioning.computeEdgeCut(ka, G)
        if cutWeight < resultWeight:
            result = ka
            resultWeight = cutWeight
    except FileNotFoundError as e:
        print("Could not find KaHiP:",e)
        print("Continuing with other partitioners.")
        
    naive = naivePartition(G, k)
    if not partitionValid(G, naive, sizelimit, isCharged):
        naive = repairPartition(G, naive, imbalance, isCharged)        
        assert(partitionValid(G, naive, sizelimit, isCharged))
    cutWeight = partitioning.computeEdgeCut(naive, G)
    if cutWeight < resultWeight:
        result = naive
        resultWeight = cutWeight
        
    return result


# In[ ]:

def spiralLayout(G, k, rowheight = 10, colwidth = 10):
    """
    Return two lists, of x and y coordinates for a spiral layout of G.

    k nodes are put in one row, keywords rowheight and colwidth determine spacing
    """
    n = G.numberOfNodes()
    z = G.upperNodeIdBound()
    x = [0 for i in range(z)]
    y = [0 for i in range(z)]
    for i in range(z):
        if G.hasNode(i):
            if int(i / k) % 2 > 0:
                x[i] = colwidth*(k-(i % k)-1)
            else:
                x[i] = colwidth*(i % k)
            
            y[i] = rowheight*int(i / k)
            
            # adapt coordinates for rounded bends
            
            ydelta = int(rowheight / 4)
            xdelta = colwidth*(1-math.cos(math.pi/3))
            rightwards = int(i / k) % 2 == 0
    
            if i % k == k-1:
                y[i] += ydelta
                x[i] = x[i] - xdelta if rightwards else x[i] + xdelta
            if i > 0 and i % k == 0:
                y[i] -= ydelta
                x[i] = x[i] - xdelta if not rightwards else x[i] + xdelta
        
    for i in range(z):
        x[i] += 1# gephi ignores coordinates with value 0
        y[i] += 1
    return x, y


# In[ ]:

def exportToGephi(G, xcoords, ycoords, part):
    """
    Export graph to Gephi, along with coordinates and partition
    """
    client = gephi.streaming.GephiStreamingClient()
    client.clearGraph()
    client.exportGraph(G)
    client.exportNodeValues(G, part, "partition")
    client.exportNodeValues(G, xcoords, 'x')
    client.exportNodeValues(G, [-elem for elem in ycoords], 'y')
    client.exportEdgeValues(G, [G.weight(u,v) for u,v in G.edges()], 'Weight')

# In[ ]:

def partitionValid(G, partition, maxBlockSize = 0, isCharged = []):
    z = G.upperNodeIdBound()
    n = G.numberOfNodes()
    if len(partition) != z:
        return False
    
    if len(isCharged) != 0 and len(isCharged) != z:
        return False
    
    if len(isCharged) == 0:
        isCharged = [False for i in range(z)]
        
    if maxBlockSize == 0:
        maxBlockSize = n
        
    chargedFragments = set()
    
    fragmentSizes = {}
    
    for v in range(G.numberOfNodes()):
        if not G.hasNode(v):
            print("Node ", v, " not in graph.")
            return False
            
        # partition invalid if two charged nodes in same fragment
        if isCharged[v]:
            if partition[v] in chargedFragments:
                print("Node", v, " is charged, but fragment", partition[v], "already has a charged node.")
                return False
            else:
                chargedFragments.add(partition[v])
    
        # partition also invalid if gaps of size 1 exist
        if G.hasNode(v+2) and partition[v+2] == partition[v] and G.hasNode(v+1) and partition[v+1] != partition[v]:
            print("Nodes", v, "and", v+2, "are in fragment", partition[v], "but", v+1, "is in fragment", partition[v+1])
            return False
    
        # partition invalid if fragment is larger than allowed
        if not partition[v] in fragmentSizes:
            fragmentSizes[partition[v]] = 1
        else:
            fragmentSizes[partition[v]] += 1
        if fragmentSizes[partition[v]] > maxBlockSize:
            print("Fragment", partition[v], "contains", fragmentSizes[partition[v]], "nodes, more than", maxBlockSize)
            return False
    
    # no reason to complain found, partition is valid
    return True


# In[ ]:

def chargesValid(G, klist, minEpsilon, isCharged):
    assert(type(klist) is list)
    assert(type(minEpsilon) is float)
    assert(type(isCharged) is list)
    assert(len(isCharged) == G.numberOfNodes())
    for k in klist:
        try:
            part = dpPartition(G, k, minEpsilon, isCharged)
        except ValueError as e:
            return False
    return True


# In[ ]:

def comparePartitionQuality(G, k, imbalance, chargedNodes = set(), silent=False):
    n = G.numberOfNodes()
    
    isCharged = [v in chargedNodes for v in range(G.numberOfNodes())]
    sizelimit = int(math.ceil(n / k)*(1+imbalance))
    if not silent:
        print("Size limit:", sizelimit)
    result = {}
    
    before = time.time()
    ml = mlPartition(G, k, imbalance, isCharged)
    timeML = time.time() - before
    if not silent:
        print("MultiLevel:", partitioning.computeEdgeCut(ml, G))
        print("Time:", timeML)
    if not partitionValid(G, ml, sizelimit, isCharged):
        ml = repairPartition(G, ml, imbalance, isCharged)
        if not silent:
            print("Repaired Multilevel:", partitioning.computeEdgeCut(ml, G))
            partitionValid(G, ml, sizelimit, isCharged)
    if not silent:
        print("Effective k", str(ml.numberOfSubsets()))
        print()
    result['ml'] = partitioning.computeEdgeCut(ml, G)
    
    before = time.time()
    greedy = greedyPartition(G, k, imbalance, isCharged)
    timeGreedy = time.time() - before
    if not silent:
        print("Greedy:", partitioning.computeEdgeCut(greedy, G))
        print("Time:", timeGreedy)
    if not partitionValid(G, greedy, sizelimit, isCharged):
        greedy = repairPartition(G, greedy, imbalance, isCharged)
        if not silent:
            print("Repaired Greedy:", partitioning.computeEdgeCut(greedy, G))
        assert(partitionValid(G, greedy, sizelimit, isCharged))
    if not silent:
        print("Effective k", str(greedy.numberOfSubsets()))
        print()
    result['greedy'] = partitioning.computeEdgeCut(greedy, G)
    
    X = int(n / k)
    tolerance = int(math.ceil(n / k)*(1+imbalance)) - X
    
    try:
        before = time.time()
        cont = dpPartition(G, k, imbalance, isCharged)
        timeDP = time.time() - before
        if not silent:
            print("Dynamic Programming:", partitioning.computeEdgeCut(cont, G))
            print("Time:", timeDP)
        if not partitionValid(G, cont, sizelimit, isCharged):
            cont = repairPartition(G, cont, imbalance, isCharged)
            if not silent:
                print("Repaired Dynamic:", partitioning.computeEdgeCut(cont, G))
            assert(partitionValid(G, cont, sizelimit, isCharged))
        result['cont'] = partitioning.computeEdgeCut(cont, G)
        if not silent:
            print("Effective k", str(cont.numberOfSubsets()))
    except ValueError as e:
        print(e)

    print()
        
    before = time.time()
    ka = kaHiPWrapper(G, k, imbalance)
    timeKa = time.time() - before
    if not silent:
        print("Raw KaHIP:", partitioning.computeEdgeCut(ka, G))
        print("Time:", timeKa)
    if not partitionValid(G, ka, sizelimit, isCharged):
        ka = repairPartition(G, ka, imbalance, isCharged)
        if not silent:
            print("Repaired KaHiP:", partitioning.computeEdgeCut(ka, G))
            partitionValid(G, ka, sizelimit, isCharged)
    if not silent:
        print("Effective k", str(ka.numberOfSubsets()))
        print()
    result['ka'] = partitioning.computeEdgeCut(ka, G)
    result['bestOfFour'] = min([result[key] for key in result])


    before = time.time()
    naive = naivePartition(G, k)
    timeNaive = time.time() - before
    if not silent:
        print("Naive:", partitioning.computeEdgeCut(naive, G))
        print("Time:", timeNaive)
    if not partitionValid(G, naive, sizelimit, isCharged):
        naive = repairPartition(G, naive, imbalance, isCharged)        
        if not silent:
            print("Repaired Naive:", partitioning.computeEdgeCut(naive, G))
        assert(partitionValid(G, naive, sizelimit, isCharged))
    if not silent:
        print("Effective k", str(naive.numberOfSubsets()))
    result['naive'] = partitioning.computeEdgeCut(naive, G)
    result['bestOfFive'] = min([result[key] for key in result])

    if not silent:
        print(str(result['bestOfFour'] / result['naive']))
    return result


# In[ ]:

def readCharges(path):
    """
    Reads file at path, returns a list of charged nodes
    """
    chargedNodes = []
    
    with open(path, 'r') as f:
        for line in f:
            chargedNodes.append(int(line)-1)
    
    return chargedNodes


# In[ ]:

def runAndPrintExperiments(epsilon = 0.2, Gnames = ["ubiquitin", "bubble", "br", "fmo", "gfp"],
                           readChargedNodes = False, pathPrefix = "../../input/", graphSuffix = "_complete.graph",
                           chargeSuffix = "_charges.resid"):
    scores = []
    initialTime = time.time()
    
    for Gname in Gnames:
        G = readGraph(pathPrefix + Gname + graphSuffix, Format.METIS)
        chargedNodes = []#readCharges(pathPrefix + Gname + chargeSuffix)
        n = G.numberOfNodes()
        graphScores = []
        
        if n > 100:
            kList = [8,12,16,20,24]
        else:
            kList = [2,4,6,8]
        
        print("Graph:", Gname, "with", n, " nodes.")
        print("chargedNodes =", chargedNodes)
        for k in kList:
            if len(chargedNodes) > k:
                continue
            print("k = ", k)
            qualities = comparePartitionQuality(G, k, epsilon, chargedNodes, False)
            scores.append(qualities)
            graphScores.append(qualities)
            print('------------------------------------------------------------------')
        graphResults = {}
        for score in graphScores:
            for key in score:
                if not key in graphResults:
                    graphResults[key] = []
                graphResults[key].append(score[key])
                
        gMeans = {}
        for method in graphResults:
            gMeans[method] = functools.reduce(operator.mul, graphResults[method], 1) ** (1/len(graphResults[method]))
            print(method, ':', str(gMeans[method]))
        print("Ratio:", gMeans['bestOfFour'] / gMeans['naive'])
        print('##################################################################')
    
    results = {}
    for score in scores:
        for key in score:
            if not key in results:
                results[key] = []
            results[key].append(score[key])
            
    print("Elapsed Time:", time.time() - initialTime)
            
    print("Geometric Means:")
    for method in results:
        print(method, ':', str(functools.reduce(operator.mul, results[method], 1) ** (1/len(results[method]))))
        
    print("Arithmetic Means:")
    for method in results:
        print(method, ':', str(sum(results[method]) / len(results[method])))


# In[ ]:

def runAndLogExperiments(runs = 1, charges = False, epsilonList=[0.1,0.2]):
    pathPrefix = "../../input/"
    graphSuffix = "_complete.graph"
    chargeSuffix = "_charges.resid"
    algoList = ['ml', 'greedy', 'ka', 'naive', 'cont']
    Gnames = ["ubiquitin", "bubble", "br", "fmo", "gfp", "fmo"]
    #kList = [2**i for i in range(1,6)]
    maxIterations = 100
    
    for Gname in Gnames:
        G = readGraph(pathPrefix + Gname + graphSuffix, Format.METIS)
        potentiallyCharged = readCharges(pathPrefix + Gname + chargeSuffix)
        
        n = G.numberOfNodes()
        kList = []
        if n > 100:
            kList = [8,12,16,20,24]
        else:
            kList = [2,4,6,8]
        
        for run in range(runs):
            with open(Gname+'-results-'+str(run)+'.dat', 'w') as f:
                f.write('\t'.join(['k']+[str(e) for e in epsilonList]+['label\n']))

                for k in kList:
                    chargedNodes = []
                    if charges:
                        valid = False
                        i = 0

                        while not valid and i < maxIterations:
                            chargedNodes = random.sample(potentiallyCharged, int(k*0.8))
                            isCharged = [v in chargedNodes for v in G.nodes()]
                            valid = chargesValid(G, [k], min(epsilonList), isCharged)
                            i += 1
                        if not valid:
                            print("No valid charges found after "+str(i)+" iterations.")
                            continue

                        print(chargedNodes)

                    # data format: one line per algorithm and k
                    lineDict = {}
                    for algo in algoList:
                        lineDict[algo] = [str(k)]

                    for epsilon in epsilonList:
                        qualities = comparePartitionQuality(G, k, epsilon, chargedNodes, True)
                        for algo in algoList:
                            if algo in qualities:
                                lineDict[algo].append(str(qualities[algo]/qualities['naive']))
                            else:
                                lineDict[algo].append('NA')
                        print('Experiments done for k=', k, ', epsilon=', epsilon)

                    for i in range(len(algoList)):
                        algo = algoList[i]
                        lineDict[algo].append(str(i))
                        f.write('\t'.join(lineDict[algo])+'\n')


# In[ ]:

def averageLogs(runs, Gnames = ["ubiquitin", "bubble", "br", "fmo", "gfp"]):
    for Gname in Gnames:

        sumEntries = []
        numEntries = []

        for run in range(runs):
            filename = Gname+'-results-'+str(run)+'.dat'
            with open(filename, 'r') as f:
                f.readline()# remove header data
                lineNumber = 0
                for line in f:
                    lineList = line.split('\t')

                    if len(sumEntries) < lineNumber+1:
                        sumEntries.append([0 for field in lineList])
                        numEntries.append([0 for field in lineList])

                    for i in range(len(sumEntries[lineNumber]), len(lineList)):
                        sumEntries[lineNumber].append(0)
                        numEntries[lineNumber].append(0)

                    for i in range(len(lineList)):
                        sumEntries[lineNumber][i] +=  float(lineList[i])
                        numEntries[lineNumber][i] += 1

                    lineNumber += 1
        assert(len(sumEntries) == len(numEntries))

        outputname = Gname+'-results-averaged.dat'
        with open(outputname, 'w') as f:
            for rowIndex in range(len(sumEntries)):
                assert(len(sumEntries[rowIndex]) == len(numEntries[rowIndex]))
                linelist = [str(sumEntries[rowIndex][colIndex] / numEntries[rowIndex][colIndex]) for colIndex in range(len(sumEntries[rowIndex])) ]
                f.write('\t'.join(linelist)+'\n')

def writePartition(part, path):
    community.PartitionWriter().write(part, path)


def runAndSavePartitions(G, Gname, k = 8, epsilon = 0.2, isCharged = []):
    n = G.numberOfNodes()
    if len(isCharged) == 0:
        isCharged = [False for v in range(G.numberOfNodes())]
    
    sizelimit = int(math.ceil(n / k)*(1+epsilon))

    ml = mlPartition(G, k, epsilon, isCharged)
    ml = repairPartition(G, ml, epsilon, isCharged)
    ml.compact()
    result = ml
    resultWeight = partitioning.computeEdgeCut(result, G)
    writePartition(ml, 'MultiLevel-k-'+str(k)+'-imbalance-'+str(epsilon)+'-'+Gname+'.part')
    print("Wrote Multilevel partition with", ml.numberOfSubsets(), " fragments and weight", resultWeight)

    greedy = greedyPartition(G, k, epsilon, isCharged)
    cutWeight = partitioning.computeEdgeCut(greedy, G)
    if cutWeight < resultWeight:
        result = greedy
        resultWeight = cutWeight
    writePartition(greedy, 'Greedy-k-'+str(k)+'-imbalance-'+str(epsilon)+'-'+Gname+'.part')
    print("Wrote Greedy partition with", greedy.numberOfSubsets(), " fragments and weight", cutWeight)

    try:
        ka = kaHiPWrapper(G, k, epsilon)
        ka = repairPartition(G, ka, epsilon, isCharged)
        ka.compact()
        cutWeight = partitioning.computeEdgeCut(ka, G)
        if cutWeight < resultWeight:
            result = ka
            resultWeight = cutWeight
        writePartition(ka, 'KaHiP-k-'+str(k)+'-imbalance-'+str(epsilon)+'-'+Gname+'.part')
        print("Wrote KaHiP partition with", ka.numberOfSubsets(), " fragments and weight", cutWeight)
    except FileNotFoundError as e:
        pass

    try:
        cont = dpPartition(G, k, epsilon, isCharged)
        cutWeight = partitioning.computeEdgeCut(cont, G)
        if cutWeight < resultWeight:
            result = cont
            resultWeight = cutWeight
        writePartition(cont, 'DP-k-'+str(k)+'-imbalance-'+str(epsilon)+'-'+Gname+'.part')
        print("Wrote DP partition with", cont.numberOfSubsets(), " fragments and weight", cutWeight)
    except ValueError as e:
        pass

    try:
        contLegacy = dpPartition(G, k, epsilon, isCharged, True)
        cutWeight = partitioning.computeEdgeCut(contLegacy, G)
        writePartition(contLegacy, 'DP-legacy-k-'+str(k)+'-imbalance-'+str(epsilon)+'-'+Gname+'.part')
        print("Wrote legacy DP partition with", contLegacy.numberOfSubsets(), " fragments and weight", cutWeight)
    except ValueError as e:
        pass
    
    naive = naivePartition(G, k)
    naive = repairPartition(G, naive, 0, isCharged)
    cutWeight = partitioning.computeEdgeCut(naive, G)
    writePartition(naive, 'Naive-k-'+str(k)+'-'+Gname+'.part')
    print("Wrote naive partition with", naive.numberOfSubsets(), " fragments and weight", cutWeight)

    cutWeight = partitioning.computeEdgeCut(result, G)
    writePartition(result, 'Best-k-'+str(k)+'-'+Gname+'.part')
    print("Wrote selected partition with", result.numberOfSubsets(), " fragments and weight", cutWeight)
    

if __name__ == '__main__':
    filename = sys.argv[1]
    k = int(sys.argv[2])
    epsilon = float(sys.argv[3])
    chargedNodes = []
    for i in range(4, len(sys.argv)):
        chargedNodes.append(int(sys.argv[i]))
    
    G = readGraph(filename, Format.METIS)
    n = G.numberOfNodes()
    isCharged = [False for i in range(n)]
    for c in chargedNodes:
        assert(c < n)
        isCharged[c] = True

    Gname = os.path.basename(filename)

    runAndSavePartitions(G, Gname, k, epsilon, isCharged)