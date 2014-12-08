from math import sqrt
from operator import itemgetter
from gzip import GzipFile

from math import degrees, atan2, asin

# Globally useful functions:

def cumulate(data):
    return reduce(lambda x, y: x + y, data)

def stats(data):
    mean = cumulate(data)/float(len(data))
    stdev = sqrt(cumulate([(i - mean)**2 for i in data])/float(len(data)))

    return (mean, stdev)

# Not needed to ****

# This function does a weighted least squares on the points x, y with
# weigths given by s. The returned data is the intercept and the
# slope (i.e., (a, b))

def wlsq(x, y, s):
    w = [1/(i**2) for i in s]
    s = cumulate(w)
    sx = cumulate([x[i] * w[i] for i in range(len(x))])
    sy = cumulate([y[i] * w[i] for i in range(len(y))])
    sxx = cumulate([x[i] * x[i] * w[i] for i in range(len(x))])
    sxy = cumulate([x[i] * y[i] * w[i] for i in range(len(x))])

    delta = s * sxx - sx**2
    return ((sxx * sy - sx * sxy)/delta, (s * sxy - sx * sy) / delta)

# *****

# This is a factory method to create Jackknife Nodes

def nodeFactory(aClass, *args):
    return apply(aClass, args)

class jkNode:
    """The class that represents all Jackknife Nodes in the tree"""

    nodeNumber = 0

    def __init__(self, level):

        self.nodeNumber = jkNode.nodeNumber
        jkNode.nodeNumber += 1

        self.level = level
        
    def getJackknife(self, source):
        pass # Must be overridden in subclass

class jkLeafNode(jkNode):
    """The class representing a Jackknife Leaf Node in the  tree"""

    id = 0
    
    def __init__(self, level, sources):

        jkNode.__init__(self, level)

        # Now set id for jackknife numerical assignment

        self.id = jkLeafNode.id        
        jkLeafNode.id += 1

#        print "%d source sin leaf node #%d" % (len(sources), self.id)

    def __repr__(self):

        output = "Leaf" 

        output += "Node Number = %d\n" % (self.nodeNumber) 
        output += "Node ID = %d\n" % (self.id) 
               
        return output

    def getJackknife(self, pcsource):
        return self.id

class jkNonLeafNode(jkNode):
    """The class representing a Jackknife Non Leaf Node in the tree"""

    def __init__(self, level, sources):
        jkNode.__init__(self, level)

        level -= 1

        splitColumn = -1
        splitValue = 0.0

        if(level > 0):
            nodeType = jkNonLeafNode
        else:
            nodeType = jkLeafNode

        splitIndex = self.__getSplitIndex(sources)
        
        self.left = nodeFactory(nodeType, level, sources[:splitIndex])
        self.right = nodeFactory(nodeType, level, sources[splitIndex:])

#        self.area = reduce(cumulate, [node.area  for node in self.nodes])

    def __repr__(self):

        output = "Non Leaf"
        output += "Node Number = %d\n" % (self.nodeNumber)
        output += "Split Column = %d\n" % (self.splitColumn)
        output += "Split Value = %12.6f\n" % (self.splitValue)
        
        output += self.left.__repr__()
        output += self.right.__repr__()

        return output

# Go through list, splitting out ra/dec, compute mean in both for one pass
# then get min./max of ra list and dec list, and go through to get sigma of each.

    def __getSplitIndex(self, sources):

        # Get Statistics

        (meanX, sigmaX) = stats(map(itemgetter(0), sources))
        (meanY, sigmaY) = stats(map(itemgetter(1), sources))
        (meanZ, sigmaZ) = stats(map(itemgetter(2), sources))

        if sigmaX > sigmaY:
            if sigmaX > sigmaZ:
                self.splitColumn = 0 #i.e., X Cartesian Coordinate
                self.splitValue = meanX
            else:
                self.splitColumn = 2 #i.e., Z Cartesian Coordinate
                self.splitValue = meanZ

        elif sigmaY > sigmaZ:
            self.splitColumn = 1 #i.e., Y Cartesian Coordinate
            self.splitValue = meanY
        else:
            self.splitColumn = 2 #i.e., Z Cartesian Coordinate
            self.splitValue = meanZ

        sources.sort(key=itemgetter(self.splitColumn)) ;

        index = int(0.5 * len(sources)) # Start at midpoint value

        if self.splitValue < sources[index][self.splitColumn]:
            
            #Start at midpoint and go down until find index for mean
            
            for source in reversed(sources[:index]):

                if source[self.splitColumn] < self.splitValue:
                    break
                else: 
                    index -= 1
        else:
            for source in sources[index:]:

                if source[self.splitColumn] > self.splitValue:
                    break
                else:
                    index += 1
        
        return index

    def getJackknife(self, pcsource):
        
        if pcsource.getCartesianValue(self.splitColumn) < self.splitValue:
            return self.left.getJackknife(pcsource)
        else:
            return self.right.getJackknife(pcsource)

class jkTree:

    """The class that represents the Jackknife Tree"""
    
    level = 0            # The level to which the Jackknife Tree should be built.
    
    def __init__(self, level, sources):

        if level > 0:

            jkTree.level = level
        
            self.root = nodeFactory(jkNonLeafNode, level, sources)
            
        else:
            self.root = None

    def __repr__(self):
        return self.root.__repr__()

    def getJackknife(self, pcsource):
        return self.root.getJackknife(pcsource)

def theMain(inFile, treeLevel, raCol, decCol):

    from common.pcsource import pcsource

    sources = pcsource.createSources(inFile, raCol, decCol)

    cSources = []

    for source in sources:
        cSources.append((source.getCartesianValue(0), 
                        source.getCartesianValue(1), 
                        source.getCartesianValue(2)))

    for source in cSources:
	print(source[0], source[1], source[2])

    jktree = jkTree(treeLevel, cSources)

#    for source in sources:
#        print source.ra, source.dec, jktree.getJackknife(source)
	
#    print jktree

#    jackknifeNumbers = [jktree.getJackknife(source) for source in sources]
#
#    for index in range(8):
#        temp = [value for value in jackknifeNumbers if value == index]
#        print index, len(temp)

if __name__ == '__main__':
    
    from time import clock

    startTime = clock()

    from sys import argv
    
    inFile = argv[1]
    treeLevel = int(argv[2])
    raCol = int(argv[3])
    decCol = int(argv[4])

    theMain(inFile, treeLevel, raCol, decCol)

#    print "Total Time        :", "%12.2f " % (clock() - startTime)

