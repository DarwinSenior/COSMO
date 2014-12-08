
from math import sin, cos, degrees, radians

# Globally useful functions:

cumulate = lambda x, y: (x + y)

dra = lambda a, b, c: (a - b) * cos(radians(c))

sparea = lambda r1, r2, d1, d2: (degrees(sin(radians(d2)) - sin(radians(d1))) * (r2 - r1))

# This is a factory method to create Sky Nodes

def nodeFactory(aClass, *args):
    return apply(aClass, args)

class skyNode:
    """The class that represents all 2d Sky Nodes in the quad tree"""

    raOffset = 0
    
    minAngle = radians(1.0/3600.0) # One arcsecond in radians
    
    def __init__(self, minra, maxra, mindec, maxdec):
        
        # get cenral ra, dec. These formulas will not work with  0/360 RA crossing.

        self.deltaRA = 0.5 * (maxra - minra)
        self.ra = maxra - self.deltaRA

        self.deltaDEC = 0.5 * (maxdec - mindec)
        self.dec = maxdec - self.deltaDEC

        self.count = 0
        self.area = 0.0 
        self.density = 0.0 
        
    # Change to def __rep__ invocation
    
    def __repr__(self):

        output = "%8d %12.6f %6.2f " % (self.count, self.area, self.density)
        output = output + "[%12.6f:%12.6f," % ((self.ra - self.deltaRA), (self.ra + self.deltaRA))
        output = output + "%12.6f:%12.6f]\n" % ((self.dec - self.deltaDEC), (self.dec + self.deltaDEC))

        return output

    def putSource(self, source):
        pass # Method should be overridden

    def findSources(self, target, angle, cosAngle):
        pass # Must be overridden

    def getStats(self):
        pass # Must be overridden

    
##  # Given a source, this function returns true if the node contains
##  # the source, false otherwise
    
##  def inNode(self, target):

##  if((target.ra + skyNode.raOffset < self.ra - self.deltaRA) or
##     (target.ra + skyNode.raOffset > self.ra + self.deltaRA) or
##     (target.dec < self.dec - self.deltaDEC) or
##     (target.dec > self.dec + self.deltaDEC)):

##    return False

##  return True

# Minus one in the init method is used to signal no number was passed
# in, which should only happen for the root node. The root node should
# neever have its number accessed so it should not be a normal value
# in order to trigger an error condition.

class skyLeafNode(skyNode):
    """The class representing a 2d sky leaf node in the quad tree"""

    def __init__(self, level, minra, maxra, mindec, maxdec, number = -1):

        skyNode.__init__(self, minra, maxra, mindec, maxdec)

        self.sources = [] # empty list of sources.
        self.number = number

        self.area = sparea(minra, maxra, mindec, maxdec)

    def __repr__(self):
        return "L " + skyNode.__repr__(self)

    def putSource(self, source):
        self.sources.append(source)

    def findSources(self, target, dangle):

        # Convert degrees angle into radians angle
        
        rangle = radians(dangle)
        
        # not test is used to eliminate testing for less than or equal
        # (which we now is a problem with floats)

        sources = []
        
        for source in self.sources:
            
            stAngle = source.pcAngularDistance(target)

            # The angle must be greater than some minimum in order to
            # prevent identity comparisons, and also must be less than
            # or equal to deasired angle.
            
            if not (stAngle > rangle) and (stAngle > skyNode.minAngle):
                sources.append(source)

        return sources

    def getStats(self):

        self.count = len(self.sources)
        self.density = float(self.count) / self.area

# The four child nodes in a non-leaf node are laid out as follows:
#
#   2   3
#   0   1
#

class skyNonLeafNode(skyNode):
    """The class representing a 2d sky non-leaf node in the quad tree"""

    def __init__(self, level, minra, maxra, mindec, maxdec, number = -1):
        skyNode.__init__(self, minra, maxra, mindec, maxdec)

        self.nodes = []
        self.number = number
        
        level -= 1
        
        if(level > 0):
            nodeType = skyNonLeafNode
        else:
            nodeType = skyLeafNode
        
        self.nodes.append(nodeFactory(nodeType, level, minra, self.ra, mindec, self.dec, 0))
        self.nodes.append(nodeFactory(nodeType, level, self.ra, maxra, mindec, self.dec, 1))
        self.nodes.append(nodeFactory(nodeType, level, minra, self.ra, self.dec, maxdec, 2))
        self.nodes.append(nodeFactory(nodeType, level, self.ra, maxra, self.dec, maxdec, 3))

        self.area = reduce(cumulate, [node.area  for node in self.nodes])

    def __repr__(self):

        output = "N " + skyNode.__repr__(self)
        
        for node in self.nodes:
            output += node.__repr__()

        return output

    def putSource(self, target):
        
        number = self.__findQuadrant(target)
	
	self.nodes[number].putSource(target)

    def findSources(self, target, angle):
        
        sources = []
       
        nodes = [] 

        # First find which quadrant contains the point. It needs to be
        # searched completely

        number = self.__findQuadrant(target)
	
	nodes.append(self.nodes[number])

        # if number < 2, we need to search dec + deltaDEC
        # otherwise its minus.

        if(number < 2):
            # we need to search number + 1 for dec only
            if (angle > (self.dec - target.dec)):
                nodes.append(self.nodes[number + 2])
        else:
            if(angle > (target.dec - self.dec)):
                nodes.append(self.nodes[number - 2])
        
        # if number % 2 we need to do ra + deltaRA, otherwise its
        # minus.

        if((number % 2) == 0):
            if (angle > dra(self.ra, target.ra + skyNode.raOffset, target.dec)):
                nodes.append(self.nodes[number + 1])
        else:
            if(angle > dra(target.ra + skyNode.raOffset, self.ra, target.dec)):
                nodes.append(self.nodes[number - 1])

        # If 3 nodes were found in search area, then fourth also must
        # be searched.  This could be tested more carefully if doing
        # circle intersection, but when using square intersection this
        # is true.
        
        if(len(nodes) == 3) :
            nodes = self.nodes

        # No need to sort
        # nodes.sort()

        # Could add test to verify node has data prior to searching it
        # Might not need += below, maybe just =
        #
        #sources += [node.findSources(target, angle, cosAngle) for node in nodes \
        #       if node.count > 0]

        for node in nodes:
            if node.count > 0:
                sources += node.findSources(target, angle)

        return sources

    def __findQuadrant(self, target):
        if(target.ra + skyNode.raOffset < self.ra):
            if(target.dec < self.dec):
                return 0
            else:
                return 2
        else:
            if(target.dec < self.dec):
                return 1
            else:
                return 3

    def getStats(self):
        
        for node in self.nodes:
            node.getStats()
            self.count += node.count
            
        self.density = float(self.count) / self.area

class skyQuadTree:

    """The class that represents all 2d sky nodes in the 2d sky quad tree"""
    
    level = 0       # The level to which the 2d sky quad tree should be built.
    
    def __init__(self, level, minra, maxra, mindec, maxdec, raOffset = 0):

        skyQuadTree.level = level
        self.minra = minra + raOffset
        self.maxra = maxra + raOffset
        self.mindec = mindec
        self.maxdec = maxdec
        
        skyNode.raOffset = raOffset
        
        self.root = nodeFactory(skyNonLeafNode, level, \
                                self.minra, self.maxra, self.mindec, self.maxdec)
    def fillTree(self, sources):
        
        for source in sources:
            self.root.putSource(source)

        self.root.getStats()

    def __repr__(self):
        return self.root.__repr__()

    def findSources(self, source, angle):

        # angle must be in degrees at this stage in order to handle
        # node tests properly. But it has to be converted to radians
        # in order to test with angularDistance on a source (or
        # dervied class base on source) class object.
        
        sources = self.root.findSources(source, angle)
        
        return sources

def theMain():
    import sys

    a = skyQuadTree(4, -40, -30, -35, -25)

    sources = []

    for line in open('/Users/rb/Desktop/example.input').readlines():
        cols = line.split()
        sources.append(source(float(cols[0]) / 60.0, float(cols[1])/60.0))

    a.fillTree(sources)

    print(a)

    limit = 1.0/60.0

    for testsource in sources:

        matches = a.findSources(testsource, limit)

        output = ""

        for s in matches:
            output = output + "%6.4f " % (s.angularDistance(testsource))

        print(len(matches), testsource.ra, testsource.dec, output)

# On My Laptop
## No raOffset in code
## Total Time     :     221.11 
#
## raOffset = 0.0 
## Total Time     :     245.58 
#
## raOffset = 90.0
## Total Time     :     238.62 

# Test Code

if __name__ == '__main__':
    
    from time import clock

    startTime = clock()

    theMain()

    print("Total Time     :", "%12.2f " % (clock() - startTime))
