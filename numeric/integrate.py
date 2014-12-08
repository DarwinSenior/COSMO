#
# Basic Extended trapezoid integrator. This class can be subclassed by other
# integration routines that rely on the basic trapezoid rule, such as the
# simpson or romberg integration routines.
#
# Robert J. Brunner
#
# March 5, 2007
#

class trapezoid:
    
    """
    This class .
    
    v1.0 Robert J. Brunner, March 5, 2007
    """

    initialResult = -1E99
    epsilon = 1E-6
    maxIterations = 25
    
    def __init__(self, integrand):
        
        self.numPoints = 0  # Number of points to sample function on next iteration
        self.result = trapezoid.initialResult  # Current value of integral

        self.integrand = integrand
        self.start = 0.0
        self.end = 0.0
        
    def integrate(self, start, end):

        """
        Calcualtes the .
    
        INPUTS: 
        
        OUTPUTS: 
        
        This method returns the 
        
        v1.0 Robert J. Brunner, March 5, 2007
                
        v1.1 Robert J. Brunner, March 14, 2007

        Added test for limits being equal
        """

        # Do nothing if limits are equal. Needed for comsological
        # calculations when z = 0
        
        if start == end: result = 0.0
        
        else:
            self.numPoints = 0   
            
            self.start = float(start)
            self.end = float(end)
            
            oldResult = trapezoid.initialResult
            
            for iteration in range(trapezoid.maxIterations):
                
                result = self.theWork(iteration)
                
                if(abs(result - oldResult) < trapezoid.epsilon * abs(oldResult)):
                    break
                
                oldResult = result
                
            else:
                
                print "Error, integration did not converge when using extended trapezoid rule."
                
        return result
    
    def theWork(self, iteration):
       
        """
        Does the real work by implementing the extended trapezoid rule.
        
        INPUTS: A 
        
        OUTPUTS: The 
        
        v1.0 Robert J. Brunner, March 5, 2007     
        """
        
        if(iteration == 0):
            self.numPoints = 1
            self.result = 0.5 * (self.end - self.start) * \
                          (self.integrand(self.start) + self.integrand(self.end))

        elif (iteration > 0):
            delta = (self.end - self.start) / self.numPoints
            
            sum = reduce(lambda x, y: x + y, \
                         [self.integrand(self.start + 0.5 * delta + i * delta) \
                          for i in range(self.numPoints)])
            
            self.result = 0.5 * (self.result + (self.end - self.start) * sum / self.numPoints)
            
            self.numPoints *= 2

        return self.result

def testfunc(x):
    return (x)

def testfunc2(x):
    from math import sqrt, pow
    return (1.0 / sqrt(0.3 * pow((1.0 + float(x)), 3) + 0.7))
            
# Test Code

if __name__ == '__main__':

    t = trapezoid(testfunc)

    print t.integrate(0.0, 1.0), "should be 0.5"

    t2 = trapezoid(testfunc2)

    print 3000.0 * (1 + float(0.558)) * t2.integrate(0.0, 0.558),
    print 2265.876909, "NW's tool",  2264.797863, "Romberg"
