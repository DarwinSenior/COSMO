# pcsource.py
#
# Robert J. Brunner
# March 26, 2007
#

# This file contains the pcsource class, which is used to hold an ra/dec
# pair, along with the precomputed quantities. It also includes a
# number of methods that read/write and create lists of sources from
# input files. These methods should work with subclasses that properly
# inherit from this class and also that override the parse class
# method.
#
# This class is a redefinition of old classes developed earlier.
#

# Note that there is a finite accuracy of the use of precomputed
# quantities in this class:

# p1 = source.create(113.663293553, 32.0010041224)
# p2 = source.create(113.663293554, 32.001004123)
# p1.pcAngularDistance(p2)
#     2.1073424255447017e-08
# p1.pcAngularDistance(p1)
#     2.1073424255447017e-08
#

# This obviously limits the accuracy of angle measurements when using
# the presource class to about 4.3 milliarcseconds.In practice it is
# even worse, which makes the haversine formula a must when getting
# the distance between close points. The limiting accuracy needs to be
# closer to 0.5 arcseconds.

from math import sin, cos, asin, acos, radians, degrees, sqrt
from gzip import GzipFile

class pcsource:
    
    """
    
    This class represents the angular coordinates for a source, and
    now includes the precomputation of sine/cosine of the two angular
    coordinates, which provides performnce boosts when doing
    correlation measurements.

    Note: No default values are allowed in the constructor anymore.
    
    v1.0 Robert J. Brunner, March 26, 2007

    v1.1 Robert J. Brunner, February 14, 2008
    
    Added cartesian coordinate getter methods.

    """

    __slots__ = ['ra', 'dec', 'sra', 'cra', 'sdec', 'cdec']
    
    def __init__(self, ra, dec, sra, cra, sdec, cdec):
        
        self.ra = float(ra)         # The longitude coordinate
        self.dec = float(dec)       # The latitude coordinate

        # Now set the precomputed quantities
        
        self.sra = float(sra)
        self.cra = float(cra)
        self.sdec = float(sdec)
        self.cdec = float(cdec)

    def __str__(self):
        """

        By placing this method in the base class, each subclass will
        inherit it and call the appropriate repr method.
        
        v1.0 Robert J. Brunner, March 26, 2007
        """
        return(self.__repr__())

    def __repr__(self):
        """
        
        Convert this object to a string. The precision has been tested
        to be sufficient for existing needs.
        
        v1.0 Robert J. Brunner, March 26, 2007
 
        """
        
        return("%13.10f % 13.10f %14.12f %14.12f %14.12f %14.12f" %
               (self.ra, self.dec, self.sra, self.cra, self.sdec, self.cdec))

    def create(ra, dec):
        """
        
        This static class method converts a ra/dec pair in a string of data
        into an actual instance of the appropriate class.
    
        v1.0 Robert J. Brunner, March 26, 2007
 
        """
        
        # Must explicitly calculate the precomputed quantities.
        
        sra = sin(radians(ra))
        cra = cos(radians(ra))
        sdec = sin(radians(dec))
        cdec = cos(radians(dec))

        return(pcsource(ra, dec, sra, cra, sdec, cdec))
        
    def parse(cls, data, racol = 0, deccol = 1):
        """
        
        This static class method converts a string of data into an actual
        instance of the appropriate class.
    
        v1.0 Robert J. Brunner, March 26, 2007
 
        """        
        ra = float(data[racol])         # The longitude coordinate
        dec = float(data[deccol])       # The latitude coordinate

        # We implicitly assume that sra, etc. follow dec column.

        sra = float(data[deccol + 1])
        cra = float(data[deccol + 2])
        sdec = float(data[deccol + 3])
        cdec = float(data[deccol + 4])

        return(cls(ra, dec, sra, cra, sdec, cdec))

    # This is the magic statement to make the method a truestatic method
    
    create = staticmethod(create)
    parse = classmethod(parse)

    def hsAngularDistance(self, aSource):
       
        """

        Calculates the angle between two sources by using the
        haversine formula.
        
        INPUTS: A second source to use for the calculation
        
        OUTPUTS: The angle in radians
        
        This method returns the angle in radians between two sources
        in the sky. The formula is derived by taking the dot product
        of the two vectors that terminate at the ra/dec pairs. To use
        this calculation, for example you might compare it to the
        ratio of the projected physical distance to the angular
        diameter distance:
        
        d/da

        Algorithmically the calculated formula is usually written as
        (note the conversion from normal Spherical coordinates to
        equatorial cordinate flips cos and sin terms when dealing with
        declination):
        
        theta = acos(cos(dec1) * cos(dec2) * cos(ra2-ra1) + sin(dec1) * sin(dec2))

        But the cos(ra2 - ra1) term can result in numerical
        instabilities for small ra differences. In this class, with no
        precoputed quantities, we actually use the Haversine formula
        (Try Google ;-) to get a numerically stable version.

        angle = 2*asin(min(1,sqrt(sin(deltaDEC)**2 + cos(dec1) *
        cos(dec2) * sin(deltaRA)**2)))
        
        v1.0 Robert J. Brunner, March 20, 2007
        
        Note that this idea goes back to much earlier code that I
        wrote for group finding in 2004.
        
        """
                
	ra1 = radians(self.ra)
	ra2 = radians(aSource.ra)
	dec1 = radians(self.dec)
	dec2 = radians(aSource.dec)
	
	deltaRA = 0.5 * (ra1 - ra2)     # delta longitude
	deltaDEC = 0.5 * (dec1 - dec2)  # delta latitude

	# Angular distance in radians
        
	angle = 2*asin(min(1.0, sqrt(sin(deltaDEC)**2 + cos(dec1) * cos(dec2) * sin(deltaRA)**2)))

        return angle
    
    def pcCosTheta(self, asource):

        """

        Calculates the cosine of the angle between two sources using
        precomputed data.        
    
        INPUTS: A second source to use for the calculation
        
        OUTPUTS: Cosine of the angle
        
        This method returns the cosine of the angle between two
        sources in the sky. The formula is derived by taking the dot
        product of the two vectors that terminate at the ra/dec
        pairs. Note that the result is left as a cosine to avoid
        repeated calls to the arcosine function. To properly handle
        this result, it is likely that any comparisons be tranformed
        into the cosine space as well. For example, by taking the
        cosine of the ratio of the projected distance to the angular
        diameter distance:
        
        cos(radians(d/da))

        Algorithmically the calculated formula is usually written as
        (note the conversion from normal spehrical coordinates to
        equatorial cordinate flips cos and sin terms when dealing with
        declination):
        
        cos(theta) = cos(dec1) * cos(dec2) * cos(ra2-ra1) + sin(dec1) * sin(dec2)

        But the cos(ra2 - ra1) term can result in numerical
        instabilities for small ra differences. So that term is
        expanded out using the cos(A - B) = cosA cosB + sinA sinB
        rule. This also allows us to use the precomputed trig function
        values thereby speeding up the calculation significantly at
        this stage, which is now simply additions and multiplications.

        NOTE: This formula is not numerically exact, such as the
        haversine formula, which can be verified by calling this
        method to find the angle between the same object. For an exact
        representation, you should use the haversine formula, or
        perhaps add a test for identity (on the source index value).
        
        v1.0 Robert J. Brunner, January 26, 2007
        
        v1.1 Robert J. Brunner, February 28, 2007
        
        FIXED: Changed cos(dec) to sin(dec) and vice versa to fix my
        spherical coordinate bug.
        
        """
        
        return(self.cdec * asource.cdec * (self.cra * asource.cra + self.sra * asource.sra) + \
               self.sdec * asource.sdec)
    
    def pcAngularDistance(self, asource):
       
        """

        Calculates the angle between two source in radians by using
        precomputed data
        
        INPUTS: A second source to use for the calculation
        
        OUTPUTS: The angle in radians
        
        To simplify the calculation this method relies on the cosTheta
        member function for the heavy lifting, where the algorithm
        details are presented.
        
        v1.0 Robert J. Brunner, March 14, 2007

        v1.1 Robert J. Brunner, March 21, 2007

        Renamed methods to capitalize on the fact we generally want
        the result in radians not degrees.

        """
        
        return acos(self.pcCosTheta(asource))

    def pcAngularDistanceDegrees(self, asource):
       
        """
        
        Calculates the angle between two sources in degrees by using
        precomputed data
        
        INPUTS: A second source to use for the calculation
        
        OUTPUTS: The angle in degrees
        
        To simplify the calculation this method relies on the cosTheta
        member function for the heavy lifting, where the algorithm
        details are presented.
        
        v1.0 Robert J. Brunner, January 26, 2007

        v1.1 Robert J. Brunner, March 14, 2007

        Modified to call new angularDistanceRadians method

        v1.2 Robert J. Brunner, March 21, 2007

        Renamed methods to capitalize on the fact we generally want
        the result in radians not degrees.
        
        """
        
        return degrees(self.pcAngularDistance(asource))

    def getCartesianValue(self, index):
        """
        
        Astronomical coordinate data can be transformed from spheircal
        coordinates to a traditional Cartesian coordinate system. This method
        utilizes the precomputed trig functions to rapidly calculate the
        appropriate Carteisan coordinate requested. Note that the X, Y, Z
        coordinates are refernced as 0, 1, 2, respectively.

        INPUTS: Coordinate index to return 
        
        OUTPUTS: Requested coordinate value.
    
        v1.0 Robert J. Brunner, February 14, 2008
        """
        
        if index == 0: # X coordinate value
            return self.cdec * self.cra

        elif index == 1: # Y coordinate value
            return self.cdec * self.sra

        elif index == 2: # Z coordinate value
            return self.sdec

        # Else we have a bad index, and must indicate an error condition

        return 1E12

    def createSources(infile, raColumn = 0, decColumn = 1, filterFunction = None):
        """
        
        Creates a new list of sources from a file containing right
        ascension, declination pairs. 
        
        INPUTS: Filename 
        
        OUTPUTS: None
    
        v1.0 Robert J. Brunner, February 23, 2007
    
        v1.1 Robert J. Brunner, March 19, 2007

        Modified to handle flexible file names and non Python 2.5 try
        blocks.

        v 1.2 Robert J. Brunner, March 26, 2007

        Modified to be static method on class source.
        
        """
    
        sources = []

        racol = int(raColumn)
        deccol = int(decColumn)
        
        # In Python 2.5 the wrapper try block is not needed, but for
        # those who have not yet upgraded, this hack is required.

        try: # For the Finally clause
            
            try: # For the except clause

                # We open the file as a gzip'd file if suffix is .gz, else
                # treat as a normal file
        
                if infile.endswith(".gz"):
                    input = GzipFile(infile, 'rb')
                else:
                    input = open(infile)
            
                for line in input:

                    if (line.startswith('#') == False):
            
                        cols = line.split()

                        if filterFunction == None or filterFunction(cols):
                            
                            ra = float(cols[racol])         # The longitude coordinate
                            dec = float(cols[deccol])       # The latitude coordinate

                            sources.append(pcsource.create(ra, dec))

            except IOError:
                print "ERROR: Could not read succesfully from ", infile
                        
        finally:
            input.close()
                
        return sources

    def readSources(cls, infile):
        """
    
        Creates a new list of source objects from specially formatted,
        gzip'd file.
        
        INPUTS: Filename containing sources with no .gz extension
        
        OUTPUTS: List of source objects
        
        v1.0 Robert J. Brunner, February 23, 2007
        
        v1.1 Robet J. Brunner, March 19, 2007
        
        Modified to handle flexible file names and non Python 2.5 try
        blocks.

        v 1.2 Robert J. Brunner, March 20, 2007

        Modified to be class method on class source.
            
        """
    
        sources = []
    
        # In Python 2.5 the wrapper try block is not needed, but for
        # those who have not yet upgraded, this hack is required.

        try: # For the Finally clause
            
            try: # For the except clause

                # We open the file as a gzip'd file if suffix is .gz,
                # else treat as a normal file
        
                if not infile.endswith(".gz"):
                    infile += ".gz"

                input = GzipFile(infile, 'rb')
        
                for line in input:

                    cols = line.split()
            
                    sources.append(cls.parse(cols))
        
            except IOError:
                print "ERROR: Could not read succesfully from ", infile

        finally:
            input.close()
             
        return(sources)
    
    def writeSources(cls, outfile, sources):
        """
    
        Writes a list of source objects to specified compressed file using
        a specially format.
        
        INPUTS: Source list and Filename (with no .gz extension)
        
        OUTPUTS: None
        
        v1.0 Robert J. Brunner, February 23, 2007
        
        v1.1 Robet J. Brunner, March 19, 2007
        
        Modified to handle flexible file names and non Python 2.5 try
        blocks.

        v 1.2 Robert J. Brunner, March 20, 2007

        Modified to be class method on class source.
                    
        """
    
        # In Python 2.5 the wrapper try block is not needed, but for
        # those who have not yet upgraded, this hack is required.

        try: # For the Finally clause
            
            try: # For the except clause
        
                # If outfile does not end with .gz, we add it on to
                # the filename
        
                if outfile.endswith(".gz") == False:
                    outfile += ".gz"
            
                output = GzipFile(outfile, 'wb')

                for source in sources:
                    output.write(repr(source) + "\n")

            except IOError:
                print "ERROR: Could not write succesfully to ", outfile

        finally:
            output.close()
    
        return None

    # Note that we can't use filter as that is a Python keyword, hence the
    # lengthy use of filterFunction.

    def precomputeSources(infile, outfile, raColumn = 0, decColumn = 1, filterFunction = None):
        """
        Creates a new file containing the precomputed data for sources
        from a file containing right ascension, declination pairs.
        
        INPUTS: Input filename, Output filename, ra/dec column ordinal
        numbers, filter function.
        
        OUTPUTS: None
        
        v1.0 Robert J. Brunner, March 19, 2007

        v1.1 Robert J. Brunner, March 21, 2007

        Added to pcsource class as static method. Modified to use parse
        class method.
        
        """
            
        racol = int(raColumn)
        deccol = int(decColumn)
          
        # In Python 2.5 the wrapper try block is not needed, but for
        # those who have not yet upgraded, this hack is required.
        
        try: # For the Finally clause
            
            try: # For the except clause

                # We open the file as a gzip'd file if suffix is .gz, else
                # treat as a normal file
        
                if infile.endswith(".gz"):
                    input = GzipFile(infile, 'rb')
                else:
                    input = open(infile)
        
                # If outfile does not end with .gz, we add it on to the filename
        
                if outfile.endswith(".gz") == False:
                    outfile += ".gz"
                
                output = GzipFile(outfile, 'wb')

                for line in input:
                
                    if (line.startswith('#') == False):
                
                        cols = line.split()
            
                        if filterFunction == None or filterFunction(cols):

                            ra = float(cols[racol])         # The longitude coordinate
                            dec = float(cols[deccol])       # The latitude coordinate

                            output.write(repr(pcsource.create(ra, dec)) + "\n")

            except IOError:
                print "ERROR: IO problem when reading or writing from specified files."
                        
        finally:
            input.close()
            output.close()

    # This is the magic statement to make the methods true "class methods"

    readSources = classmethod(readSources)
    writeSources = classmethod(writeSources)
    
    # This is the magic statement to make the methods true "static methods"

    createSources = staticmethod(createSources)
    precomputeSources = staticmethod(precomputeSources)
    
# End of pcsource class definition

# Test Code

if __name__ == '__main__':

    a = pcsource.create(10, 10)
    b = pcsource.create(11,11)

    print a.hsAngularDistance(b)
    print a.pcCosTheta(b)
    print a.pcAngularDistance(b)
    print a
    print b
    
    # To sort a list of source objects, simply call the list's sort
    # method but pass in a lambda function defined such that
    #
    # myList.sort(lambda x, y: cmp(x.dec,y.dec))
    #
    # which will sort the list's elements by declination, which is suitable
    # for most applications as it avoids the cos(dec) term that is necessary
    # for ra comparisons.

    p1 = pcsource.create(113.663293553, 32.0010041224)
    p2 = pcsource.create(113.663293554, 32.001004123)
    print "%20.16E" % (p1.pcAngularDistance(p2))
    print "%20.16E" % (p1.pcAngularDistance(p1))
    print "%20.16E" % (p1.hsAngularDistance(p1))
