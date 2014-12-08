# pczsource.py
#
# Robert J. Brunner
# March 21, 2007
#
# This file contains the zsource class, which is used to hold an
# ra/dec/z tuple. It also includes a method that creates lists of
# zsources from input files, the read/write methods are inherited from
# the source class . These methods should work with subclasses that
# properly inherit from this class and also that override the parse
# class method.
#

from math import sin, cos, asin, acos, radians, degrees, sqrt
from gzip import GzipFile

from common.pcsource import pcsource
from cdistance import cdistance

class pczsource(pcsource):
    """

    This class represents the spatial coordinates for a source, and as
    such it inherits the angular coordinate data from the source class.
    
    v1.0 Robert J. Brunner, January 26, 2007

    v1.1 Robert J. Brunner, March 5, 2007

    Dropped angular diameter distance attribute as it is unnecessary.

    v1.2 Robert J. Brunner, March 15, 2007

    Changed precompute to use line-of-sight comoving coordinate only.

    v1.3 Robert J. Brunner, March 22, 2007

    Removed precomputed quantities from this class. spatial
    precomputes are now in sprecompute.py

    v1.4 Robert J. Brunner, March 27, 2007

    Added precomputed quantities back into this class. Since Natalie
    is no longer needing the Python code base, I can focus it back
    onto the task of quickly estimating correlations.

    """

  # General representation for calculating cosmological distances for all sources

    cd = cdistance(0.3, 0.7, 1.0)
    
    __slots__ = ['ra', 'dec', 'z', 'sra', 'cra', 'sdec', 'cdec', 'loscd']
    
    def __init__(self, ra, dec, z, sra, cra, sdec, cdec, loscd):

        pcsource.__init__(self, float(ra), float(dec), float(sra), float(cra), float(sdec), float(cdec))
        
        self.z = float(z)
        
        # Now set the precomputed quantities
        
        self.loscd = float(loscd)

    def create(ra, dec, z):
        """
        
        This static class method converts a ra/dec/z tuple in a string
        of data into an actual instance of the appropriate class.
    
        v1.0 Robert J. Brunner, March 26, 2007
 
        """

        # Must explicitly calculate the precomputed quantities.
        
        sra = sin(radians(ra))
        cra = cos(radians(ra))
        sdec = sin(radians(dec))
        cdec = cos(radians(dec))

        loscd = pczsource.cd.dm(z)
        
        return(pczsource(ra, dec, z, sra, cra, sdec, cdec, loscd))

    def parse(cls, data, racol = 0, deccol = 1, zcol = 2):
        """
        
        This classmethod converts a string of data into an actual
        instance of the appropriate class. By using a class method
        here, we allow subclasses to override this method, but still
        use the helper methods such as read/write sources.
    
        v1.0 Robert J. Brunner, March 21, 2007

        v1.1 Robert J. Brunner, March 27, 2007

        Changed to static method, and explictly extracted all columns.
        """
        
        ra = float(data[racol])         # The longitude coordinate
        dec = float(data[deccol])       # The latitude coordinate
        z = float(data[zcol])
        
        # We implicitly assume that sra, etc. follow dec column.

        sra = float(data[zcol + 1])
        cra = float(data[zcol + 2])
        sdec = float(data[zcol + 3])
        cdec = float(data[zcol + 4])        

        loscd = float(data[zcol + 5])

        return(cls(ra, dec, z, sra, cra, sdec, cdec, loscd))

    # This is the magic statement to make the methods true static methods
    
    create = staticmethod(create)
    parse = classmethod(parse)
    
    def __repr__(self):
        """

        In this method, we first use the string representation of the
        base class, to which we add the string representation of this
        class.
        
        The precomputed trig quanities need at least ten significant
        figures to maintain numerical stability. In this case,
        however, I am going with simplicity and dumping twelve, which
        is the default value of the trig calculations as floats.
        
        v1.0 Robert J. Brunner, January 26, 2007

        v1.1 Robert J. Brunner, March 21, 2007
        
        Modified to use parent source class repr method.

        v1.2 Robert J. Brunner, March 22, 2007
        
        Removed precomputed quantities from output.

        """
        
        return("%13.10f % 13.10f %6.4f %14.12f %14.12f %14.12f %14.12f %15.10f" %
               (self.ra, self.dec, self.z,
                self.sra, self.cra, self.sdec, self.cdec, self.loscd))

    def pcDeltaDistance(self, asource):
        """

        Calculates the line-of-sight comoving distance between two
        sources.
        
        INPUTS: A second source to use for the calculation
        
        OUTPUTS: The comoving line-of-sight distance between the two
        sources.
        
        v1.0 Robert J. Brunner, March 22, 2007
        """

        return abs(self.pcDistance() - asource.pcDistance())

    def pcDistance(self):
        """

        Returns the line-of-sight comoving distance for the given
        source.
        
        INPUTS: None
        
        OUTPUTS: The comoving line-of-sight distance to the source
        
        v1.0 Robert J. Brunner, March 22, 2007
        """

        return self.loscd

    def pcComovingDistance(self, asource):
        """

        Calculates the line-of-sight comoving distance between two
        sources.
        
        INPUTS: A second source to use for the calculation
        
        OUTPUTS: The comoving line-of-sight distance between the two
        sources.

        Note: This method simply takes the square root of the
        pcComovingDistance2 output.
        
        v1.0 Robert J. Brunner, March 22, 2007
        """

        return sqrt(self.pcComovingDistance2(asource))
    
    def pcComovingDistance2(self, asource):
        """

        Calculates the square of the line-of-sight comoving distance
        between two sources. This method is used since we can save
        square root calculations during correlation measurements by
        squaring the bin edges, and simply leaving the disstance
        squared.
        
        INPUTS: A second source to use for the calculation
        
        OUTPUTS: The comoving line-of-sight distance squared between
        the two sources.

        Note: This method uses all precomputed quantites and assumes a
        flat universe. For higher accuracy, the unprecomputed
        angularDistance method could be used to get Theta.
        
        v1.0 Robert J. Brunner, March 22, 2007
        """

        return (self.pcDistance()**2 + asource.pcDistance()**2 -
                2.0 * self.pcDistance() * asource.pcDistance() * self.pcCosTheta(asource))

    def createSources(infile, raColumn = 0, decColumn = 1, zColumn = 2, filterFunction = None):
        """
        
        Creates a new list of sources from a file containing right
        ascension, declination pairs. 
        
        INPUTS: Filename 
        
        OUTPUTS: None
    
        v1.0 Robert J. Brunner, February 23, 2007
    
        v1.1 Robert J. Brunner, March 19, 2007

        Modified to handle flexible file names and non Python 2.5 try
        blocks.

        v 1.2 Robert J. Brunner, March 20, 2007

        Modified to be class method on class zsource.
        
        """
    
        sources = []

        racol = int(raColumn)
        deccol = int(decColumn)
        zcol = int(zColumn)
        
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
                            z = float(cols[zcol])
                            
                            sources.append(pczsource.create(ra, dec, z))

            except IOError:
                print("ERROR: Could not read succesfully from ", infile)
                        
        finally:
            input.close()
                
        return sources
    
    # Note that we can't use filter as that is a Python keyword, hence the
    # lengthy use of filterFunction.

    def precomputeSources(infile, outfile, raColumn = 0, decColumn = 1, zColumn = 2,
                          filterFunction = None):
        """
        Creates a new file containing the precomputed data for sources
        from a file containing right ascension, declination pairs.
        
        INPUTS: Input filename, Output filename, ra/dec/z column
        ordinal numbers, filter function.
        
        OUTPUTS: None
        
        v1.0 Robert J. Brunner, March 19, 2007

        v1.1 Robert J. Brunner, March 21, 2007

        Added to pczsource class as class method. Modified to use parse
        class method.
        
        """
            
        racol = int(raColumn)
        deccol = int(decColumn)
        zcol = int(zColumn)
          
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
                            z = float(cols[zcol])
                            
                            output.write(repr(pczsource.create(ra, dec, z)) + "\n")

            except IOError:
                print "ERROR: IO problem when reading or writing from specified files."
                        
        finally:
            input.close()
            output.close()

    # Define as static methods.
    
    createSources = staticmethod(createSources)
    precomputeSources = staticmethod(precomputeSources)

# End of pczsource class definition
                 
# Test Code

if __name__ == '__main__':

    a = pczsource.create(10, 10, 0.1)
    b = pczsource.create(11, 11, 0.2)

    print(a.pcComovingDistance(b))
    print(a.pcAngularDistance(b))
    print(a.pcDeltaDistance(b))
    print("Distance to a = ", a.pcDistance())
    print(a)
    print(b)
