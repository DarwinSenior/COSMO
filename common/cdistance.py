#
# Robert J. Brunner
#
# A collection of cosmology distance calculations based on the details
# presented in the unpublished paper by David Hogg, astro-ph/9905116
#

from math import sqrt, pi
from numeric.integrate import trapezoid

class cdistance:

    def __init__(self, om = 0.3, ol = 0.7, h0 = 1.0):
        self.om = float(om)
        self.ol = float(ol)
        self.h0 = float(h0)

        self.iRoutine = trapezoid(self.oneOverEvolution)
          
    def evolution(self, z):
        """Calculates E(z) for a given redshift assuming standard 
        concordance cosmology
        
        INPUTS: Redshift, z
        
        OUTPUTS: E(z)
        
        v1.0 Robert J. Brunner, Jan. 18, 2007
        """
        
        return sqrt(self.om * pow((1.0 + float(z)), 3) + self.ol)

    def oneOverEvolution(self, z):
        """Calculates 1.0/E(z), which is needed for distance calculations
        
        INPUTS: Redshift, z
        
        OUTPUTS: 1.0/E(z)
        
        v1.0 Robert J. Brunner, Jan. 18, 2007
        """
        
        return 1.0/self.evolution(float(z))

    def dm(self, z):
        """Calculates line-of-sight comoving distance
        
        INPUTS: Redshift, z
        
        OUTPUTS: D_m(z) in units of h^{-1} Mpc
        
        v1.0 Robert J. Brunner, Jan. 18, 2007
        """
        
        return ((3000.0 / self.h0) * self.iRoutine.integrate(0.0, float(z)))
    
    def da(self, z):
        """Calculates angular-diamter distance
        
        INPUTS: Redshift, z
        
        OUTPUTS: D_a(z) in units of h^{-1} Mpc
        
        v1.0 Robert J. Brunner, Jan. 18, 2007
        """
        
        return (self.dm(float(z)) / (1.0 + float(z)))
        
    def dl(self, z):
        """Calculates luminosity distance
        
        INPUTS: Redshift, z
        
        OUTPUTS: D_l(z) in units of h^{-1} Mpc
        
        v1.0 Robert J. Brunner, Jan. 18, 2007
        """
        
        return (self.dm(float(z)) * (1.0 + float(z)))
        
# Test Code

if __name__ == '__main__':

    c = cdistance()

    print("Luminosity distance at a redshift of 1.0 = ", c.dl(4.17), " h^{-1} Mpc")
    print("Angular Diameter distance at a redshift of 1.0 = ", c.da(1.0), "  h^{-1} Mpc")

    print("Luminosity distance at a redshift of .843 = ", c.dl(0.843), " h^{-1} Mpc", "RI gives 3748.841008")
