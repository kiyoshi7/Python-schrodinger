'''
Created on Feb 7, 2019

@author: Daniel
'''
from pylab import *
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

class TimeIndependent(): 
    
    ## dimentions         _
    ## |       :       |   | Vmax
    ## |    ___:___    |   |        _
    ## |___|       |___|  _| Vmin   _|- Vb
    ##  -B  -A 0 A   B*
    ## \______/|\______/
    ##     Po       Pf

    # where to start shooting
    option = 1
    Range = 2
    ## values for step Function
    Vmax = 20
    Vmin = 0
    Vb = 20
    SA = 1
    SB = 2
    SPo = -SA-SB
    SPf = SA+SB
    
    ## Values for Quadratic Function
    QB1 = 7.64774E-17
    QB2 = -0.78704
    QB3 = -2.08692E-17
    QB4 = 0.02572
    QB0 = 5
    
    NSteps = 1000
    Xaxis = []
    PotentialFunction = []


    ## Shrodinger
    psi = np.zeros([NSteps,2])     # Wave function values and its derivative (psi and psi')
    psi0 = array([1,0])   # Wave function initial states Asimetrical
    E = 0
    
    def __init__(self):
        #self.Vmax = 10
        pass
    
    ## caculate the potential STEP Function saves it as a list and returns a list
    def V(self, x):
        if self.option == 1:
            if -self.SA <= x <= self.SA:
                return self.Vb
            elif x <= self.SPo:
                return self.Vmax
            elif x >= self.SPf:
                return self.Vmax
            else:
                return self.Vmin
        if self.option == 2:
            return ((self.QB1*x**1)+(self.QB2*x**2)+(self.QB3*x**3)+(self.QB4*x**4)+self.QB0)
        else:
            return(x)
            
    def PotentialStep(self):
        self.SPo = -self.SA-self.SB
        self.SPf = self.SA+self.SB
        self.Xaxis = linspace(-self.Range, self.Range, self.NSteps)
        self.option = 1
        Potential=[]
        X = []
        for i in self.Xaxis:
            Potential.append(self.V(i))
            X.append(i)
        self.PotentialFunction = Potential
        return ([X, Potential])
       
    ## caculate the potential Quadratic Function saves it as a list and returns a list
    def QuadraticFunction(self):
        self.option = 2
        self.Xaxis = linspace(-self.Range, self.Range, self.NSteps)
        Potential=[]
        for i in self.Xaxis:
            Potential.append(self.V(i))
        self.PotentialFunction = Potential
        return (self.Xaxis, Potential)


    ###################
    ##  Schrodinger  ##
    ###################
    
    def SchrodingerEquation(self, x, p):
        """
        Returns derivatives for the 1D schrodinger eq.
        Requires global value E to be set somewhere. State0 is
        first derivative of the wave function psi, and state1 is
        its second derivative.
        """
        state0 = p[1]
        state1 = 2.0*(self.V(x) - self.E)*p[0]
        return array([state0, state1])

    def Wave_function(self, energy):
        """
        Calculates wave function psi for the given value
        of energy E and returns value at point b
        """
        self.E = energy
        psi = solve_ivp(self.SchrodingerEquation,
                        [self.Xaxis[0], self.Xaxis[-1]],
                        self.psi0,
                        max_step = 1*(self.Xaxis[-1]-self.Xaxis[0])/self.NSteps,
                        method = 'LSODA'
                        )
        #psi = odeint( self.SchrodingerEquation, self.psi0, self.Xaxis)
        return psi.t , psi.y
    
    def wavefunctionZero(self, x):
        a = self.Wave_function(x)
        return a[1][0][-1]
    #
    def calculate(self):
        ## create energy by psi[-1]
        y=[]
        x=[]
        e = []
        print(self.Vmin," ",self.Vmax)
        ## generate list of the last item in list of caculated psi, shot at different energies
        for i in np.arange(self.Vmin, self.Vmax, 0.1):                                        ## change here to set arb. energy values
            a=self.wavefunctionZero(i)
            y.append(a)
            x.append(i)
        ## Find where the list above crosses zero and generate the acepted energies
        last = True
        current = True
        for i in range(len(y)):
            if i != 0:
                if i <= 0:
                    current = False
                if i > 0:
                    current = True

                if y[i-1]*y[i] <= 0:
                    ## find zero
                    xe = brentq(self.wavefunctionZero, x[i], x[i-1])
                    ye = self.Wave_function(xe)[1][0][-1]
                    e.append(xe) ## accepted energy
                last = current
        ## Plot where psi crosses 0
        #plt.plot(x,y)
        #plt.show()

        PlotX = []
        PlotY = []
        for i in e:
            a=m.Wave_function(i)
            y=[]
            ShotEnergy = []
            Peak = max(map(abs, a[1][0]))
            for j in a[1][0]:
                y.append(i+j/Peak)
                ShotEnergy.append(i)
            PlotX.append(a[0])
            PlotY.append(ShotEnergy)
            #plt.plot(list(a[0]), list(y))
            #plt.plot(list(a[0]), ShotEnergy, color='black')
        return(PlotX, PlotY)

    def find_analytic_energies(self):
        """
        Calculates Energy values for the finite square well using analytical
        model (Griffiths, Introduction to Quantum Mechanics, 1st edition, page 62.)
        """
        k=0
        if self.Vmin == 0:
            k = 0.01
        en = linspace(self.Vmin+k, self.Vmax, 100)
        z = sqrt(2*en)
        z0 = sqrt(2*self.Vb)
        
        z_zeroes = []
        f_sym = lambda z: tan(z)-sqrt((z0/z)**2-1)      # Formula 2.138, symmetrical case
        print(f_sym(z))
        f_asym = lambda z: -1/tan(z)-sqrt((z0/z)**2-1)  # Formula 2.138, antisymmetrical case
        found = []
        # first find the zeroes for the symmetrical case
        s = sign(f_sym(z))
        for i in range(len(s)-1):   # find zeroes of this crazy function
           if s[i]+s[i+1] == 0:
               zero = brentq(f_sym, z[i], z[i+1])
               z_zeroes.append(zero)
        for i in range(0, len(z_zeroes),2):   # discard z=(2n-1)pi/2 solutions cause that's where tan(z) is discontinous
            found.append(z_zeroes[i]**2/2)

        
        # Now for the asymmetrical
        z_zeroes = []
        s = sign(f_asym(z))
        for i in range(len(s)-1):   # find zeroes of this crazy function
           if s[i]+s[i+1] == 0:
               zero = brentq(f_asym, z[i], z[i+1])
               z_zeroes.append(zero)
        for i in range(0, len(z_zeroes),2):   # discard z=npi solutions cause that's where ctg(z) is discontinous
            found.append(z_zeroes[i]**2/2)
        return found

    


m=TimeIndependent()
m.QuadraticFunction()
#a=m.calculate()
#a=m.find_analytic_energies()
