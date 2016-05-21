import numpy as np
import scipy.integrate as integrate
from scipy.optimize import minimize

from section import Price, Cost, Material, Geometry
from beam import Beam
import sys

class Column(Beam):
    def getcompression(self):
        if self.geo.rtype.lower == 'rc':
            # consider reinforcement contribution
            Nr = np.sum(self.geo.Ar*self.mat.fr)
            Nc = self.geo.Ac*self.mat.fc
        else:
            Nr = 0.0
            Nc = self.geo.Ac*self.mat.fc
        return Nr+Nc

    def getcapacity(self, N):
        """ nh: number of discretiztions along beam (effective) depth d"""
        h = self.geo.h
        xr = self.geo.xr
        Ar = self.geo.Ar
        d = np.max(xr)
        bdist = self.geo.bdist    # distribution of concrete width
        ecu = self.mat.ecu
        eru = self.mat.eru
        if self.geo.rtype.lower() == 'rc':
            negeru = -eru
        else:
            negeru = -eru
        concss=self.mat.concss    #concrete stress-strain relation
        reinss = self.mat.reinss    # reinforcement stress-strain relation

        rfmin = 1.e16
        def objfunc(c):
            # assume concrete crush, determine steel/frpbar strain
            er = (d-c)/c*ecu
            if er>eru or er<negeru:    # control by reinforcement failure
                ec = c/(d-c)*eru
            else:
                ec = ecu
            ecdist = lambda x: -ec/c*x+ec
            erdist = lambda x: (x-c)/c*ec
            Fc,err = integrate.quad(lambda x: concss(ecdist(x))*bdist(x), 0, c)
            Fr = np.sum(reinss(erdist(xr))*Ar)
            rf = np.abs(Fc-Fr-N)
            return rf
        c0 = 0.5*d
        res = minimize(objfunc, c0, bounds=((0.1*d,None),))
        csol = res.x

        er = (d-csol)/csol*ecu
        if er>eru or er<negeru:    # control by reinforcement failure
            ecsol = csol/(d-csol)*eru
        else:
            ecsol = ecu
        ecdistsol = lambda x: -ecsol/csol*x+ecsol
        erdistsol = lambda x: (x-csol)/csol*ecsol

        # calculate the capacity
        Fc,err = integrate.quad(lambda x: concss(ecdistsol(x))*bdist(x), 0, csol)
        Fr = np.sum(reinss(erdistsol(xr))*Ar)
        Ncomp = np.abs(Fc-Fr)
        Mc,err = integrate.quad(lambda x: concss(ecdistsol(x))*bdist(x)*x, 0, csol)
        Mr = np.sum(reinss(erdistsol(xr))*Ar*xr)
        M = np.abs(Mc-Mr-Ncomp*h/2.)
        return Ncomp,M,csol


if __name__ == '__main__':
    price = Price(1.0, 1.0, 1.0, 1.0)
    cost = Cost(price)
    ## RC example ex7.5 7.6 7.7 (Shu Shilin concrete textbook) checked
    # define material
    fc = 9.6
    fr = 300; Er = 200e3
    ecu = 0.0033
    eco = 0.002
    eru = fr/Er*100
    def concss(ec):    # must be scalar function
        if ec<=eco:
            sc = fc*(1.-(1.-ec/eco)**(2-20./60))
        else:
            sc = fc
        return sc
    # reinss must be a vector function
    reinss = lambda x: np.maximum(np.minimum(x*Er, fr), -fr)
    mat = Material(fc, fr, eco, ecu, eru, concss, reinss)
    # define geometry
    h = 500.
    Ac = 300.*500.
    Afb = np.array([0.])
    Aft = 0.
    As = np.array([1012., 308.])
    xf = np.array([0.])
    xs = np.array([40.,460.])
    def bdist(x):
        if x<0 or x>500:
            b = 0.
        else:
            b = 300.
        return b
    geo = Geometry('rc', h, Ac, Afb, Aft, As, xf, xs, bdist)
    # define beam
    column = Column(geo=geo, mat=mat, cost=cost)
    Nmax = column.getcompression()
    print 'axial compression = {} kN'.format(Nmax/1e3)
    N,M,csol = column.getcapacity(1150e3)
    print 'N = {} kN; M = {} kN-m'.format(N/1e3, M/1e6)
