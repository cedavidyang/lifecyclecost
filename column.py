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
        res = minimize(objfunc, c0, method='L-BFGS-B', bounds=((0.1*d,None),))
        # res = minimize(objfunc, c0, bounds=((0.1*d,None),))
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

def getrccolumn(price):
    cost = Cost(price)
    # define material
    fc = 30.0
    fr = 300; Er = 200e3
    ecu = 0.0038
    eru = fr/Er*100
    Ec = 4700*np.sqrt(fc)
    fcp = 0.9*fc
    eco = 1.8*fcp/Ec
    def concss(ec):
        if ec<=eco:
            sc = fcp*(2*ec/eco-(ec/eco)**2)
        else:
            sc = fcp-0.15*fcp*(ec-eco)/(0.0038-eco)
        return sc
    # reinss must be a vector function
    reinss = lambda x: np.maximum(np.minimum(x*Er, fr), -fr)
    mat = Material(fc, fr, eco, ecu, eru, concss, reinss)
    # define geometry
    h = 700.
    Ac = np.pi/4*(h)**2
    Afb = np.array([0.])
    Aft = 0.
    As = np.array([645., 1290, 1290, 1290, 1290, 1290, 645])
    xf = np.array([0.])
    xs = np.array([65., 103, 207.5, 350, 492.5, 597, 635])
    def bdist(x):
        if x<0 or x>h/2.:
            b = 0.
        else:
            b = 2.*np.sqrt((h/2.)**2-(x-(h/2.))**2)
        return b
    geo = Geometry('rc', h, Ac, Afb, Aft, As, xf, xs, bdist)
    # define beam
    column = Column(geo=geo, mat=mat, cost=cost)
    Nmax = column.getcompression()
    print 'axial compression = {} kN'.format(Nmax/1e3)
    N,M,csol = column.getcapacity(6000e3)
    print 'N = {} kN; M = {} kN-m; c={}'.format(N/1e3, M/1e6, csol)
    return column


def getfrpcolumn(price):
    cost = Cost(price)
    # define geometry
    h = 700.
    Ac = np.pi/4*h**2
    Afb = np.array([0.])
    tf = 6.0;
    Aft = np.pi*h*tf
    As = np.array([387, 774.,774.,774.,774.,774.,387.])
    xf = np.array([0.])
    xs = np.array([40., 81.5, 195, 350, 505, 618.5, 660])
    def bdist(x):
        if x<0 or x>h/2.:
            b = 0.
        else:
            b = 2.*np.sqrt((h/2.)**2-(x-(h/2.))**2)
        return b
    geo = Geometry('rc', h, Ac, Afb, Aft, As, xf, xs, bdist)
    # define material
    fc = 30.0
    fr = 500; Er = 45e3
    eco = 0.002
    psif=0.95; ke=0.55; ka=1.0; kb=1.0
    efe = ke*fr/Er
    fl = 2*Er*tf*efe/h
    print 'confinement ratio = {}'.format(fl/fc)
    fcc = fc+psif*3.3*ka*fl
    eccu = eco*(1.50+12*kb*fl/fc*(efe/eco)**0.45)
    print 'fcc = {}'.format(fcc)
    print 'ecc = {}'.format(eccu)
    if eccu>0.01: eccu=0.01
    eru = fr/Er*100
    Ec = 4700*np.sqrt(fc)
    def concss(ec):
        E2 = (fcc-fc)/eccu
        etp = 2*fc/(Ec-E2)
        if ec<=etp:
            sc = Ec*ec-(Ec-E2)**2/(4*fc)*ec**2
        elif ec>etp and ec<eccu:
            sc = fc+E2*ec
        else:
            sc = 0
        return sc
    ecu = eccu
    fc = fcc
    # reinss must be a vector function
    # reinss = lambda x: np.maximum(np.minimum(x*Er, fr), -fr)
    def reinss(x):
        sr = np.maximum(np.minimum(x*Er, fr), -fr)
        if isinstance(x,float) and sr<0:
            sr = 1e-6*sr
        else:
            sr[sr<0] = 1e-6*sr[sr<0]
        return sr
    mat = Material(fc, fr, eco, ecu, eru, concss, reinss)
    # define beam
    column = Column(geo=geo, mat=mat, cost=cost)
    Nmax = column.getcompression()
    print 'axial compression = {} kN'.format(Nmax/1e3)
    N,M,csol = column.getcapacity(7500.e3)
    print 'N = {} kN; M = {} kN-m; c={}'.format(N/1e3, M/1e6, csol)

    return column

if __name__ == '__main__':
    matprice = {'Cconc': 104.57*8.72, 'Csteel': 4871.,
            'Cfb': 3.90e3*7.77, 'Cft': 3.90e3*7.77}
    price = Price(matprice=matprice)
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
