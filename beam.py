import numpy as np
import scipy.integrate as integrate

from section import Price, Cost, Material, Geometry
import sys

class Beam(object):
    def __init__(self, geo=None, mat=None, cost=None):
        self.geo=geo
        self.mat=mat
        self.cost = cost

    def setgeo(self, geo):
        self.geo=geo

    def setmat(self, mat):
        self.mat=mat

    def setcost(self, cost):
        self.cost=cost

    def matcost(self):
        if self.geo is None or self.mat is None:
            print "geometry and material are not defined"
            sys.exit(1)
        if self.cost.price is None:
            print "Cost is not initiated properly"
            sys.exit(1)
        conccost = self.geo.Ac*self.cost.price.conc
        steelcost = self.geo.As*self.cost.price.steel
        fbcost = self.geo.Afb*self.cost.price.fb
        ftcost = self.geo.Aft*self.cost.price.ft
        self.cost.setmatcost(conccost, steelcost, fbcost, ftcost)
        return self.cost.matcost

    def transcost(self, distance):
        if self.geo is None or self.mat is None:
            print "geometry and material are not defined"
            sys.exit(1)
        if self.cost.price is None:
            print "Cost is not initiated properly"
            sys.exit(1)
        conccost = self.geo.Ac*self.cost.price.transcost('concrete', distance)
        steelcost = self.geo.As*self.cost.price.transcost('steel', distance)
        fbcost = self.geo.Afb*self.cost.price.transcost('fb', distance)
        ftcost = self.geo.Aft*self.cost.price.transcost('ft', distance)
        self.cost.settranscost(conccost, steelcost, fbcost, ftcost)
        return self.cost.matcost

    def getmoment(self, nd=100, carray=None):
        """ nh: number of discretiztions along beam (effective) depth d
            carray: given carray
        """
        xr = self.geo.xr
        Ar = self.geo.Ar
        dmin = np.min(xr)
        dmax = np.max(xr)
        bdist = self.geo.bdist    # distribution of concrete width
        ecu = self.mat.ecu
        eru = self.mat.eru
        concss=self.mat.concss    #concrete stress-strain relation
        reinss = self.mat.reinss    # reinforcement stress-strain relation

        if carray is None:
            carray = np.linspace(dmin/nd, dmin, num=nd)
        rfmin = 1.e16
        for c in carray:
            # assume concrete crush, determine steel/frpbar strain
            er = (dmax-c)/c*ecu
            if er>eru:    # control by reinforcement failure
                ec = c/(dmax-c)*eru
            else:
                ec = ecu
            ecdist = lambda x: -ec/c*x+ec
            erdist = lambda x: (x-c)/c*ec
            Fc,err = integrate.quad(lambda x: concss(ecdist(x))*bdist(x), 0, c)
            Fr = np.sum(reinss(erdist(xr))*Ar)
            rf = np.abs(Fc-Fr)
            if rf<rfmin:
                rfmin = rf
                csol = c
                ecsol = ec
                ecdistsol = lambda x: -ecsol/csol*x+ecsol
                erdistsol = lambda x: (x-csol)/csol*ecsol

        # calculate the capacity
        Mc,err = integrate.quad(lambda x: concss(ecdistsol(x))*bdist(x)*x, 0, csol)
        Mr = np.sum(reinss(erdistsol(xr))*Ar*xr)
        M = np.abs(Mc-Mr)
        self.capacity = M
        return M,csol


def getrcbeam(price):
    cost = Cost(price)
    fc = 20.7
    ecu = 0.0038
    fr = 496.; Er = 200000.
    eru = fr/Er*100.
    Ec = 4700*np.sqrt(fc)
    fcp = 0.9*fc
    eco = 1.8*fcp/Ec
    def concss(ec):
        if ec<=eco:
            sc = fcp*(2*ec/eco-(ec/eco)**2)
        else:
            sc = fcp-0.15*fcp*(ec-eco)/(0.0038-eco)
        return sc
    reinss = lambda x: np.minimum(x*Er, fr)
    mat = Material(fc, fr, eco, ecu, eru, concss, reinss)
    # define geometry
    h=790.
    Ac = 400.*(790.-190)+2600.*190.
    Afb = np.array([0.])
    Aft = 0.
    As = np.array([8*645.])
    xf = np.array([0.])
    xs = np.array([687.3+2*(36-29)])
    def bdist(x):
        if x<0 or x>790:
            b = 0.
        elif x>=0 and x<190:
            b = 2600.
        else:
            b = 400.
        return b
    geo = Geometry('rc', h, Ac, Afb, Aft, As, xf, xs, bdist)
    # define beam
    beam = Beam(geo=geo, mat=mat, cost=cost)
    return beam


def getfrpbeam(price):
    cost = Cost(price)
    fc = 20.7
    ecu = 0.0038
    fr = 496.; Er = 44800.
    eru = fr/Er
    Ec = 4700*np.sqrt(fc)
    fcp = 0.9*fc
    eco = 1.8*fcp/Ec
    def concss(ec):
        if ec<=eco:
            sc = fcp*(2*ec/eco-(ec/eco)**2)
        else:
            sc = fcp-0.15*fcp*(ec-eco)/(0.0038-eco)
        return sc
    reinss = lambda x: np.minimum(x*Er, fr)
    mat = Material(fc, fr, eco, ecu, eru, concss, reinss)
    # define geometry
    h=790.
    Ac = 400.*(790.-190)+2600.*190.
    Afb = np.array([8*645.])
    Aft = 0.
    As = np.array([0.])
    xf = np.array([687.3+2*(36-29)])
    xs = np.array([0.])
    def bdist(x):
        if x<0 or x>790:
            b = 0.
        elif x>=0 and x<190:
            b = 2600.
        else:
            b = 400.
        return b
    geo = Geometry('frp', h, Ac, Afb, Aft, As, xf, xs, bdist)
    # define beam
    beam = Beam(geo=geo, mat=mat, cost=cost)
    print 'material cost = {}'.format(beam.matcost())
    M,c = beam.getmoment()
    print 'capacity = {} kN-m; neutral axis at x={}'.format(M/1e6, c)
    return beam


if __name__ == '__main__':
    price = Price(1.0, 1.0, 1.0, 1.0)
    cost = Cost(price)
    ## RC example ex4.1, 4.2, 4.10, 4.11, 4.12 (Shu Shilin concrete textbook) checked
    # define material
    fc = 11.9
    fr = 300; Er = 200e3
    ecu = 0.0033
    eco = 0.002
    eru = fr/Er*100
    def concss(ec):
        if ec<=eco:
            sc = fc*(1.-(1.-ec/eco)**(2-25./50))
        else:
            sc = fc
        return sc
    reinss = lambda x: np.minimum(x*Er, fr)
    mat = Material(fc, fr, eco, ecu, eru, concss, reinss)
    # define geometry
    h=500.
    Ac = 200.*(500.-120)+500.*120.
    Afb = np.array([0.])
    Aft = 0.
    As = np.array([2454.])
    xf = np.array([0.])
    xs = np.array([442.5])
    def bdist(x):
        if x<0 or x>500:
            b = 0.
        elif x>=0 and x<120:
            b = 500.
        else:
            b = 200.
        return b
    geo = Geometry('rc', h, Ac, Afb, Aft, As, xf, xs, bdist)
    # define beam
    beam = Beam(geo=geo, mat=mat, cost=cost)
    print 'material cost = {}'.format(beam.matcost())
    M,c = beam.getmoment()
    print 'capacity = {} kN-m; neutral axis at x={}'.format(M/1e6, c)


    ## FRP example example beam in ACI440.1R-06 checked (concrete stress-strain curve
    ## in Wight and MacGregor Textbook has been used)
    # define material
    fc = 27.6
    fr = 496.; Er = 44800.
    ecu = 0.003
    eru = fr/Er
    Ec = 4700*np.sqrt(fc)
    fcp = 0.9*fc
    eco = 1.8*fcp/Ec
    def concss(ec):
        if ec<=eco:
            sc = fcp*(2*ec/eco-(ec/eco)**2)
        else:
            sc = fcp-0.15*fcp*(ec-eco)/(0.0038-eco)
        return sc
    reinss = lambda x: np.minimum(x*Er, fr)
    mat = Material(fc, fr, eco, ecu, eru, concss, reinss)
    # define geometry
    h=305.
    Ac = 305.*178.
    Afb = np.array([400.])
    Aft = 0.
    As = np.array([0.])
    xf = np.array([250.])
    xs = np.array([0.])
    def bdist(x):
        if x<0 or x>305:
            b = 0.
        else:
            b = 178.
        return b
    geo = Geometry('frp', h, Ac, Afb, Aft, As, xf, xs, bdist)
    # define beam
    beam = Beam(geo=geo, mat=mat, cost=cost)
    print 'material cost = {}'.format(beam.matcost())
    M,c = beam.getmoment()
    print 'capacity = {} kN-m; neutral axis at x={}'.format(M/1e6, c)
