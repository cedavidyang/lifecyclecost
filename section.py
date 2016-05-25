import numpy as np
import scipy.integrate as integrate
import sys

RHO_STEEL = 7800.    #kg/m3
RHO_FRP = 1800.
RHO_CONC = 2500.

class Geometry(object):
    def __init__(self, rtype, h, Ac, Afb, Aft, As, xf, xs, bdist):
        self.rtype=rtype
        self.h = h
        self.Ac = Ac
        self.Afb = Afb
        self.Aft = Aft
        self.As = As
        self.xf = xf
        self.xs = xs
        self.bdist = bdist
        if rtype.lower() == 'frp':
            self.Ar = Afb
            self.xr = xf
        elif rtype.lower() == 'rc':
            self.Ar = As
            self.xr = xs
        else:
            print "unknown type, must by frp or rc"
            sys.exit(1)

    def addprop(self,propdict):
        try:
            self.ds = propdict['ds']
        except KeyError:
            pass
        try:
            self.Arepair = propdict['Arepair']
        except KeyError:
            pass
        try:
            self.cover = propdict['cover']
        except KeyError:
            pass


class Material(object):
    def __init__(self, fc, fr, eco, ecu, eru, concss, reinss):
        self.fc = fc
        self.fr = fr
        self.eco = eco
        self.ecu=ecu
        self.eru=eru
        self.concss = concss
        self.reinss = reinss


class Cost(object):
    def __init__(self, price):
        self.price = price

    def setmatcost(self, conccost, steelcost, fbcost, ftcost):
        self.matcost = conccost+steelcost+fbcost+ftcost

    def settranscost(self, conccost, steelcost, fbcost, ftcost):
        self.transcost = conccost+steelcost+fbcost+ftcost


class Price(object):
    def __init__(self, matprice=None, mntprice=None, failprice=None):
        self.matprice = matprice
        self.mntprice = mntprice
        self.failprice = failprice

    def setmntprice(self, mntprice):
        self.mntprice = mntprice

    def setmatprice(self, matprice):
        self.matprice = matprice

    def setfailprice(self, failprice):
        self.failprice = failprice

    def transprice(self, cargo, distance):
        if cargo.lower() == 'sand':
            transprice = 0.035*8*distance
        elif cargo.lower() == 'water':
            transprice = 0.035*8*distance
        elif cargo.lower() == 'concrete':
            transprice = 0.035*8*distance
        elif cargo.lower() == 'steel':
            transprice = 0.035*8*distance
        elif cargo.lower() == 'fb':
            transprice = 0.035*8*distance
        elif cargo.lower() == 'ft':
            transprice = 0.035*8*distance
        else:
            print "unknow type of cargo"
            sys.exit(1)
        return transprice


def frpdegrade(t):
    """ t in week"""
    r = (1-0.8)*np.exp(-t/216.1)+0.8
    return r
