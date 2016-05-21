import numpy as np
import scipy.integrate as integrate
import sys

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
        self.conccost = conccost
        self.steelcost = steelcost
        self.fbcost = fbcost
        self.ftcost = ftcost
        self.matcost = conccost+steelcost+fbcost+ftcost

    def settranscost(self, conccost, steelcost, fbcost, ftcost):
        self.conccost = conccost
        self.steelcost = steelcost
        self.fbcost = fbcost
        self.ftcost = ftcost
        self.transcost = conccost+steelcost+fbcost+ftcost


class Price(object):
    def __init__(self, conc, steel, fb, ft):
        self.conc = conc
        self.steel = steel
        self.fb = fb
        self.ft = ft

    def gettranscost(self, cargo, distance):
        if cargo.lower() == 'concrete':
            cost = 1
        elif cargo.lower() == 'steel':
            cost = 2
        elif cargo.lower() == 'fb':
            cost = 3
        elif cargo.lower() == 'ft':
            cost = 4
        else:
            print "unknow type of cargo"
            sys.exit(1)
        return cost

    def setmntcost(self, directcost, indirectcost, discount):
        self.directcost = directcost
        self.indirectcost = indirectcost
        self.discount = discount
