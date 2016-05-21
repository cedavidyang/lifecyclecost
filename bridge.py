import numpy as np
import scipy.integrate as integrate

from section import Price, Cost, Material, Geometry
from beam import Beam, getrcbeam, getfrpbeam
from column import Column

import sys

class Bridge(object):
    def __init__(self, beam, beamno, column, colno, deck,
            span=None, width=None, height=None, distance=None):
        self.beam = beam
        self.column = column
        self.beamno = beamno
        self.colno = colno
        self.deck = deck
        self.span = span
        self.width = width
        self.height = height
        self.distance = distance

    def setcost(self, cost):
        self.cost = cost

    def setdistance(self, distance):
        seld.distance = distance

    def matcost(self):
        beamno = self.beamno
        colno = self.colno
        span = self.span
        height = self.height
        beam = self.beam
        column = self.column
        cost = beamno*span*beam.matcost() + colno*height*column.matcost()
        return cost

    def transcost(self):
        beamno = self.beamno
        colno = self.colno
        span = self.span
        height = self.height
        beam = self.beam
        column = self.column
        distance = self.distance
        cost = beamno*span*beam.transcost(distance) + colno*height*column.transcost(distance)
        return cost

    def mntcost(self, mntplan):
        beamno = self.beamno
        colno = self.colno
        span = self.span
        width = self.width
        height = self.height
        beam = self.beam
        column = self.column

        bplan = mntplan['beam']
        dplan = mntplan['deck']
        cplan = mntplan['column']

        # direct cost
        bcostpm = self.beam.directmntcost(bplan)
        dcostpm = self.deck.directmntcost(dplan)
        ccostpm = self.column.directmntcost(cplan)
        directcost = beamno*span*bcostpm + width*span*dcostpm + colno*height*ccostpm

        # indirect cost
        pindirect = self.beam.cost.price.mntprice['indirect']
        v = self.beam.cost.price.mntprice['discount']
        allmt = np.hstack((bplan['Mt'],dplan['Mt'],cplan['Mt']))
        allmt = np.unique(allmt)
        nmnt = allmt.size
        # price data
        Cruncar = pindirect['Cruncar']
        Cruntruck = pindirect['Cruntruck']
        Cwage = pindirect['Cwage']
        Cdriver = pindirect['Cdriver']
        Ccargo = pindirect['Ccargo']
        Ttruck = pindirect['Ttruck']
        Adt = pindirect['ADT']
        Ldetour = pindirect['Ldetour']
        Tdetour = pindirect['Tdetour']
        speed = pindirect['S']
        Ocar = pindirect['Ocar']
        Otruck = pindirect['Otruck']

        # indirect cost
        indirectcost = 0.
        for t in allmt:
            Crun = (Cruncar*(1-Ttruck)+Cruntruck*Ttruck)*Ldetour*Tdetour*Adt/(1+v)**t
            Ctime = (Cwage*Ocar*(1-Ttruck)+(Cdriver*Otruck+Ccargo)*Ttruck)*Ldetour*Tdetour*Adt/speed/(1+v)**t
            indirectcost += (Crun+Ctime)

        return (directcost, indirectcost)


if __name__ == '__main__':
    import copy
    # RC bridge
    # price information
    matprice = {'Cconc': 104.57*8.72, 'Csteel': 1.16e3*8.72,
            'Cfb': 3.90e3*7.77, 'Cft': 3.90e3*7.77}
    price = Price(matprice=matprice)
    rcbeam = getrcbeam(price)
    rcdeck = copy.deepcopy(rcbeam)
    rccolumn = copy.deepcopy(rcbeam)
    rcbridge = Bridge(rcbeam, 10, rccolumn, 0, rcdeck,
            span=9.1, width=11.6, height=6., distance=10.)

    # initial cost
    rcmatcost = rcbridge.matcost()
    rctranscost = rcbridge.transcost()
    print 'RC material cost = {}'.format(rcmatcost)
    print 'RC transportation cost = {}'.format(rctranscost)

    # maintenance cost
    rcbridge.beam.geo.addprop({'Arepair':400*1000.,'cover':40.})
    rcbridge.column.geo.addprop({'Arepair':400*1000.,'cover':40.})
    Arepair = 1000.*1000.; cover=30.
    rcbridge.deck.geo.addprop({'Arepair':Arepair,'cover':cover})

    discount = 0.02
    mntpricedirect = {'Cwb': 11.5*9, 'Crm': 43.28*9, 'Cbc': 16.41*9, 'Ctp': 53.71*9}
    mntpriceindirect = {'Cruncar': 0.07*9, 'Cruntruck': 0.34*9, 'Cwage': 20.77*9, 'Cdriver': 24.54*9,
            'Ccargo': 3.64*9, 'Ttruck': 0.12, 'ADT': 8500., 'Ldetour': 2.9, 'Tdetour': 7, 'S': 50,
            'Ocar': 1.5, 'Otruck': 1.05}
    mntprice = {'direct':mntpricedirect, 'indirect':mntpriceindirect, 'discount':discount}

    rcbridge.beam.cost.price.setmntprice(mntprice)
    rcbridge.column.cost.price.setmntprice(mntprice)
    rcbridge.deck.cost.price.setmntprice(mntprice)

    mntplan = {'Mt': np.arange(20., 101., 30.), 'Me': 20.}
    bridgemntplan = {'beam':mntplan, 'column':mntplan, 'deck':mntplan}
    mntplanN = {'Mt': None, 'Me': 10.}
    rcRarray = rcbridge.beam.lifecycleR(75, mntplan)
    rcRarrayN = rcbridge.beam.lifecycleR(75, mntplanN)
    rcmntcost = rcbridge.mntcost(bridgemntplan)
    print 'RC maintenace cost = {}'.format(rcmntcost)

    # # FRP bridge
    # # price information
    # matprice = {'Cconc': 104.57*8.72-100, 'Csteel': 1.16e3*8.72,
            # 'Cfb': 3.90e3*7.77, 'Cft': 3.90e3*7.77}
    # price = Price(matprice=matprice)
    # frpbeam = getfrpbeam(price)
    # # initial cost
    # frpmatcost = frpbeam.matcost()
    # frptranscost = frpbeam.transcost(10.)
    # print 'FRP material cost = {}'.format(frpmatcost)
    # print 'FRP transportation cost = {}'.format(frptranscost)
    # # maintenance cost
    # mntplanN = {'Mt': None, 'Me': 10.}
    # frpRarray = frpbeam.lifecycleR(75, mntplanN)
