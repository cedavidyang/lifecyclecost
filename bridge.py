import numpy as np
import scipy.integrate as integrate

from section import Price, Cost, Material, Geometry
from beam import Beam
from column import Column

import sys

class Bridge(object):
    def __init__(self, beam, beamno, column, colno, span, height, distance):
        self.beam = beam
        self.column = column
        self.beamno = beamno
        self.colno = colno
        self.span = span
        self.height = height
        self.distance = distance

    def matcost(self):
        beamno = self.beamno
        colno = self.colno
        span = self.span
        heigh = self.height
        beam = self.beam
        column = self.column
        cost = beamno*span*beam.matcost() + colno*height*column.matcost()
        return cost

    def transcost(self):
        beamno = self.beamno
        colno = self.colno
        span = self.span
        heigh = self.height
        beam = self.beam
        column = self.column
        distance = self.distance
        cost = beamno*span*beam.transcost(distance) + colno*height*column.transcost(distance)
        return cost


if __name__ == '__main__':
    # price information
    price = Price(1.0, 1.0, 1.0, 1.0)
    # # initial cost
    # initalcost = bridge.matcost()+bridge.transcost()
