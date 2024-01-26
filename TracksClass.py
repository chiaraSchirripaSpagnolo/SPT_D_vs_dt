# -*- coding: utf-8 -*-
"""
Last Modified 25 Jan 2024

@author: Chiara
"""

class Tracks:
    def __init__ (self, idxInMovie, coord, stoich, events ):
        self.idxInMovie=idxInMovie #list of lists: one list for each subtrajectory
        self.coord=coord #list of lists of lists: one list for each subtrajectory
                        #in each one 2 lists: one for x, one for y, for each frame
        self.stoich=stoich #list of lists of arrays: one list for each subtrajectory
                            #in each one a sequence of arrays of stoichiometries for each frame
        self.events=events #list of lists: each list like one row in seqOfEvents