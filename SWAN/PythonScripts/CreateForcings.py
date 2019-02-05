# -*- coding: utf-8 -*-
"""
Created on Mon May 14 15:17:06 2018

@author: Nate VanArendonk
"""

speed = 30
direction = range(0,350,10)
tide = [-2,5.5]

forcing_list = []

## Create a set of lists inside of a list for each forcing set
for ii in range(0,len(tide)):
    for jj in range(0,len(direction)):
        forcing_list.append([30,direction[jj],tide[ii]])
        
        



