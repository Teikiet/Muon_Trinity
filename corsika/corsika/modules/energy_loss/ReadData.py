#!/usr/bin/python

import os
import sys
import string


data = open("Properties8.dat", "r")

nEl = 0
nOpt = 0
count = 0 # the data comes with flexible row-blocks, need to count
for line in data.readlines():
    line = line.rstrip()
    if line == "" or line[0] == "#":
        continue

#    print str(count) + " " + line
    
    if count == 0:
        index = int(line[0:4].strip()) # sternheimers index
        name = line[5:11].strip() # tag
        sig = int(line[11:14].strip()) # significant figures of atomic mass if element
        weight = float(line[16:27].strip()) # atomic weight, if element
        weight_error = int((line[31:32].strip() or "0"))# error in last place
        ZoverA = float(line[34:42].strip()) # Z/A
        rho = float(line[42:52].strip())          # Sternheimers density
        rho_corr = float(line[52:63].strip())# Corrected density
        state = line[64].strip() # Solid Liquid Gas Diatomicgas
        nEl = int(line[66:69].strip()) # number of elements
        nAtom = int(line[70:72].strip()) # atoms of el. 1
        nOpt = int(line[73:74].strip()) # number of optional lines
        type = line[75] # type: Element, Radiactive element, Inorganic compound, Organic compound, Polymer, Mixture, Biological

    elif count == 1:
        short = line[0:4].strip() # if element
        long = line[4:100].strip()
                   
    elif count == 2:
        desc = line # description and formula
                   
    elif count == 3:
        Ieff = float(line[0:10].strip()) # ieff
        Cbar = float(line[11:19].strip()) # cbar
        x0 = float(line[21:28].strip()) # x0
        x1 = float(line[30:37].strip()) # x1
        aa = float(line[39:46].strip()) # a / aa
        sk = float(line[48:55].strip()) # k / sk
        delta0 = float(line[59:].strip()) # delta0  / dlt0
                   
    elif count == 4 and count < 4+nEl:
        elZ = int(line[0:11].strip()) # Z
        # num frac.
        # weigh frac.
        pass
                   
    elif count>=4+nEl and count<=4+nEl+nOpt:
        # print "opt : " + str(count) + " " + str(4+nEl+nOpt) + " " + line
        # property A25
        # value F20
        # if 1:5 is "Note:" following lines are extra comments
        if (count == 4+nEl+nOpt):
            count = -1

            print str(elZ) + " " + long + " " + str(Ieff) + " " + str(ZoverA) + " " + state + " " + type + " " + str(rho_corr) + " " + str(delta0) + " " + str(x0)  + " " + str(x1)  + " " + str(aa)  + " " + str(sk) + " " + str(Cbar)

            
    count += 1

                  
if (count != 0):
    print ("Error " + str(count))
                
