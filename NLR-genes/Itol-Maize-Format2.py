#!/usr/bin/env python
# coding=utf-8

#Usage: python Itol-Maize-Format2.py <list-of-tree-Maize-IDs> > <output_file>
 
import sys
import csv
import re
b73='(.*)1eb(.*)'
b97='(.*)18ab(.*)'
cml103='(.*)21ab(.*)'
cml228='(.*)22ab(.*)'
cml247='(.*)23ab(.*)'
cml277='(.*)24ab(.*)'
cml322='(.*)25ab(.*)'
cml333='(.*)26ab(.*)'
cml52='(.*)19ab(.*)'
cml69='(.*)20ab(.*)'

hp3='(.*)27ab(.*)'
il14h='(.*)28ab(.*)'
ki11='(.*)30ab(.*)'
ki3='(.*)29ab(.*)'
ky21='(.*)31ab(.*)'
m162w='(.*)33ab(.*)'
m37w='(.*)32ab(.*)'
mo18w='(.*)34ab(.*)'
ms71='(.*)35ab(.*)'
nc350='(.*)36ab(.*)'
nc358='(.*)37ab(.*)'
oh43='(.*)39ab(.*)'
oh7b='(.*)38ab(.*)'
p39='(.*)40ab(.*)'
tx303='(.*)41ab(.*)'
tzi8='(.*)42ab(.*)'
Sobic='(.*)Sobic(.*)'

with open(sys.argv[1]) as fi:

    for line in fi:
        if re.match(b73, line):
            print(line.rstrip()+',range,#FFC125,b73')
        if re.match(b97, line):
            print(line.rstrip()+',range,#4169E1,b97')
        if re.match(ky21, line):
            print(line.rstrip()+',range,#4169E1,ky21')
        if re.match(m162w, line):
            print(line.rstrip()+',range,#4169E1,m162w')
        if re.match(oh43, line):
            print(line.rstrip()+',range,#4169E1,oh43')
        if re.match(oh7b, line):
            print(line.rstrip()+',range,#4169E1,oh7b')
        if re.match(ms71, line):
            print(line.rstrip()+',range,#4169E1,ms71')
        
        if re.match(cml103, line):
            print(line.rstrip()+',range,#32CD32,cml103')
        if re.match(cml228, line):
            print(line.rstrip()+',range,#32CD32,cml228')
        if re.match(cml247, line):
            print(line.rstrip()+',range,#32CD32,cml247')
        if re.match(cml277, line):
            print(line.rstrip()+',range,#32CD32,cml277')
        if re.match(cml322, line):
            print(line.rstrip()+',range,#32CD32,cml322')
        if re.match(cml333, line):
            print(line.rstrip()+',range,#32CD32,cml333')
        if re.match(cml52, line):
            print(line.rstrip()+',range,#32CD32,cml52')
        if re.match(cml69, line):
            print(line.rstrip()+',range,#32CD32,cml69')
        if re.match(ki11, line):
            print(line.rstrip()+',range,#32CD32,ki11')
        if re.match(ki3, line):
            print(line.rstrip()+',range,#32CD32,ki3')
        if re.match(nc350, line):
            print(line.rstrip()+',range,#32CD32,nc350')
        if re.match(nc358, line):
            print(line.rstrip()+',range,#32CD32,nc358')
        if re.match(tzi8, line):
            print(line.rstrip()+',range,#32CD32,tzi8')


        if re.match(hp3, line):
            print(line.rstrip()+',range,#DA70D6,hp301')

        if re.match(il14h, line):
            print(line.rstrip()+',range,#FF4500,il14h')
        if re.match(p39, line):
            print(line.rstrip()+',range,#FF4500,p39')

        if re.match(m37w, line):
            print(line.rstrip()+',range,#787878,m37w')
        if re.match(mo18w, line):
            print(line.rstrip()+',range,#787878,mo18w')
        if re.match(tx303, line):
            print(line.rstrip()+',range,#787878,tx303')


        if re.match(Sobic, line):
            print(line.rstrip()+',range,#94714f,Sobic')

#if re.match(, line):
            #print(line)
