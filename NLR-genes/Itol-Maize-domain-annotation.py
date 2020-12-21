#!/usr/bin/env python
# coding=utf-8
#Usage: python Itol-Maize-Domain-Annotation.py <NLR-ID.txt file from Plant_R_Genes > > <output_file>
#NB. Any domains not in list below will not be annotated
import sys
import csv
import re

A='(.*)ABC_trans_N(.*)'
B='(.*)BUD22(.*)'
D='(.*)DUF212(.*)'
F='(.*)F-box(.*)'
NAM='(.*)NAM-associated(.*)'
Nop='(.*)Nop14(.*)'
PA='(.*)PAH(.*)'
P2='(.*)PP2\((.*)'
P2C='(.*)PP2C(.*)'
PK='(.*)Pkinase(.*)'
PKT='(.*)Pkinase_Tyr(.*)'
S='(.*)SDA1(.*)'
TF='(.*)ThiF(.*)'
TX='(.*)Thioredoxin(.*)'
TA='(.*)Transpos_assoc(.*)'
U='(.*)UBA_e1_thiolCys(.*)'
UP='(.*)UPRTase(.*)'
WR='(.*)WRKY(.*)'
Z='(.*)zf-RVT(.*)'

"""
1.ABC_trans_N
2.BUD22
3.NAM
4.PAH
5.PKinase
6.PP2
7.Pkinase_Tyr
8.SDA1
9.Thioredoxin
10.WRKY
11.zf-RVT
"""


with open(sys.argv[1]) as fi:
    for row in csv.reader(fi, delimiter=','):
        if re.match(A, row[2]):
            print(row[0]+',#829191,ABC_trans_N')
        if re.match(B, row[2]):
            print(row[0]+',#F2B7C6,BUD22')
        if re.match(D, row[2]):
            print(row[0]+',#694F5D,DUF212')
        if re.match(F, row[2]):
            print(row[0]+',#69656A,F-BOX')
        if re.match(NAM, row[2]):
            print(row[0]+',#744253,NAM')
        if re.match(Nop, row[2]):
            print(row[0]+',#000000,Nop')
        if re.match(PA, row[2]):
            print(row[0]+',#D2BF55,PAH')
        if re.match(P2, row[2]):
            print(row[0]+',#E8DB7D,PP2')
        if re.match(P2C, row[2]):
            print(row[0]+',#000000,PP2C') 
        if re.match(PK, row[2]):
            print(row[0]+',#EDB458,PKinase')       
        if re.match(PKT, row[2]):
            print(row[0]+',#F4F1BB,Pkinase_Tyr')
        if re.match(S, row[2]):
            print(row[0]+',#B3EFB2,SDA1')
        if re.match(TF, row[2]):
            print(row[0]+',#020402,ThiF')
        if re.match(TX, row[2]):
            print(row[0]+',#D9F8C1,Thioredoxin')
        if re.match(U, row[2]):
            print(row[0]+',#000000,UBA')
        if re.match(TA, row[2]):
            print(row[0]+',#000000,Transp')
        if re.match(UP, row[2]):
            print(row[0]+',#020402,UPRTase')
        if re.match(WR, row[2]):
            print(row[0]+',#C5D86D,WRKY')
        if re.match(Z, row[2]):
            print(row[0]+',#657153,zf-RVT')
"""
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
            print(line.rstrip()+',range,#787878,Sobic')

#if re.match(, line):
            #print(line)
"""
