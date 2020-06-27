#! /usr/bin/env python3
import numpy as np
import EANN
EANN.init_pes()
coor=np.zeros((3,EANN.nnmod.numatom),dtype=np.float64,order="F")
force=np.zeros((EANN.nnmod.numforce),dtype=np.float64,order="F")
y=np.zeros(1,dtype=np.float64)
table=0
start_force=1
with open('1','r') as f1:
   while True:
      string=f1.readline()
      if not string: break
      for i in range(EANN.nnmod.numatom):
         string=f1.readline()
         m=list(map(float,string.split()))
         for k in range(3):
            coor[k][i]=m[k]
      EANN.eann_out(table,start_force,coor,y,force)
      print(y)
      print(force)
EANN.deallocate_all()
