#! /usr/bin/env python3
import numpy as np
import PEANN
PEANN.init_pes()
coor=np.zeros((3,PEANN.sharedmod.numatom),dtype=np.float64,order="F")
force=np.zeros((PEANN.sharedmod.numforce),dtype=np.float64)
y=np.zeros(1,dtype=np.float64)
table=0
start_force=1
with open('1','r') as f1:
   while True:
      string=f1.readline()
      if not string: break
      for i in range(PEANN.sharedmod.numatom):
         string=f1.readline()
         m=list(map(float,string.split()))
         for k in range(3):
            coor[k][i]=m[k]
      PEANN.peann_out(table,start_force,coor,y,force)
      print(y)
      print(force)
PEANN.deallocate_all()
