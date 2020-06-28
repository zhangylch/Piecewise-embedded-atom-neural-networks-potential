# Piesewise embedded atom neural network package

Piesewise embedded atom neural network (PEANN) is a Fortran package to fit highly efficient and accurate neural network models for interatomic potentials and molecular properties of complex systems. This project is still under active development. The current version only includes the function of using a trained model. The full package, including automatic generation of dataset, automatic filter of descriptors, as well as training new models, will be released soon.<br>

## 1. Compling a model:
----------------------
### 1. Direct employment in Fortran:
------------------
At the beginning, users need to compile the Fortran files to obtain an executable file by the script "makefile". Then users can call the subroutine ‘"init_pes" to do some initializations. There is a "test.f90" which is an example to calculate the energy and force of a given configuration.Then users can pass the coordinates of configurations by calling the subroutine.<br>
>"peann_out(table, start_force, coor, y, force)".<br>
>>1. "table": The form of the coordinates. "0": Cartesian coordinates. "1": Fraction coordinates; <br>
>>2. "start_force": "1": energies and force vectors output with the variable "y" and "force" in the subroutine respectively; "0": only energies obtained;<br>
>>3. "coor": the coordinates of the configurations provided for the calculating;<br>
>>4. "y": the array holding the energy of the corresponding configurations; <br>
>>5. "force": the array holding the atomic force of the corresponding configurations with its dimension numforce (numforce is the number of atomic forces).<br>
Finally, "deallocate_pes" was employed to release the allocated memories.<br>

### 2. Fortran to Python interface:
-----------------------
These Fortran files can be wrapped into the PEANN.cpython-37m-x86_64-linux-gnu.so by runing the script "run". The interface is similar with the Fortran interface. "test.py" is an example.<br>
```python
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
```
At the beginning, users should import the module "PEANN". Then the function named "init_pes" should be called to do some initializations including load the model parameters,allocate the allocatable arrays, etc. Next, the energy (y) and force vectors (force) of corresponding configuration (coor) can be calculated by calling "peann_out". Arguments in the "peann_out" are same with that in Fortran. Finally, "deallocate_pes" was employed to release the allocated memories.<br>
With a compiled executable file, users need to copy the trained model parameters in the floder "potential" to the folder "para" in the work directory, including "input_density", "input", "cell", "atom" and paramters of neural networks. "input_density" and "input_nn" contain the hyper-parameters of descriptors and neural networks. "cell" contains the lattice parameters of the system. "atom" contains the element types in the order of the input geometry.<br>

## 3. Embedded atom neural network:
-----------------------------------
Using an EANN model is similar as using a PEANN model. The only difference is to use "EANN_out" instead of "PEANN_out".

## 4. Ab initio data source:
---------------------------------------------
### Bulk water:
The bulk water data were calculated by Bingqing Cheng *et.al.* *Proceedings of the National Academy of Sciences* 116, 1110-1115 (2019). Original data can be found in https://archive.materialscloud.org/2018.0020/v1.<br>
### Metal systems:
Cu, Ge, Mo systems were reported by Yunxing Zuo *et.al.* *J. Phys. Chem. A* 124, 731-745 (2020). Original data can be found in https://github.com/materialsvirtuallab/mlearn

## 5. References:
The original EANN model: Yaolong Zhang *et. al.* *J. Phys. Chem. Lett.* 10, 4962-4967 (2019).

The PEANN model: 

The EANN model for electronic friction tensor: Yaolong Zhang *et. al.* *J. Phys. Chem. C* 124, 186 (2020).

The EANN model for dipole/transition dipole/polarizability: arXiv:2004.13605
