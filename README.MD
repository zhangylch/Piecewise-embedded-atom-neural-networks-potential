Piesewise embedded atom neural network (PEANNN) and embedded atom neural network(EANN)
===========
1.Employ the model:
---------------
With a trained model, users can employ the model to obtain the energies and atomic force vectors of configurations. Furthermore, the models can be interfaced into the MD package for atomic simulation. In order to employ the model, you need to copy the trained model parameters in the floder "potential" to the folder "para" in the work directory.  including "input","cell", "atom" and paramters of neural networks. "cell" contains the lattice parameter of the system. "atom" contains the atom names with the order of the input geometry.<br>

2.Fortran interface:
------------------
At the beginning, users need to compile the Fortran files to obtain an executable file by the script "makefile". Then users need to call the subroutine ‘"init_pes" to do some initializations. There is a "test.f90" which is an example to calculate energy and force of corresponding configurations.Then you can pass the coordinates of configurations by calling the subroutine.<br>
>>"peann_out(table, start_force, coor, y, force)".<br>
>>>>1. "table": The form of the coordinates. "0": Cartesian coordinates. "1": Fraction coordinates; <br>
>>>>2. "start_force": "1": energies and force vectors output with the variable ‘"y" and ‘"force" in the subroutine respectively; "0": only energies obtained;<br>
>>>>3. "coor": the coordinates of the configurations provided for the calculating;<br>
>>>>4. "y": the array holding the energy of the corresponding configurations; <br>
>>>>5. "force": the array holding the atomic force of the corresponding configurations with its dimension numforce (numforce is the number of atomic forces).<br>

3.Python interface:
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
At the beginning, users should import the module "PEANN". Then the function named "init_pes" should be called to do some initialization including load the modelparameters,allocate the allocatable arrays, etc. Next, the energy (y) and force vectors (force) of corresponding configuration (coor) can be calculated by calling "peann_out" Variables in the "peann_out" are same with that in Fortran interface. Finally, "deallocate_pes" was employed to release the allocated memories.
