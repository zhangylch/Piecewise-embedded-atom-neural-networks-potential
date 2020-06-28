#Macro definition

FC = ifort
# -pg

 FFLAGS = -O5  -xhost -ipo #   -CB -g -traceback -debug parallel
 LDFLAGS = -L /share/software/compiler/intel/intel-compiler-2017.5.239/mkl -lmkl_rt

#FFLAGS = -O2 -openmp -static
#LDFLAGS = -L/work1/soft/intel/mkl/10.1.2.024/lib/em64t -lmkl_em64t -lmkl_lapack -lguide

#MKLROOT = /opt/intel/composer_xe_2015.3.187/mkl
#FFLAGS = -O3 -fast -openmp -static -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
#LDFLAGS =  $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
#           $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a \
#           -Wl,--start-group \
#           $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
#           $(MKLROOT)/lib/intel64/libmkl_core.a \
#           $(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
#           -Wl,--end-group -lpthread -lm

#FFLAGS = -O2 -openmp
#MKLROOT = /work1/private/chenjun/soft/intel/mkl
#LDFLAGS = -L$(MKLROOT)/lib/intel64 -lmkl_rt -Wl,-rpath $(MKLROOT)/lib/intel64 

#end of Macro definition


OBJ_NN = mod.o test.o init_pes.o readinput.o get_index.o allocate_all.o allocate_wave.o  cart_to_frac.o readnn_w.o EANN_out.o get_energy.o get_force.o singleNN.o get_wave.o get_wave_force.o period.o deallocate_all.o der_NN.o
nn.exe: $(OBJ_NN)
	${FC} ${FFLAGS} $^ -o $@ ${LDFLAGS}
test.o: test.f90  
	${FC} ${FFLAGS} -c test.f90 
mod.o: mod.f90
	${FC} ${FFLAGS} -c mod.f90
readinput.o: readinput.f90    
	${FC} ${FFLAGS} -c readinput.f90  
get_index.o: get_index.f90    
	${FC} ${FFLAGS} -c get_index.f90  
allocate_all.o: allocate_all.f90
	${FC} ${FFLAGS} -c allocate_all.f90
allocate_wave.o: allocate_wave.f90
	${FC} ${FFLAGS} -c allocate_wave.f90
init_pes.o: init_pes.f90      
	${FC} ${FFLAGS} -c init_pes.f90   
cart_to_frac.o: cart_to_frac.f90
	${FC} ${FFLAGS} -c cart_to_frac.f90
readnn_w.o: readnn_w.f90
	${FC} ${FFLAGS} -c readnn_w.f90
EANN_out.o: EANN_out.f90
	${FC} ${FFLAGS} -c EANN_out.f90
get_energy.o: get_energy.f90
	${FC} ${FFLAGS} -c get_energy.f90
get_force.o: get_force.f90              
	${FC} ${FFLAGS} -c get_force.f90     
#get_hess.o: get_hess.f90              
#	${FC} ${FFLAGS} -c get_hess.f90     
der_NN.o: der_NN.f90              
	${FC} ${FFLAGS} -c der_NN.f90     
singleNN.o: singleNN.f90                
	${FC} ${FFLAGS} -c singleNN.f90      
get_wave.o: get_wave.f90
	${FC} ${FFLAGS} -c get_wave.f90
get_wave_force.o: get_wave_force.f90
	${FC} ${FFLAGS} -c get_wave_force.f90
#get_wave_hess.o: get_wave_hess.f90
#	${FC} ${FFLAGS} -c get_wave_hess.f90
period.o: period.f90                
	${FC} ${FFLAGS} -c period.f90       
deallocate_all.o: deallocate_all.f90
	${FC} ${FFLAGS} -c deallocate_all.f90
modify_input.o: modify_input.f90
	${FC} ${FFLAGS} -c modify_input.f90




clean:
	rm -f nn.exe nn.o nn_w.o *.mod

#rm -f *.exe *.mod *.o fort.* *.out
rm:
	rm -f M?? W??.txt

# a.out; gprof a.out gmon.out > out

