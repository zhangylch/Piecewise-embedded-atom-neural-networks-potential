#Macro definition

FC = ifort

 FFLAGS = -O5 -pg  #  -CB -g -traceback -debug parallel
 LDFLAGS = -lmkl_rt

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


OBJ_NN = constant.o representation.o sharedmod.o test.o get_index_rep.o init_pes.o period_dis_rep.o readinput_rep.o orbital_rep.o get_wave_rep.o get_pdensity_rep.o  allocate_all_rep.o get_wave_force_rep.o get_pdensity_rep.o readnn.o allocate_shared.o get_energy.o get_force.o PEANN_out.o cart_to_frac.o singleNN.o der_NN.o deallocate_all_rep.o piece.o  get_density_rep.o
nn.exe: $(OBJ_NN)
	${FC} ${FFLAGS} $^ -o $@ ${LDFLAGS}
constant.o: constant.f90
	${FC} ${FFLAGS} -c constant.f90
sharedmod.o: sharedmod.f90          
	${FC} ${FFLAGS} -c sharedmod.f90      
representation.o: representation.f90
	${FC} ${FFLAGS} -c representation.f90
test.o: test.f90
	${FC} ${FFLAGS} -c test.f90
init_pes.o: init_pes.f90
	${FC} ${FFLAGS} -c init_pes.f90
readinput_rep.o: readinput_rep.f90
	${FC} ${FFLAGS} -c readinput_rep.f90
get_index_rep.o: get_index_rep.f90
	${FC} ${FFLAGS} -c get_index_rep.f90
period_dis_rep.o: period_dis_rep.f90
	${FC} ${FFLAGS} -c period_dis_rep.f90
orbital_rep.o: orbital_rep.f90
	${FC} ${FFLAGS} -c orbital_rep.f90
get_wave_rep.o: get_wave_rep.f90
	${FC} ${FFLAGS} -c get_wave_rep.f90
get_density_rep.o: get_density_rep.f90
	${FC} ${FFLAGS} -c get_density_rep.f90
allocate_all_rep.o: allocate_all_rep.f90
	${FC} ${FFLAGS} -c allocate_all_rep.f90
get_wave_force_rep.o: get_wave_force_rep.f90
	${FC} ${FFLAGS} -c get_wave_force_rep.f90
get_pdensity_rep.o: get_pdensity_rep.f90
	${FC} ${FFLAGS} -c get_pdensity_rep.f90
readnn.o: readnn.f90
	${FC} ${FFLAGS} -c readnn.f90
allocate_shared.o: allocate_shared.f90
	${FC} ${FFLAGS} -c allocate_shared.f90
get_energy.o: get_energy.f90
	${FC} ${FFLAGS} -c get_energy.f90
get_force.o: get_force.f90
	${FC} ${FFLAGS} -c get_force.f90
PEANN_out.o: PEANN_out.f90
	${FC} ${FFLAGS} -c PEANN_out.f90
cart_to_frac.o: cart_to_frac.f90
	${FC} ${FFLAGS} -c cart_to_frac.f90
singleNN.o: singleNN.f90
	${FC} ${FFLAGS} -c singleNN.f90
der_NN.o: der_NN.f90
	${FC} ${FFLAGS} -c der_NN.f90
deallocate_all_rep.o: deallocate_all_rep.f90
	${FC} ${FFLAGS} -c deallocate_all_rep.f90
piece.o: piece.f90
	${FC} ${FFLAGS} -c piece.f90
