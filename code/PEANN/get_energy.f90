subroutine get_energy(dire,cart,y)
     use constant
     use representation
     use sharedmod
     implicit none
     integer(kind=intype) :: natom,ntype,ikpoint
     real(kind=typenum) :: dire(atomdim,numatom,length),cart(atomdim,numatom),y
     real(kind=typenum) :: density(norbit),tmp
!f2py integer(kind=intype),intent(aux) :: numatom,length,norbit,mnl,mnhid
       do ikpoint=1,nkpoint
         tmp=0d0
!$omp parallel default(shared)
!$omp do private(natom,ntype,density) firstprivate(z,dire) reduction(+:tmp)
         do natom=1,numatom 
           ntype=index_ele(natom)
           call get_wave(natom,dire,cart(:,natom),density)
           z(1:norbit,0)=density
           call NN(mnl,nhid(ntype),nl(0:nhid(ntype)+1,ntype),w(0:mnl,1:mnl,1:nhid(ntype)+1,ntype,ikpoint),z(0:mnl,0:nhid(ntype)+1))
           tmp=tmp+z(1,1+nhid(ntype))           
         end do
!$omp end do
!$omp end parallel
         y=tmp
       end do
     return
end subroutine
