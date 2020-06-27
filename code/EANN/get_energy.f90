subroutine get_energy(coor,dire,y)
     use constant
     use nnmod
     implicit none
     real(kind=typenum) :: coor(atomdim,numatom),dire(atomdim,numatom,length)
     real(kind=typenum) :: y
     real(kind=typenum) :: z(0:mnl,0:mnhid+1),tmp,wf(norbit)
     real(kind=typenum) :: vec(atomdim,cutnum),dis(cutnum)
     integer(kind=intype) :: natom,ikpoint,ntype,scutnum,index_dis(cutnum)
!f2py integer(kind=intype),intent(aux) :: numatom,length,outputneuron,nkpoint,cutnum,norbit,mnl,mnhid
       z(0,:)=1d0
       y=0d0
       do ikpoint=1,nkpoint
         tmp=0d0
!$omp parallel default(shared)
!$omp do private(natom,ntype,scutnum,vec,index_dis,dis,wf)  reduction(+:tmp) firstprivate(z,dire) schedule(dynamic,1)
         do natom=1,numatom
           ntype=index_ele(natom)
           call period(natom,ntype,coor(:,natom),dire,scutnum,vec,dis,index_dis)
           call get_wave(ntype,scutnum,index_dis,vec,dis,wf)
           z(1:norbit,0)=wf
           call NN(mnl,nhid(ntype),nl(0:nhid(ntype)+1,ntype),w(0:mnl,1:mnl,1:nhid(ntype)+1,ntype,ikpoint),z(0:mnl,0:nhid(ntype)+1))
           tmp=tmp+z(1,nhid(ntype)+1)
         end do
!$omp end do
!$omp end parallel
         y=tmp
       end do
     return
end subroutine
