subroutine get_force(dire,cart,y,force)
    use constant
    use representation
    use sharedmod
    implicit none
    integer(kind=intype) natom,ntype,m,i,j,k,fneuron,num,ikpoint,outneuron
    integer(kind=intype) index_natom(2,wavecutnum),index_nwave(0:nwave-1)
    real(kind=typenum) doutdg(norbit,outputneuron),dire(atomdim,numatom,length),cart(atomdim,numatom)
    real(kind=typenum) y,force(numforce)
    real(kind=typenum) tmp,tmp1(numforce),tmp2(atomdim),tmp3(atomdim)
    real(kind=typenum) density(norbit),pdensity(atomdim,wavecutnum)
!f2py integer(kind=intype),intent(aux) :: norbit,wavecutnum,outputneuron,nwave,numatom,length,numforce
      do ikpoint=1,nkpoint
        tmp=0d0
        tmp1=0d0
!$omp parallel default(shared)
!$omp do private(natom,ntype,doutdg,m,i,j,fneuron,num,outneuron,k,index_natom, &
!$omp index_nwave,tmp2,tmp3,pdensity,density) firstprivate(z,dire) reduction(+:tmp,tmp1)
        do natom=1,numatom
          ntype=index_ele(natom)
          call get_wave_force(natom,dire,cart(:,natom),index_natom,index_nwave,density,pdensity)
          z(1:norbit,0)=density
          call NN(mnl,nhid(ntype),nl(0:nhid(ntype)+1,ntype),w(0:mnl,1:mnl,1:nhid(ntype)+1,ntype,ikpoint),z(0:mnl,0:nhid(ntype)+1))
          tmp=tmp+z(1,nhid(ntype)+1)
          call der_NN(norbit,outputneuron,mnl,nhid(ntype),nl(:,ntype),z(:,0:nhid(ntype)+1),    &
          w(:,:,1:nhid(ntype)+1,ntype,ikpoint),doutdg)
          do outneuron=1,outputneuron
            tmp2=0d0
            num=0
            do m=0,ipsin
              do i=0,nwave-1
                do fneuron=1,index_nwave(i)
                  num=num+1
                  j=index_natom(1,num)
                  k=index_natom(2,num)
                  if(j<=neff) then
                    tmp3=-doutdg(k,outneuron)*pdensity(:,num)
                    tmp1(j*3-2:j*3)=tmp1(j*3-2:j*3)+tmp3
                    tmp2=tmp2-tmp3
                  end if
                end do
              end do
            end do
            tmp1(natom*3-2:natom*3)=tmp1(natom*3-2:natom*3)+tmp2  
          end do
        end do
!$omp end do
!$omp end parallel
        y=tmp
        force=tmp1
      end do
    return
end subroutine
