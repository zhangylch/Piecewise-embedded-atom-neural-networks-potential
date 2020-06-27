subroutine get_force(cart,dire,y,force)
    use constant
    use nnmod
    implicit none
    integer(kind=intype) natom,ntype,m,i,j,k,fneuron,num,ikpoint,outneuron,scutnum
    integer(kind=intype) index_dis(cutnum)
    real(kind=typenum) doutdg(norbit,outputneuron),dire(atomdim,numatom,length),cart(atomdim,numatom)
    real(kind=typenum) y,force(numforce),z(0:mnl,0:mnhid+1)
    real(kind=typenum) tmp,tmp1(numforce)
    real(kind=typenum) density(norbit),pdensity(atomdim,pcutnum),pcenden(atomdim,norbit) 
    real(kind=typenum) vec(atomdim,cutnum),dis(cutnum)
    real(kind=typenum) t6,t7
!f2py integer(kind=intype),intent(aux) :: norbit,cutnum,outputneuron,numatom,length,nkpoint,numforce,mnl,mnhid,pcutnum,
      z(0,:)=1d0
      do ikpoint=1,nkpoint
        tmp=0d0
        tmp1=0d0
!$omp parallel default(shared)
!$omp do private(natom,ntype,doutdg,m,i,j,fneuron,num,outneuron,k,scutnum,index_dis, &
!$omp vec,dis,pdensity,density,pcenden) firstprivate(z,dire) reduction(+:tmp,tmp1)
        do natom=1,numatom
          ntype=index_ele(natom)
          call period(natom,ntype,cart(:,natom),dire,scutnum,vec,dis,index_dis)
          call get_wave_force(natom,ntype,scutnum,index_dis,vec,dis,density,pdensity,pcenden)
          z(1:norbit,0)=density
          call NN(mnl,nhid(ntype),nl(0:nhid(ntype)+1,ntype),w(0:mnl,1:mnl,1:nhid(ntype)+1,ntype,ikpoint),z(0:mnl,0:nhid(ntype)+1))
          tmp=tmp+z(1,nhid(ntype)+1)
          call der_NN(norbit,outputneuron,mnl,nhid(ntype),nl(:,ntype),z(:,0:nhid(ntype)+1),    &
          w(:,:,1:nhid(ntype)+1,ntype,ikpoint),doutdg)
          do outneuron=1,outputneuron
            i=0
            do m=1,norbit
              do num=1,scutnum
                i=i+1
                j=index_dis(num)
                if(j<=neff) then
                  tmp1(j*3-2:j*3)=tmp1(j*3-2:j*3)-doutdg(m,outneuron)*pdensity(:,i)
                end if
              end do
            end do
            if (natom<=neff) then
              j=natom*3-2
              do k=1,norbit
                tmp1(j:j+2)=tmp1(j:j+2)-doutdg(k,outneuron)*pcenden(:,k)
              end do
            end if
          end do
        end do
!$omp end do
!$omp end parallel
        y=tmp
        force=tmp1
      end do
    return
end subroutine
