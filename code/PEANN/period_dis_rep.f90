subroutine period_dis(natom,dire,cart,scutnum,index_dis,dis,effvec)
     use constant
     use representation
     use sharedmod
     implicit none
     integer(kind=intype) natom,ninit
     integer(kind=intype) i,j,num,i1,i2,i3,scutnum,index_dis(cutnum),sca(3),boder(2,3)
     real(kind=typenum) tmp,tmp1(atomdim)
     real(kind=typenum) cart(atomdim),dire(atomdim,numatom,length)
     real(kind=typenum) effvec(atomdim,cutnum),dis(cutnum)
!f2py integer(kind=intype),intent(aux) :: numatom,length,cutnum
       sca=ceiling(cart/dier)
       ninit=(length+1)/2
       dire(:,natom,ninit)=100d0
       do i=1,3
         boder(1,i)=max(1,sca(i)-interaction)
         boder(2,i)=min(numrs(i),sca(i)+interaction)
       end do
       scutnum=0
       do i3=boder(1,3),boder(2,3)
         do i2=boder(1,2),boder(2,2)
           do i1=boder(1,1),boder(2,1)
             do i=1,index_rs(i1,i2,i3)
               j=index_numrs(1,i,i1,i2,i3)
               num=index_numrs(2,i,i1,i2,i3)
               tmp1=dire(:,j,num)-cart
               tmp=dot_product(tmp1,tmp1)
               if(tmp<=rcsq) then
                 scutnum=scutnum+1
                 effvec(:,scutnum)=tmp1
                 dis(scutnum)=dsqrt(tmp)
                 index_dis(scutnum)=j
               end if
             end do
           end do
         end do
       end do
       dire(:,natom,ninit)=cart
     return
end subroutine
