subroutine PEANN_out(table,start_force,coor,y,force)
     use constant
     use representation
     use sharedmod
     implicit none
     integer(kind=intype) :: table,start_force,i,j,num,sca(atomdim)
     real(kind=typenum) :: coor(atomdim,numatom),dire(atomdim,numatom,length),oriminv(3),bond
     real(kind=typenum) :: cart(atomdim,numatom),fcoor(atomdim,numatom)
     real(kind=typenum) :: y,force(numforce)
!f2py integer(kind=intype),intent(aux,hide) :: numatom,numforce
!f2py real(kind=typenum),intent(in),check(0==0) :: coor
!f2py real(kind=typenum),intent(inout) :: force
!f2py integer(kind=intype),intent(in) :: start_force,table
!f2py real(kind=typenum),intent(inout) :: y
       if(table==1) then
         fcoor=coor
         cart=matmul(scalmatrix,coor)
       else 
         cart=coor
         fcoor=matmul(inv_matrix,coor)
       end if
!------------------------shift the atom 1 to the origin---------------------------------------
       oriminv=1d5
       do j=1,atomdim
         if(cart(j,1)<oriminv(j)) oriminv(j)=cart(j,1)
       end do
       do i=2,numatom
         do j=1,atomdim
           bond=fcoor(j,i)-fcoor(j,1)
           sca(j)=nint(bond)
         end do
         cart(:,i)=cart(:,i)-sca(1)*scalmatrix(:,1)-sca(2)*scalmatrix(:,2)-sca(3)*scalmatrix(:,3)
         do j=1,atomdim
           if(cart(j,i)<oriminv(j)) oriminv(j)=cart(j,i)
         end do
       end do
       do num=1,length
         do i=1,numatom
           dire(:,i,num)=cart(:,i)+vector(:,num)
         end do
       end do
       oriminv=oriminv-hinc-rc
       do i=1,numatom
         cart(:,i)=cart(:,i)-oriminv
       end do
       index_rs=0
       do num=1,length
         do i=1,numatom
           dire(:,i,num)=dire(:,i,num)-oriminv
           if(dire(1,i,num)>0d0.and.dire(1,i,num)<rangecoor(1).and.dire(2,i,num)>0d0.and.dire(2,i,num)<rangecoor(2)   &
           .and.dire(3,i,num)>0d0.and.dire(3,i,num)<rangecoor(3)) then
             sca=ceiling(dire(:,i,num)/dier)
             index_rs(sca(1),sca(2),sca(3))=index_rs(sca(1),sca(2),sca(3))+1
             index_numrs(:,index_rs(sca(1),sca(2),sca(3)),sca(1),sca(2),sca(3))=[i,num]
           end if
         end do
       end do
!       write(*,*) index_rs
!------------------------periodic image and interaction cell--------------------------------------
       if(start_force==0) then
         call get_energy(dire,cart,y)
       else if(start_force==1) then
         call get_force(dire,cart,y,force)
       end if
     return
end subroutine
