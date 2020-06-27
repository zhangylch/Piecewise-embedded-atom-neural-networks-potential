subroutine EANN_out(table,start_force,coor,y,force)
     use constant
     use nnmod
     implicit none
     integer(kind=intype) :: num,i,j,k,table,start_force,sca(3)
     real(kind=typenum) :: oriminv(3),fcoor(atomdim,numatom),bond
     real(kind=typenum) :: coor(atomdim,numatom),cart(atomdim,numatom),dire(atomdim,numatom,length)
     real(kind=typenum) :: y,force(numforce),hess(numforce,numforce)
!f2py integer(kind=intype),intent(aux,hide) :: numatom,numforce
!f2py real(kind=typenum),intent(in),check(0==0) :: coor
!f2py real(kind=typenum),intent(inout) :: force
!f2py integer(kind=intype),intent(in) :: start_force,table
!f2py real(kind=typenum),intent(inout) :: y
       if(table==0) then
         cart=coor
         fcoor=matmul(inv_matrix,coor)
       else 
         fcoor=coor
         cart=matmul(scalmatrix,coor)
       end if
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
       oriminv=oriminv-hinc-maxrc
       do i=1,numatom
         cart(:,i)=cart(:,i)-oriminv
       end do
       index_rs=0
       do num=1,length
         do i=1,numatom
           dire(:,i,num)=dire(:,i,num)-oriminv
           if(dire(1,i,num)>0d0.and.dire(1,i,num)<rangecoor(1).and.dire(2,i,num)>0d0.and.dire(2,i,num)<rangecoor(2)  &
           .and.dire(3,i,num)>0d0.and.dire(3,i,num)<rangecoor(3)) then
             sca=ceiling(dire(:,i,num)/dier)
             index_rs(sca(1),sca(2),sca(3))=index_rs(sca(1),sca(2),sca(3))+1
             index_numrs(:,index_rs(sca(1),sca(2),sca(3)),sca(1),sca(2),sca(3))=[i,num] 
           end if
         end do
       end do
       if(start_force==0) then
         call get_energy(cart,dire,y)
       else if(start_force==1) then
         call get_force(cart,dire,y,force)
!       else 
!         call get_hess(cart,dire,y,force,hess)
       end if
     return
end subroutine
   
