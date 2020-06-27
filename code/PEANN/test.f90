program co2
     use constant
     use representation
     implicit none
     integer(kind=intype) :: i,j
     real(kind=typenum) :: coor(3,192),y,force(576)
       call init_pes
       open(100,file='1')
       do i=1,2
         read(100,*)
         do j=1,192
           read(100,*) coor(:,j)
         end do
         call PEANN_out(0,1,coor,y,force)
         write(*,*) y
         write(*,*) force
       end do
       close(100)
       call deallocate_all
end 
