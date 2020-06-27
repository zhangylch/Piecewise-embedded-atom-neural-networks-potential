subroutine get_index_rep
     use constant
     use representation
     implicit none
     integer(kind=intype) :: i,k,m,n,num1
       num1=0
       inv_power=1
       do i=0,ipsin
         do k=0,i
           do m=0,i-k
             n=i-k-m
             num1=num1+1
             inv_power(n,m,k)=num1
             index_power(1:3,num1)=[n,m,k]
           end do
         end do
       end do
     return
end subroutine
