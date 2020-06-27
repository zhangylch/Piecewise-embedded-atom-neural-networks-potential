subroutine get_index
     use constant
     use nnmod
     implicit none
     integer(kind=intype) :: i,j,m,ntype,k,num,l,n
       num=0
       do i=0,ipsin
         do j=0,maxnwave
           num=num+1
           index_orbit((num-1)*numtype+1:num*numtype)=num
         end do
       end do
       l=0
       inv_power=1
       do i=0,ipsin
         do k=0,i
           do m=0,i-k
             n=i-k-m
             l=l+1
             inv_power(n,m,k)=l
             index_power(1:3,l)=[n,m,k]
           end do
         end do
       end do
       do i=1,numatom
         do j=1,numtype
           if(atom(i)==atomtype(j)) index_ele(i)=j
         end do
       end do
     return
end subroutine
