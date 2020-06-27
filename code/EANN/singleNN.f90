subroutine NN(mnl,nhid,nl,w,z)
     use constant
     implicit none
     integer(kind=intype) :: mnl,nhid,nl(0:nhid+1)
     real(kind=typenum) :: w(0:mnl,mnl,nhid+1),z(0:mnl,0:nhid+1)
     real(kind=typenum),external :: ddot
     integer(kind=intype) :: ilayer,neu,num
       do ilayer=1,nhid+1
         do neu=1,nl(ilayer) 
           num=nl(ilayer-1)+1
!           z(neu,ilayer)=dot_product(z(0:nl(ilayer-1),ilayer-1),w(0:nl(ilayer-1),neu,ilayer))
           z(neu,ilayer)=ddot(num,z(0:nl(ilayer-1),ilayer-1),1,w(0:nl(ilayer-1),neu,ilayer),1)
           if(ilayer<nhid+1) z(neu,ilayer)=dtanh(z(neu,ilayer))
         end do
       end do
     return
end subroutine
