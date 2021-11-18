subroutine der_NN(norbit,outputneuron,mnl,nhid,nl,zt,w,dada)
     use constant
     implicit none
     integer(kind=intype) :: nhid,mnl,norbit,outputneuron
     integer(kind=intype) :: nl(0:nhid+1)
     integer(kind=intype) :: flayer,fneuron,neuron,outneu
     real(kind=typenum) :: zt(0:mnl,0:nhid+1),w(0:mnl,mnl,nhid+1)
     real(kind=typenum) :: dada(norbit,outputneuron)
     real(kind=typenum) :: tmp(norbit,mnl),tmp1(norbit,mnl)
!f2py integer(kind=intype),intent(in,hide) :: mnl,nhid,norbit
!f2py integer(kind=intype),intent(in) :: nl
!f2py real(kind=typenum),intent(in,out) :: zt
!f2py real(kind=typenum),intent(out) :: dada
!f2py real(kind=typenum),intent(in) :: w
     tmp1(:,1:nl(1))=w(1:nl(0),1:nl(1),1) 
     do flayer=1,nhid,1 
       do fneuron=1,nl(flayer)
         tmp(:,fneuron)=tmp1(:,fneuron)*(1d0-zt(fneuron,flayer)*zt(fneuron,flayer))
       end do
       call dgemm('N','N',nl(0),nl(flayer+1),nl(flayer),1d0,tmp(:,1:nl(flayer)),nl(0),  &
       w(1:nl(flayer),1:nl(flayer+1),flayer+1),nl(flayer),0d0,tmp1,nl(0))
     end do
     dada=tmp1(:,1:nl(nhid+1))
     return
end subroutine
subroutine secder_NN(norbit,outputneuron,mnl,nhid,nl,zt,w,dada,d2ada)
     use constant
     implicit none
     integer(kind=intype) :: nhid,mnl,norbit,outputneuron
     integer(kind=intype) :: nl(0:nhid+1)
     integer(kind=intype) :: flayer,fneuron,neuron,outneu
     real(kind=typenum) :: zt(0:mnl,0:nhid+1),w(0:mnl,mnl,nhid+1)
     real(kind=typenum) :: dflda(mnl,outputneuron,nhid),dadil(norbit,mnl,nhid)
     real(kind=typenum) :: dwadil(norbit,mnl,nhid+1),deract(mnl,nhid)
     real(kind=typenum) :: dada(norbit,outputneuron),d2ada(norbit,norbit,outputneuron)
     real(kind=typenum) :: tmp_array(mnl,outputneuron),tmp(norbit,mnl)
!f2py integer(kind=intype),intent(in,hide) :: mnl,nhid,norbit
!f2py integer(kind=intype),intent(in) :: nl
!f2py real(kind=typenum),intent(in,out) :: zt
!f2py real(kind=typenum),intent(out) :: dada,d2ada
!f2py real(kind=typenum),intent(in) :: w
     dwadil(:,:,1)=w(1:nl(0),:,1) 
     do flayer=1,nhid,1 
       do fneuron=1,nl(flayer)
         deract(fneuron,flayer)=1d0-zt(fneuron,flayer)*zt(fneuron,flayer)
         dadil(:,fneuron,flayer)=deract(fneuron,flayer)*dwadil(:,fneuron,flayer)
       end do
       call dgemm('N','N',nl(0),nl(flayer+1),nl(flayer),1d0,dadil(:,1:nl(flayer),flayer),nl(0),  &
       w(1:nl(flayer),1:nl(flayer+1),flayer+1),nl(flayer),0d0,dwadil(:,1:nl(flayer+1),flayer+1),nl(0))
     end do
     dada=dwadil(:,1:nl(nhid+1),nhid+1)
     dflda(1:nl(nhid),:,nhid)=w(1:nl(nhid),:,nhid+1)
     do flayer=nhid,2,-1
       do neuron=1,nl(nhid+1)
         tmp_array(:,neuron)=dflda(1:nl(flayer),neuron,flayer)*deract(1:nl(flayer),flayer)
       end do
       call dgemm('N','N',nl(flayer-1),nl(nhid+1),nl(flayer),1d0,w(1:nl(flayer-1),1:nl(flayer),flayer),&
       nl(flayer-1),tmp_array(1:nl(flayer),:),nl(flayer),  &
       0d0,dflda(1:nl(flayer-1),:,flayer-1),nl(flayer-1))
     end do
     d2ada=0d0
     do flayer=1,nhid
       do outneu=1,nl(nhid+1)
         do fneuron=1,nl(flayer)
           tmp_array(fneuron,outneu)=-2d0*zt(fneuron,flayer)*dflda(fneuron,outneu,flayer)
         end do
       end do
       do outneu=1,nl(nhid+1)
         do fneuron=1,nl(flayer)
           tmp(:,fneuron)=tmp_array(fneuron,outneu)*dadil(:,fneuron,flayer)
         end do
         call dgemm('N','t',nl(0),nl(0),nl(flayer),1d0,tmp(:,1:nl(flayer)),&
         nl(0),dwadil(:,1:nl(flayer),flayer),nl(0),1d0,d2ada(:,:,outneu),nl(0))
       end do
     end do
     return
end subroutine
