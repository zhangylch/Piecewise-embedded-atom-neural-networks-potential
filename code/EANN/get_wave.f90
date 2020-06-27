subroutine get_wave(ntype1,scutnum,index_dis,vec,dis,wf)
     use constant
     use nnmod
     implicit none
     integer(kind=intype) :: ntype1,scutnum,index_dis(cutnum)
     integer(kind=intype) :: l,num,i,j,k,m,ntype,countnum,numpara,numorbit,numorbit1,numwave
     real(kind=typenum) :: vec(atomdim,cutnum),dis(cutnum),power_vec(scutnum,totpara)
     real(kind=typenum) :: tmp1,tmpwf(numtype)
     real(kind=typenum) :: wf(norbit),gauss_wave(scutnum,0:maxnwave)
!f2py integer(kind=intype),intent(aux) :: cutnum,totpara,totpara,numtype,maxnwave,norbit
       call gauss_orbital(ntype1,scutnum,index_dis,dis,gauss_wave)
       call wave_factor(scutnum,vec,power_vec)
       numwave=0
       numpara=0
       numorbit=0
       do m=0,ipsin
         numpara=numpara+1
         numorbit1=numorbit
         countnum=numwave
         do k=0,maxnwave
           tmpwf=0d0
           do num=1,scutnum
             j=index_dis(num)
             ntype=index_ele(j)
             if(k<=nwave(ntype)) then
               tmpwf(ntype)=tmpwf(ntype)+gauss_wave(num,k)*power_vec(num,numpara)
             end if
           end do
           numorbit1=numorbit1+1
           tmp1=tmpwf(1)*weight_wave(countnum+1,ntype1)
           do i=2,numtype
             tmp1=tmp1+tmpwf(i)*weight_wave(countnum+i,ntype1)
           end do
           wf(numorbit1)=tmp1*tmp1*factor_wave(numpara)
           countnum=countnum+numtype
         end do
         do l=2,npara(m)
           numpara=numpara+1
           numorbit1=numorbit
           countnum=numwave
           do k=0,maxnwave
             tmpwf=0d0
             do num=1,scutnum
               j=index_dis(num)
               ntype=index_ele(j)
               if(k<=nwave(ntype)) then
                 tmpwf(ntype)=tmpwf(ntype)+gauss_wave(num,k)*power_vec(num,numpara)
               end if
             end do
             numorbit1=numorbit1+1
             tmp1=tmpwf(1)*weight_wave(countnum+1,ntype1)
             do i=2,numtype
               tmp1=tmp1+tmpwf(i)*weight_wave(countnum+i,ntype1)
             end do
             wf(numorbit1)=wf(numorbit1)+tmp1*tmp1*factor_wave(numpara)
             countnum=countnum+numtype
           end do
         end do
         numwave=numwave+numtype*(maxnwave+1)
         numorbit=numorbit+maxnwave+1
       end do
     return
end subroutine
