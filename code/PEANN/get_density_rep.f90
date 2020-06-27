subroutine get_density(natom,scutnum,index_dis,dis,effvec,density)
     use constant
     use representation
     use sharedmod
     implicit none
     integer(kind=intype) natom,scutnum
     integer(kind=intype) i,l,num,j,k,n,m,ntype,num1,num2,num3,num4,num5,num6,ntype1
     integer(kind=intype) index_nwave(2,0:nwave-1),index_dis(cutnum),index_numdis(pcutnum),index_numdis1(pcutnum)
     real(kind=typenum) power_vec(scutnum,totpara),dis(cutnum)
     real(kind=typenum) effvec(atomdim,cutnum)
     real(kind=typenum) wf(maxnumtype),effect(pcutnum),tmp
     real(kind=typenum) density(norbit)
!f2py integer(kind=intype),intent(aux) :: cutnum,totpara,maxnumtype,nwave,norbit,pcutnum
       ntype1=index_ele(natom)
       call cutoff(scutnum,ntype1,index_nwave,index_numdis,index_numdis1,dis,effect)
       call wave_factor(scutnum,effvec,power_vec)
       num=0
       num2=0
       num3=0
       do m=0,ipsin
         num3=num3+1
         num4=num2
         num1=num+1
         num5=0
         num6=0
         do k=0,nwave-1
           wf=0d0
           do i=1,index_nwave(1,k)
             num6=num6+1
             l=index_numdis1(num6)
             j=index_dis(l)
             ntype=index_ele(j)
             wf(ntype)=wf(ntype)+power_vec(l,num3)
           end do
           do i=1,index_nwave(2,k)
             num5=num5+1
             l=index_numdis(num5)
             j=index_dis(l)
             ntype=index_ele(j)
             wf(ntype)=wf(ntype)+effect(num5)*power_vec(l,num3)
           end do
           tmp=wf(1)*weight_wave(num4+1,ntype1)
           do ntype=2,maxnumtype
             tmp=tmp+wf(ntype)*weight_wave(num4+ntype,ntype1)
           end do
           density(num1+k)=tmp*tmp*factor_wave(num3)
           num4=num4+maxnumtype
         end do
         do n=2,npara(m)
           num3=num3+1
           num4=num2
           num1=num+1
           num5=0
           num6=0
           do k=0,nwave-1
             wf=0d0
             do i=1,index_nwave(1,k)
               num6=num6+1
               l=index_numdis1(num6)
               j=index_dis(l)
               ntype=index_ele(j)
               wf(ntype)=wf(ntype)+power_vec(l,num3)
             end do
             do i=1,index_nwave(2,k)
               num5=num5+1
               l=index_numdis(num5)
               j=index_dis(l)
               ntype=index_ele(j)
               wf(ntype)=wf(ntype)+effect(num5)*power_vec(l,num3)
             end do
             tmp=wf(1)*weight_wave(num4+1,ntype1)
             do ntype=2,maxnumtype
               tmp=tmp+wf(ntype)*weight_wave(num4+ntype,ntype1)
             end do
             density(num1+k)=density(num1+k)+tmp*tmp*factor_wave(num3)
             num4=num4+maxnumtype
           end do
         end do
         num=num+nwave
         num2=num2+nwave*maxnumtype
       end do
     return
end subroutine
