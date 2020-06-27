subroutine get_pdensity(natom,scutnum,index_dis,index_nwave,index_numdis,index_numdis1,   &
     effect,dedr,effvec,index_natom,index_tnwave,density,pdensity)
     use constant
     use representation
     use sharedmod
     implicit none
     integer(kind=intype) natom,scutnum
     integer(kind=intype) i,l,j,k,m,n,ntype,ntype1,count1
     integer(kind=intype) nnorbit,inorbit,fnnorbit,nnwave,fnnwave,nangular
     integer(kind=intype) npdensity,ipdensity,jpdensity,kpdensity,ncut1,ncut2
     integer(kind=intype) index_dis(cutnum),index_nwave(2,0:nwave-1),index_tnwave(0:nwave-1)
     integer(kind=intype) index_numdis(pcutnum),index_numdis1(pcutnum),index_natom(2,wavecutnum)
     real(kind=typenum) power_vec(scutnum,totpara),der_power_vec(atomdim,scutnum,totpara)
     real(kind=typenum) effvec(atomdim,cutnum)
     real(kind=typenum) wf(maxnumtype),tmp,tmp1(atomdim),tmp2
     real(kind=typenum) density(norbit)
     real(kind=typenum) dwfdcoor(atomdim,scutnum)
     real(kind=typenum) pdensity(atomdim,wavecutnum)
     real(kind=typenum) effect(pcutnum),dedr(atomdim,pcutnum)
!f2py integer(kind=intype),intent(aux) :: cutnum,totpara,maxnumtype,nwave,norbit,pcutnum,wavecutnum
       call der_wave_factor(scutnum,effvec,power_vec,der_power_vec)
       do k=0,nwave-1
         index_tnwave(k)=index_nwave(1,k)+index_nwave(2,k)
       end do
       nnorbit=0
       nnwave=0
       nangular=0
       npdensity=0
       ipdensity=0
       jpdensity=0
       ntype1=index_ele(natom)
       do m=0,ipsin
         nangular=nangular+1
         fnnorbit=nnorbit+1
         fnnwave=nnwave
!--------------------------for avoid the init of pdensity-------------------------
         ncut1=0
         ncut2=0
         do k=0,nwave-1
           inorbit=fnnorbit+k
           wf=0d0
           do i=1,index_nwave(1,k)
             ncut2=ncut2+1
             ipdensity=ipdensity+1
             l=index_numdis1(ncut2)
             j=index_dis(l)
             index_natom(:,ipdensity)=[j,inorbit]
             ntype=index_ele(j)
             wf(ntype)=wf(ntype)+power_vec(l,nangular)
             tmp1=der_power_vec(:,l,nangular)*weight_wave(fnnwave+ntype,ntype1)
             dwfdcoor(:,i)=tmp1
           end do
           count1=index_nwave(1,k)+1
           do i=count1,index_tnwave(k)
             ncut1=ncut1+1
             ipdensity=ipdensity+1
             l=index_numdis(ncut1)
             j=index_dis(l)
             index_natom(:,ipdensity)=[j,inorbit]
             ntype=index_ele(j)
             wf(ntype)=wf(ntype)+effect(ncut1)*power_vec(l,nangular)
             tmp1=(effect(ncut1)*der_power_vec(:,l,nangular)+dedr(:,ncut1)*power_vec(l,nangular))*weight_wave(fnnwave+ntype,ntype1)
             dwfdcoor(:,i)=tmp1
           end do
           tmp=wf(1)*weight_wave(fnnwave+1,ntype1)
           do ntype=2,maxnumtype
             tmp=tmp+wf(ntype)*weight_wave(fnnwave+ntype,ntype1)
           end do
           tmp2=tmp*factor_wave(nangular)
           density(inorbit)=tmp*tmp2
           do i=1,index_tnwave(k)
             jpdensity=jpdensity+1
             pdensity(:,jpdensity)=2d0*tmp2*dwfdcoor(:,i)
           end do
           fnnwave=fnnwave+maxnumtype
         end do
!-------------------------------------------------------------------------------------------------
         do n=2,npara(m)
           nangular=nangular+1
           kpdensity=npdensity
           fnnorbit=nnorbit+1
           fnnwave=nnwave
           ncut1=0
           ncut2=0
           do k=0,nwave-1
             inorbit=fnnorbit+k
             wf=0d0
             do i=1,index_nwave(1,k)
               ncut2=ncut2+1
               l=index_numdis1(ncut2)
               j=index_dis(l)
               ntype=index_ele(j)
               wf(ntype)=wf(ntype)+power_vec(l,nangular)
               tmp1=der_power_vec(:,l,nangular)*weight_wave(fnnwave+ntype,ntype1)
               dwfdcoor(:,i)=tmp1
             end do
             count1=index_nwave(1,k)+1
             do i=count1,index_tnwave(k)
               ncut1=ncut1+1
               l=index_numdis(ncut1)
               j=index_dis(l)
               ntype=index_ele(j)
               wf(ntype)=wf(ntype)+effect(ncut1)*power_vec(l,nangular)
               tmp1=(effect(ncut1)*der_power_vec(:,l,nangular)+dedr(:,ncut1)*power_vec(l,nangular))  &
               *weight_wave(fnnwave+ntype,ntype1)
               dwfdcoor(:,i)=tmp1
             end do
             tmp=wf(1)*weight_wave(fnnwave+1,ntype1)
             do ntype=2,maxnumtype
               tmp=tmp+wf(ntype)*weight_wave(fnnwave+ntype,ntype1)
             end do
             tmp2=tmp*factor_wave(nangular)
             density(inorbit)=density(inorbit)+tmp*tmp2
             do i=1,index_tnwave(k)
               kpdensity=kpdensity+1
               pdensity(:,kpdensity)=pdensity(:,kpdensity)+2d0*tmp2*dwfdcoor(:,i)
             end do
             fnnwave=fnnwave+maxnumtype
           end do
         end do
         nnorbit=nnorbit+nwave
         nnwave=nnwave+nwave*maxnumtype
         npdensity=ipdensity
       end do
     return
end subroutine
