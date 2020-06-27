subroutine get_wave_force(natom,ntype1,scutnum,index_dis,vec,dis,density,pdensity,pcenden)
     use constant
     use nnmod
     implicit none
     integer(kind=intype) natom,ntype1,scutnum
     integer(kind=intype) i,l,num,j,k,m,n,ntype,numpcut,numpcut1
     integer(kind=intype) nnorbit,inorbit,fnnorbit,nnwave,fnnwave,nangular
     integer(kind=intype) index_dis(cutnum)
     real(kind=typenum) power_vec(scutnum,totpara),der_power_vec(atomdim,scutnum,totpara)
     real(kind=typenum) vec(atomdim,cutnum),dis(cutnum)
     real(kind=typenum) wf(numtype),tmp,tmp1(atomdim),tmp2
     real(kind=typenum) gauss_wave(scutnum,0:maxnwave),der_gauss_wave(atomdim,scutnum,0:maxnwave)
     real(kind=typenum) density(norbit)
     real(kind=typenum) dwfdcoor(atomdim,scutnum),cendwf(atomdim)
     real(kind=typenum) pdensity(atomdim,pcutnum),pcenden(atomdim,norbit)
     real(kind=typenum) t6,t7
!f2py integer(kind=intype),intent(aux) :: cutnum,totpara,numtype,maxnwave,norbit,pcutnum
       call der_gauss(ntype1,scutnum,index_dis,vec,dis,gauss_wave,der_gauss_wave)
       call der_wave_factor(scutnum,vec,power_vec,der_power_vec)
       nnorbit=0
       nnwave=0
       nangular=0
       numpcut=0
       do m=0,ipsin
         nangular=nangular+1
         fnnorbit=nnorbit+1
         fnnwave=nnwave
         numpcut1=numpcut
!--------------------------for avoid the init of pdensity-------------------------
         do k=0,maxnwave
           inorbit=fnnorbit+k
           wf=0d0
           cendwf=0d0
           do l=1,scutnum
             j=index_dis(l)
             ntype=index_ele(j)
             if(k<=nwave(ntype)) then
               wf(ntype)=wf(ntype)+gauss_wave(l,k)*power_vec(l,nangular)
               tmp1=(gauss_wave(l,k)*der_power_vec(:,l,nangular)+der_gauss_wave(:,l,k)*power_vec(l,nangular))   &
               *weight_wave(fnnwave+ntype,ntype1)
               dwfdcoor(:,l)=tmp1
               cendwf=cendwf-tmp1
             end if
           end do
           tmp=wf(1)*weight_wave(fnnwave+1,ntype1)
           do ntype=2,numtype
             tmp=tmp+wf(ntype)*weight_wave(fnnwave+ntype,ntype1)
           end do
           tmp2=tmp*factor_wave(nangular)
           density(inorbit)=tmp*tmp2
           do l=1,scutnum
             numpcut1=numpcut1+1
             pdensity(:,numpcut1)=2d0*tmp2*dwfdcoor(:,l)
           end do
           pcenden(:,inorbit)=2d0*tmp2*cendwf
           fnnwave=fnnwave+numtype
         end do
!-------------------------------------------------------------------------------------------------
         do n=2,npara(m)
           nangular=nangular+1
           fnnorbit=nnorbit+1
           fnnwave=nnwave
           numpcut1=numpcut
           do k=0,maxnwave
             inorbit=fnnorbit+k
             wf=0d0
             cendwf=0d0
             do l=1,scutnum
               j=index_dis(l)
               ntype=index_ele(j)
               if(k<=nwave(ntype)) then
                 wf(ntype)=wf(ntype)+gauss_wave(l,k)*power_vec(l,nangular)
                 tmp1=(gauss_wave(l,k)*der_power_vec(:,l,nangular)+der_gauss_wave(:,l,k)*power_vec(l,nangular))   &
                 *weight_wave(fnnwave+ntype,ntype1)
                 dwfdcoor(:,l)=tmp1
                 cendwf=cendwf-tmp1
               end if
             end do
             tmp=wf(1)*weight_wave(fnnwave+1,ntype1)
             do ntype=2,numtype
               tmp=tmp+wf(ntype)*weight_wave(fnnwave+ntype,ntype1)
             end do
             tmp2=tmp*factor_wave(nangular)
             density(inorbit)=density(inorbit)+tmp*tmp2
             do i=1,scutnum
               numpcut1=numpcut1+1
               pdensity(:,numpcut1)=pdensity(:,numpcut1)+2d0*tmp2*dwfdcoor(:,i)
             end do
             pcenden(:,inorbit)=pcenden(:,inorbit)+2d0*tmp2*cendwf
             fnnwave=fnnwave+numtype
           end do
         end do
         nnorbit=nnorbit+maxnwave+1
         nnwave=nnwave+(maxnwave+1)*numtype
         numpcut=numpcut+(maxnwave+1)*scutnum
       end do
     return
end subroutine
