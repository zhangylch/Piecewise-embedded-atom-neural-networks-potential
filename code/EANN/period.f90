subroutine period(natom,ntype,cart,dire,scutnum,vec,dis,index_dis)
     use constant 
     use nnmod
     implicit none 
     integer(kind=intype) :: natom,sca(3),boder(2,3) 
     integer(kind=intype) :: i,j,i1,i2,i3,num,ninit,scutnum 
     integer(kind=intype) :: ntype,index_dis(cutnum)
     real(kind=typenum) :: cart(atomdim),dire(atomdim,numatom,length) 
     real(kind=typenum) :: tmp1(atomdim),vec(atomdim,cutnum) 
     real(kind=typenum) :: tmp,dis(cutnum) 
!f2py integer(kind=intype),intent(aux) :: numatom,length,cutnum
       sca=ceiling(cart/dier)
       ninit=(length+1)/2
       dire(:,natom,ninit)=100d0
       do i=1,3
         boder(1,i)=max(1,sca(i)-interaction(ntype))
         boder(2,i)=min(numrs(i),sca(i)+interaction(ntype))
       end do
       scutnum=0
       do i3=boder(1,3),boder(2,3)
         do i2=boder(1,2),boder(2,2)
           do i1=boder(1,1),boder(2,1)
             do i=1,index_rs(i1,i2,i3)
               j=index_numrs(1,i,i1,i2,i3)
               num=index_numrs(2,i,i1,i2,i3)
               tmp1=dire(:,j,num)-cart
               tmp=dot_product(tmp1,tmp1)
               if(tmp<=rcsq(ntype)) then
                 scutnum=scutnum+1
                 vec(:,scutnum)=tmp1
                 dis(scutnum)=dsqrt(tmp)
                 index_dis(scutnum)=j
               end if
             end do
           end do
         end do
       end do
       dire(:,natom,ninit)=cart
     return
end subroutine
!====================================================================================
subroutine effectcos(ntype,scutnum,dis,effect)
     use constant
     use nnmod
     implicit none
     real(kind=typenum) :: dis(cutnum)
     real(kind=typenum) :: effect(cutnum),tmp
     integer(kind=intype) :: ntype,j,num,l,scutnum!i is the centre atom
!f2py integer(kind=intype),intent(aux) :: cutnum
       do num=1,scutnum
         tmp=0.5d0*dcos(pi*dis(num)/rc(ntype))+0.5d0
         effect(num)=tmp*tmp
       end do
     return
end subroutine
!==================================================================================================
subroutine deffect(ntype,scutnum,dis,effect,dedr)
     use constant
     use nnmod
     implicit none
     real(kind=typenum) :: dis(cutnum),tmp
     real(kind=typenum) :: dedr(cutnum),effect(cutnum)
     integer(kind=intype) :: ntype,scutnum
     integer(kind=intype) :: j,l,num !k is the centre atom
!f2py integer(kind=intype),intent(aux) :: cutnum
       do num=1,scutnum
         tmp=0.5d0*dcos(pi*dis(num)/rc(ntype))+0.5d0
         effect(num)=tmp*tmp
         dedr(num)=-tmp*dsin(pi*dis(num)/rc(ntype))*pi/rc(ntype)
       end do
     return
end subroutine
!=============================================================================================================
subroutine dbondlendx(scutnum,vec,dis,drdx)
     use constant
     use nnmod
     implicit none
     real(kind=typenum) :: dis(cutnum)
     real(kind=typenum) :: drdx(atomdim,cutnum),vec(atomdim,cutnum)
     integer(kind=intype) :: scutnum
     integer(kind=intype) :: j,l,num! dr/dxi
!f2py integer(kind=intype),intent(aux) :: cutnum
       do num=1,scutnum
         drdx(:,num)=vec(:,num)/dis(num)
       end do
     return
end subroutine
!----------------------------------------------------------------------
subroutine gauss_orbital(ntype1,scutnum,index_dis,dis,gauss_wave)
     use constant
     use nnmod
     implicit none
     integer(kind=intype) :: scutnum,j,l,num,k,ntype,index_dis(cutnum),ntype1
     real(kind=typenum) :: dis(cutnum),effect(cutnum),tmp
     real(kind=typenum) :: gauss_wave(scutnum,0:maxnwave),vec(atomdim,cutnum)
!f2py integer(kind=intype),intent(aux) :: cutnum,maxnwave
       call effectcos(ntype1,scutnum,dis,effect)
       do k=0,maxnwave
         do num=1,scutnum
           j=index_dis(num)
           ntype=index_ele(j)
           if(k<=nwave(ntype)) then
             tmp=dis(num)-rs(ntype,k)
             gauss_wave(num,k)=dexp(-inta(ntype,k)*tmp*tmp)*effect(num)
           end if
         end do
       end do  
     return
end subroutine
!----------------------------------------------------------------------
subroutine der_gauss(ntype1,scutnum,index_dis,vec,dis,gauss_wave,der_gauss_wave)
     use constant
     use nnmod
     implicit none
     integer(kind=intype) :: scutnum,index_dis(cutnum),ntype1
     real(kind=typenum) :: vec(atomdim,cutnum),dis(cutnum)
     real(kind=typenum) :: effect(cutnum),dedr(cutnum),drdx(atomdim,cutnum)
     real(kind=typenum) :: gauss_wave(scutnum,0:maxnwave),der_gauss_wave(atomdim,scutnum,0:maxnwave)
     real(kind=typenum) :: tmp,tmp1,tmp2
     integer(kind=intype) :: j,l,num,k,ntype
!f2py integer(kind=intype),intent(aux) :: cutnum,maxnwave
       call deffect(ntype1,scutnum,dis,effect,dedr)
       call dbondlendx(scutnum,vec,dis,drdx)
       do k=0,maxnwave
         do num=1,scutnum
           j=index_dis(num)
           ntype=index_ele(j)
           if(k<=nwave(ntype)) then
             tmp=dis(num)-rs(ntype,k)
             tmp1=-inta(ntype,k)*tmp
             tmp2=dexp(tmp1*tmp)
             gauss_wave(num,k)=tmp2*effect(num)
             der_gauss_wave(:,num,k)=(2d0*tmp1*gauss_wave(num,k)+tmp2*dedr(num))*drdx(:,num)
           end if
         end do
       end do  
     return
end subroutine
!----------------------------------------------------------------------------
subroutine wave_factor(scutnum,vec,power_vec)
     use constant
     use nnmod
     implicit none
     integer(kind=intype) :: scutnum,i,j,k,l,m,n,p,num
     real(kind=typenum) :: vec(atomdim,cutnum),tmp_vec(atomdim,scutnum,0:ipsin)
     real(kind=typenum) :: power_vec(scutnum,totpara)
!f2py integer(kind=intype),intent(aux) :: cutnum,ipsin,totpara
       tmp_vec(:,:,0)=1d0
       do m=1,ipsin
         do num=1,scutnum
           do k=1,atomdim
             tmp_vec(k,num,m)=tmp_vec(k,num,m-1)*vec(k,num)
           end do
         end do
       end do  
       power_vec(:,1)=1d0
       num=1
       do i=1,ipsin
         do j=1,npara(i)
           num=num+1
           p=index_power(1,num)
           n=index_power(2,num)
           m=index_power(3,num)
           do l=1,scutnum
             power_vec(l,num)=tmp_vec(1,l,p)*tmp_vec(2,l,n)*tmp_vec(3,l,m)
           end do
         end do
       end do
     return
end subroutine
!----------------------------------------------------------------------------
subroutine der_wave_factor(scutnum,vec,power_vec,der_power_vec)
     use constant
     use nnmod
     implicit none
     integer(kind=intype) i,j,k,l,m,n,p,num,num1,num2,num3,scutnum
     real(kind=typenum) power_vec(scutnum,totpara),der_power_vec(atomdim,scutnum,totpara),vec(atomdim,cutnum)
     real(kind=typenum) tmp_vec(atomdim,scutnum,0:ipsin)
!f2py integer(kind=intype),intent(aux) :: totpara,cutnum,ipsin
       tmp_vec(:,:,0)=1d0
       do m=1,ipsin
         do l=1,scutnum
           do k=1,atomdim
             tmp_vec(k,l,m)=tmp_vec(k,l,m-1)*vec(k,l)
           end do
         end do
       end do
       power_vec(:,1)=1d0
       der_power_vec(:,:,1)=0d0
       num=1
       do i=1,ipsin
         do j=1,npara(i)
           num=num+1
           p=index_power(1,num)
           n=index_power(2,num)
           m=index_power(3,num)
           do l=1,scutnum
             num1=inv_power(p-1,n,m)
             num2=inv_power(p,n-1,m)
             num3=inv_power(p,n,m-1)
             power_vec(l,num)=tmp_vec(1,l,p)*tmp_vec(2,l,n)*tmp_vec(3,l,m)
             der_power_vec(1,l,num)=dble(p)*power_vec(l,num1)
             der_power_vec(2,l,num)=dble(n)*power_vec(l,num2)
             der_power_vec(3,l,num)=dble(m)*power_vec(l,num3)
           end do
         end do
       end do
     return
end subroutine
!-----------------------------------------------------------------------------------
!subroutine secdbondlendx(i,scutnum,vec,dis,index_dis,secdrdx)
!      use constant
!      use nnmod
!      implicit none
!      integer(kind=intype) :: i,j,l,m,n,k! dr/dxi
!      integer(kind=intype) :: num,inum,jnum
!      integer(kind=intype) :: scutnum,index_dis(cutnum)
!      real(kind=typenum) :: dis(cutnum),secdrdx(numforce,numforce,cutnum)
!      real(kind=typenum) :: vec(atomdim,cutnum),tmp
!        secdrdx=0d0
!        inum=(i-1)*atomdim
!        do num=1,cutnum
!          j=index_dis(num)
!          if(i<=neff.or.j<=neff) then
!            tmp=dis(num)*dis(num)*dis(num)
!            if(i<=neff) then
!              inum=(i-1)*atomdim
!              do n=1,atomdim
!                secdrdx(inum+n,inum+n,num)=(dis(num)*dis(num)-vec(n,num)*vec(n,num))/tmp
!              end do
!              do n=1,atomdim-1
!                do k=n+1,atomdim
!                  secdrdx(inum+k,inum+n,num)=-vec(n,num)*vec(k,num)/tmp
!                end do
!              end do
!              do n=2,atomdim
!                do k=1,n-1
!                  secdrdx(inum+k,inum+n,num)=secdrdx(inum+n,inum+k,num)
!                end do
!              end do 
!            end if
!            if(j<=neff) then
!              jnum=(j-1)*atomdim
!              do n=1,atomdim
!                secdrdx(jnum+n,jnum+n,num)=(dis(num)*dis(num)-vec(n,num)*vec(n,num))/tmp
!              end do
!              do n=1,atomdim-1
!                do k=n+1,atomdim
!                  secdrdx(jnum+k,jnum+n,num)=-vec(n,num)*vec(k,num)/tmp
!                end do
!              end do
!              do n=2,atomdim
!                do k=1,n-1
!                  secdrdx(jnum+k,jnum+n,num)=secdrdx(jnum+n,jnum+k,num)
!                end do
!              end do 
!            end if
!            if(i<=neff.and.j<=neff) then
!              inum=(i-1)*atomdim
!              jnum=(j-1)*atomdim
!              secdrdx(inum+1:inum+atomdim,jnum+1:jnum+atomdim,num)=-secdrdx(inum+1:inum+atomdim,inum+1:inum+atomdim,num)
!              secdrdx(jnum+1:jnum+atomdim,inum+1:inum+atomdim,num)=secdrdx(inum+1:inum+atomdim,jnum+1:jnum+atomdim,num)
!            end if
!          end if
!        end do
!      return
!end subroutine
!!==========================================================================================================================
!subroutine d2bondlendx(i,vec,dis,scutnum,index_dis,drdx,d2rdx)
!     use constant
!     use nnmod
!     implicit none
!     integer(kind=intype) :: i,scutnum,index_dis(cutnum)
!     real(kind=typenum) :: dis(cutnum),d2rdx(numforce,numforce,cutnum)
!     real(kind=typenum) :: drdx(atomdim,cutnum),vec(atomdim,cutnum)
!     integer(kind=intype) :: j,l,num,m,n,inum,jnum,p,k! dr/dxi
!       d2rdx=0d0
!       do num=1,length
!         j=index_dis(num)
!         drdx(:,num)=vec(:,num)/dis(num)
!         if(i<=neff) then
!           inum=(i-1)*atomdim
!           do k=1,atomdim
!             do p=1,k
!               d2rdx(inum+p,inum+k,num)=drdx(k,num)*drdx(p,num)
!             end do
!           end do
!           do k=1,atomdim-1
!             do p=k,atomdim
!               d2rdx(inum+p,inum+k,num)=d2rdx(inum+k,inum+p,num)
!             end do
!           end do
!         end if
!         if(j<=neff) then
!           jnum=(j-1)*atomdim
!           do k=1,atomdim
!             do p=1,k
!               d2rdx(jnum+p,jnum+k,num)=drdx(k,num)*drdx(p,num)
!             end do
!           end do
!           do k=1,atomdim-1
!             do p=k+1,atomdim
!               d2rdx(jnum+p,jnum+k,num)=d2rdx(jnum+k,jnum+p,num)
!             end do
!           end do
!         end if
!         if(j<=neff.and.i<=neff) then 
!           jnum=(j-1)*atomdim
!           inum=(i-1)*atomdim
!           d2rdx(jnum+1:jnum+atomdim,inum+1:inum+atomdim,num)=-d2rdx(inum+1:inum+atomdim,inum+1:inum+atomdim,num)
!           d2rdx(inum+1:inum+atomdim,jnum+1:jnum+atomdim,num)=d2rdx(jnum+1:jnum+atomdim,inum+1:inum+atomdim,num)
!         end if
!       end do
!     return
!end subroutine
!!==================================================================================================
!subroutine secdeffect(ntype,dis,scutnum,effect,dedr,secdedr)
!     use constant
!     use nnmod
!     implicit none
!     integer(kind=intype) :: scutnum,ntype
!     real(kind=typenum) :: dis(cutnum),secdedr(cutnum),tmp,tmp1,tmp2,tmp3
!     real(kind=typenum) :: dedr(cutnum),effect(cutnum)
!     integer(kind=intype) :: j,l,num !k is the centre atom
!       do num=1,scutnum
!         tmp=pi/rc(ntype)
!         tmp2=dcos(tmp*dis(num))
!         tmp3=dsin(tmp*dis(num))
!         tmp1=0.5d0*tmp2+0.5d0
!         effect(num)=tmp1*tmp1
!         dedr(num)=-tmp1*tmp3*tmp
!         secdedr(num)=(-tmp1*tmp2+0.5d0*tmp3*tmp3)*tmp*tmp
!       end do
!     return
!end subroutine
!!----------------------------------------------------------------------
!subroutine secder_gauss(i,coor,dire,scutnum,index_dis,vec,gauss_wave,der_gauss_wave,secder_gauss_wave)
!     use constant
!     use nnmod
!     implicit none
!     integer(kind=intype) :: i
!     integer(kind=intype) :: scutnum,index_dis(cutnum)
!     real(kind=typenum) :: coor(atomdim),dire(atomdim,numatom,length)
!     real(kind=typenum) :: dis(cutnum),effect(cutnum)
!     real(kind=typenum) :: dedr(cutnum),drdx(atomdim,cutnum),secdedr(cutnum)
!     real(kind=typenum) :: gauss_wave(cutnum,0:maxnwave),vec(atomdim,cutnum)
!     real(kind=typenum) :: der_gauss_wave(atomdim,cutnum,0:maxnwave)
!     real(kind=typenum) :: d2rdx(numforce,numforce,cutnum)
!     real(kind=typenum) :: secder_gauss_wave(numforce,numforce,cutnum,0:maxnwave)
!     real(kind=typenum) :: secdrdx(numforce,numforce,cutnum)
!     real(kind=typenum) :: tmp,tmp1,tmp2,tmp3,tmp4,tmp5
!     integer(kind=intype) :: j,l,num,m,k,ntype,p,jnum,inum
!       ntype=index_ele(i)
!       call period(i,ntype,coor,dire,scutnum,vec,dis,index_dis)
!       call secdeffect(ntype,dis,scutnum,effect,dedr,secdedr)
!       call d2bondlendx(i,vec,dis,scutnum,index_dis,drdx,d2rdx)
!       call secdbondlendx(i,scutnum,vec,dis,index_dis,secdrdx)
!       do k=0,maxnwave
!         do num=1,scutnum
!           j=index_dis(num)
!           ntype=index_ele(j)
!           if(k<=nwave(ntype)) then
!             tmp=dis(num)-rs(ntype,k)
!             tmp1=-inta(ntype,k)*tmp
!             tmp2=dexp(tmp1*tmp)
!             gauss_wave(num,k)=tmp2*effect(num)
!             tmp4=2d0*tmp1*gauss_wave(num,k)
!             tmp5=tmp4+dedr(num)*tmp2
!             der_gauss_wave(:,num,k)=tmp5*drdx(:,num)
!             if(i<=neff.or.j<=neff) then
!               jnum=(j-1)*atomdim
!               inum=(i-1)*atomdim
!               tmp3=2d0*(-inta(ntype,k)*gauss_wave(num,k)+tmp1*tmp5)+      &
!               2d0*tmp1*tmp2*dedr(num)+tmp2*secdedr(num)
!               secder_gauss_wave(:,:,num,k)=tmp3*d2rdx(:,:,num)+tmp5*secdrdx(:,:,num)
!             end if
!           end if
!         end do
!       end do
!     return
!end subroutine
!!----------------------------------------------------------------------------
!subroutine secder_wave_factor(iatom,scutnum,index_dis,vec,power_vec,der_power_vec,secder_power_vec)
!     use constant
!     use nnmod
!     implicit none
!     integer(kind=intype) :: iatom,scutnum,index_dis(cutnum)
!     real(kind=typenum) :: vec(atomdim,cutnum),tmp_vec(atomdim,0:ipsin,cutnum)
!     real(kind=typenum) :: power_vec(-2:ipsin,-2:ipsin,-2:ipsin,cutnum)
!     real(kind=typenum) :: der_power_vec(atomdim,0:ipsin,0:ipsin,0:ipsin,cutnum)
!     real(kind=typenum) :: secder_power_vec(numforce,numforce,0:ipsin,0:ipsin,0:ipsin,cutnum)
!     integer(kind=intype) :: i,j,k,l,m,n,p,num,inum,jnum
!       tmp_vec(:,0,:)=1d0
!       do num=1,scutnum
!         j=index_dis(num)
!         do m=1,ipsin
!           do k=1,atomdim
!             tmp_vec(k,m,num)=tmp_vec(k,m-1,num)*vec(k,num)
!           end do
!         end do
!       end do  
!       power_vec=1d0
!       secder_power_vec=0d0
!       do i=0,ipsin 
!         do num=1,scutnum
!           j=index_dis(num)
!           do m=0,nipsin(i)
!             do n=0,nipsin(i)-m
!               p=nipsin(i)-m-n
!               power_vec(p,n,m,num)=tmp_vec(1,p,num)*tmp_vec(2,n,num)*tmp_vec(3,m,num)
!               der_power_vec(1,p,n,m,num)=dble(p)*power_vec(p-1,n,m,num)
!               der_power_vec(2,p,n,m,num)=dble(n)*power_vec(p,n-1,m,num)
!               der_power_vec(3,p,n,m,num)=dble(m)*power_vec(p,n,m-1,num)
!               if(iatom<=neff.or.j<=neff) then
!                 if(j<=neff) then
!                   jnum=(j-1)*atomdim
!                   secder_power_vec(jnum+1,jnum+1,p,n,m,num)=dble(p-1)*dble(p)*power_vec(p-2,n,m,num)
!                   secder_power_vec(jnum+2,jnum+1,p,n,m,num)=dble(n)*dble(p)*power_vec(p-1,n-1,m,num)
!                   secder_power_vec(jnum+3,jnum+1,p,n,m,num)=dble(m)*dble(p)*power_vec(p-1,n,m-1,num)
!                   secder_power_vec(jnum+2,jnum+2,p,n,m,num)=dble(n-1)*dble(n)*power_vec(p,n-2,m,num)
!                   secder_power_vec(jnum+3,jnum+2,p,n,m,num)=dble(m)*dble(n)*power_vec(p,n-1,m-1,num)
!                   secder_power_vec(jnum+3,jnum+3,p,n,m,num)=dble(m-1)*dble(m)*power_vec(p,n,m-2,num)
!                   secder_power_vec(jnum+1,jnum+2,p,n,m,num)=secder_power_vec(jnum+2,jnum+1,p,n,m,num)
!                   secder_power_vec(jnum+1,jnum+3,p,n,m,num)=secder_power_vec(jnum+3,jnum+1,p,n,m,num)
!                   secder_power_vec(jnum+2,jnum+3,p,n,m,num)=secder_power_vec(jnum+3,jnum+2,p,n,m,num)
!                 end if
!                 if(iatom<=neff) then
!                   inum=(iatom-1)*atomdim
!                   secder_power_vec(inum+1,inum+1,p,n,m,num)=dble(p-1)*dble(p)*power_vec(p-2,n,m,num)
!                   secder_power_vec(inum+2,inum+1,p,n,m,num)=dble(n)*dble(p)*power_vec(p-1,n-1,m,num)
!                   secder_power_vec(inum+3,inum+1,p,n,m,num)=dble(m)*dble(p)*power_vec(p-1,n,m-1,num)
!                   secder_power_vec(inum+2,inum+2,p,n,m,num)=dble(n-1)*dble(n)*power_vec(p,n-2,m,num)
!                   secder_power_vec(inum+3,inum+2,p,n,m,num)=dble(m)*dble(n)*power_vec(p,n-1,m-1,num)
!                   secder_power_vec(inum+3,inum+3,p,n,m,num)=dble(m-1)*dble(m)*power_vec(p,n,m-2,num)
!                   secder_power_vec(inum+1,inum+2,p,n,m,num)=secder_power_vec(inum+2,inum+1,p,n,m,num)
!                   secder_power_vec(inum+1,inum+3,p,n,m,num)=secder_power_vec(inum+3,inum+1,p,n,m,num)
!                   secder_power_vec(inum+2,inum+3,p,n,m,num)=secder_power_vec(inum+3,inum+2,p,n,m,num)
!                 end if
!                 if(iatom<=neff.and.j<=neff) then
!                   inum=(iatom-1)*atomdim
!                   jnum=(j-1)*atomdim
!                   secder_power_vec(inum+1:inum+3,jnum+1:jnum+3,p,n,m,num)=                  &
!                   -secder_power_vec(jnum+1:jnum+3,jnum+1:jnum+3,p,n,m,num)
!                   secder_power_vec(jnum+1:jnum+3,inum+1:inum+3,p,n,m,num)=                  &
!                   secder_power_vec(inum+1:inum+3,jnum+1:jnum+3,p,n,m,num)
!                 end if
!               end if
!             end do
!           end do
!         end do
!       end do
!     return
!end subroutine
