subroutine get_wave_hess(i,coor,dire,wf,dwfdcoor,wf_hess)
     use constant
     use nnmod
     implicit none
!     integer(kind=intype) :: i
!     integer(kind=intype) :: l,num,j,k,m,n,p,ntype,countnum,num1,h,iforce,jforce,inum,jnum,num2,ntype1
!     integer(kind=intype) :: scutnum,index_dis(cutnum)
!     real(kind=typenum) :: coor(atomdim),dire(atomdim,numatom,length)
!     real(kind=typenum) :: vec(atomdim,cutnum)
!     real(kind=typenum) :: tmp(atomdim),tmp1,tmp2,tmp4(atomdim),tmp5(numforce,numforce)
!     real(kind=typenum) :: wf(norbit),dwfdcoor(numforce,norbit),wf_hess(numforce,numforce,norbit)
!     real(kind=typenum) :: tmpwf(maxnpara),tmpdwf(3*numatom,maxnpara)
!     real(kind=typenum) :: tmphess(numforce,numforce,maxnpara)
!     real(kind=typenum) :: gauss_wave(cutnum,0:maxnwave),der_gauss_wave(atomdim,cutnum,0:maxnwave)
!     real(kind=typenum) :: power_vec(-2:ipsin,-2:ipsin,-2:ipsin,cutnum)
!     real(kind=typenum) :: der_power_vec(atomdim,0:ipsin,0:ipsin,0:ipsin,cutnum)
!     real(kind=typenum) :: secder_gauss_wave(numforce,numforce,cutnum,0:maxnwave)
!     real(kind=typenum) :: secder_power_vec(numforce,numforce,0:ipsin,0:ipsin,0:ipsin,cutnum)
!       call secder_gauss(i,coor,dire,scutnum,index_dis,vec,gauss_wave,der_gauss_wave,secder_gauss_wave)
!       call secder_wave_factor(i,scutnum,index_dis,vec,power_vec,der_power_vec,secder_power_vec)
!       ntype1=index_ele(i)
!       wf=0d0
!       dwfdcoor=0d0
!       wf_hess=0d0
!       num2=0
!       countnum=0
!       do m=0,ipsin
!         do k=0,maxnwave
!           num2=num2+1
!           tmpwf(1:npara(m))=0d0
!           tmpdwf(:,1:npara(m))=0d0
!           tmphess(:,:,1:npara(m))=0d0
!           do num=1,scutnum
!             j=index_dis(num)
!             ntype=index_ele(j)
!             if(k<=nwave(ntype)) then 
!               tmp=der_gauss_wave(:,num,k)
!               tmp1=gauss_wave(num,k)
!               tmp2=weight_wave(countnum+ntype,ntype1)
!               num1=0
!               do n=0,nipsin(m)
!                 do p=0,nipsin(m)-n
!                   num1=num1+1
!                   h=nipsin(m)-n-p
!                   tmpwf(num1)=tmpwf(num1)+tmp1*power_vec(h,p,n,num)*tmp2   
!                   tmp4=(tmp*power_vec(h,p,n,num)+tmp1*der_power_vec(:,h,p,n,num))*tmp2
!!----------------------------------------the first derivate----------------------------------------
!                   tmpdwf((j-1)*3+1:j*3,num1)=tmpdwf((j-1)*3+1:j*3,num1)+tmp4
!                   tmpdwf((i-1)*3+1:i*3,num1)=tmpdwf((i-1)*3+1:i*3,num1)-tmp4
!                   tmp5=0d0
!                   if(j<=neff.or.i<=neff) then
!                     inum=(i-1)*atomdim
!                     jnum=(j-1)*atomdim
!                     if(j<=neff) then
!                       do iforce=1,atomdim
!                         do jforce=1,atomdim
!                           tmp5(jnum+jforce,jnum+iforce)=tmp(jforce)*der_power_vec(iforce,h,p,n,num)
!                         end do
!                       end do
!                     end if
!                     if(i<=neff) then
!                       do iforce=1,atomdim
!                         do jforce=1,atomdim
!                           tmp5(inum+jforce,inum+iforce)=tmp(jforce)*der_power_vec(iforce,h,p,n,num)
!                         end do
!                       end do
!                     end if
!                     if(i<=neff.and.j<=neff) then
!                       tmp5(inum+1:inum+atomdim,jnum+1:jnum+atomdim)=   &
!                       -tmp5(jnum+1:jnum+atomdim,jnum+1:jnum+atomdim)
!                       tmp5(jnum+1:jnum+atomdim,inum+1:inum+atomdim)=       &
!                       tmp5(inum+1:inum+atomdim,jnum+1:jnum+atomdim)
!                     end if
!                     do jforce=1,numforce
!                       do iforce=1,jforce
!                         tmp5(iforce,jforce)=tmp5(iforce,jforce)+tmp5(jforce,iforce)
!                       end do
!                     end do
!                     do jforce=1,numforce-1
!                       do iforce=jforce+1,numforce
!                         tmp5(iforce,jforce)=tmp5(jforce,iforce)
!                       end do
!                     end do
!                     tmphess(:,:,num1)=tmphess(:,:,num1)+(tmp5+tmp1*secder_power_vec(:,:,h,p,n,num) &
!                     +secder_gauss_wave(:,:,num,k)*power_vec(h,p,n,num))*tmp2
!                   end if
!                 end do
!               end do 
!             end if
!           end do 
!           do num1=1,index_para(num2)
!             h=index_power(1,num1,num2)
!             p=index_power(2,num1,num2)
!             n=index_power(3,num1,num2)
!             tmp2=2d0*tmpwf(num1)
!             wf(num2)=wf(num2)+tmpwf(num1)*tmpwf(num1)*factor_wave(h,p,n)
!             dwfdcoor(:,num2)=dwfdcoor(:,num2)+tmp2*tmpdwf(1:numforce,num1)*factor_wave(h,p,n)
!             do jforce=1,numforce
!               do iforce=1,jforce
!                 wf_hess(iforce,jforce,num2)=wf_hess(iforce,jforce,num2)+factor_wave(h,p,n)   &
!                 *(tmp2*tmphess(iforce,jforce,num1)+2d0*tmpdwf(iforce,num1)*tmpdwf(jforce,num1))
!               end do
!             end do
!           end do
!           do jforce=1,numforce
!             do iforce=1,jforce
!               wf_hess(jforce,iforce,num2)=wf_hess(iforce,jforce,num2)
!             end do
!           end do
!           countnum=countnum+numtype
!         end do
!       end do
!     return
!end subroutine
