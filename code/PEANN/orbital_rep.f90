!=================================================================================================
subroutine wave_factor(scutnum,effvec,power_vec)
     use constant
     use representation
     use sharedmod
     implicit none
     integer(kind=intype) i,j,k,l,m,n,p,num,scutnum
     real(kind=typenum) power_vec(scutnum,totpara),effvec(atomdim,cutnum)    
     real(kind=typenum) tmp_vec(atomdim,scutnum,0:ipsin)
!f2py integer(kind=intype),intent(aux) :: cutnum,ipsin,totpara
       tmp_vec(:,:,0)=1d0
       do m=1,ipsin
         do l=1,scutnum
           do k=1,atomdim
             tmp_vec(k,l,m)=tmp_vec(k,l,m-1)*effvec(k,l)
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
subroutine der_wave_factor(scutnum,effvec,power_vec,der_power_vec)
     use constant
     use representation
     use sharedmod
     implicit none
     integer(kind=intype) i,j,k,l,m,n,p,num,num1,num2,num3,scutnum
     real(kind=typenum) power_vec(scutnum,totpara),der_power_vec(atomdim,scutnum,totpara),effvec(atomdim,cutnum)
     real(kind=typenum) tmp_vec(atomdim,scutnum,0:ipsin)
!f2py integer(kind=intype),intent(aux) :: totpara,cutnum,ipsin
       tmp_vec(:,:,0)=1d0
       do m=1,ipsin
         do l=1,scutnum
           do k=1,atomdim
             tmp_vec(k,l,m)=tmp_vec(k,l,m-1)*effvec(k,l)
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
