subroutine cutoff(scutnum,ntype,index_nwave,index_numdis,index_numdis1,dis,effect)
     use constant
     use representation
     use sharedmod
     implicit none
     real(kind=typenum) dis(cutnum),effect(pcutnum),tmp
     integer(kind=intype) ntype,num,scutnum,k,j,m,index_nwave(2,0:nwave-1)
     integer(kind=intype) index_numdis(pcutnum),index_numdis1(pcutnum)
!f2py integer(kind=intype),intent(aux) :: cutnum,pcutnum,nwave
       index_nwave=0
       j=0
       m=0
       do k=0,nwave-1
         do num=1,scutnum
           if(dis(num)<=initrs(k,ntype)) then
             m=m+1
             index_numdis1(m)=num
             index_nwave(1,k)=index_nwave(1,k)+1
           else if(dis(num)<cenrs(k,ntype)) then
             j=j+1
             index_numdis(j)=num
             index_nwave(2,k)=index_nwave(2,k)+1
             tmp=(dis(num)-initrs(k,ntype))*(dis(num)-finalrs(k,ntype))
             effect(j)=dexp(tmp*tmp*normrs(k,ntype))*expon(k,ntype)+expalpha(k,ntype)
           end if
         end do
       end do
     return
end subroutine
!======================================================================
subroutine der_cutoff(scutnum,ntype,index_nwave,index_numdis,index_numdis1,effvec,dis,effect,dedr)
     use constant
     use representation
     use sharedmod
     implicit none
     integer(kind=intype) num,scutnum,k,j,l,m,ntype
     integer(kind=intype) index_nwave(2,0:nwave-1),index_numdis(pcutnum),index_numdis1(pcutnum)
     real(kind=typenum) dis(cutnum),effect(pcutnum),dedr(atomdim,pcutnum),tmp,tmp1,tmp2,tmp3
     real(kind=typenum) drdx(atomdim,cutnum),effvec(atomdim,cutnum)
!f2py integer(kind=intype),intent(aux) :: cutnum,pcutnum,nwave
       do l=1,scutnum
         drdx(:,l)=effvec(:,l)/dis(l)
       end do
       index_nwave=0
       j=0
       m=0
       do k=0,nwave-1
         do num=1,scutnum
           if(dis(num)<=initrs(k,ntype)) then
             m=m+1
             index_numdis1(m)=num
             index_nwave(1,k)=index_nwave(1,k)+1
           else if(dis(num)<cenrs(k,ntype)) then
             j=j+1
             index_nwave(2,k)=index_nwave(2,k)+1
             index_numdis(j)=num
             tmp=dis(num)-initrs(k,ntype)
             tmp1=dis(num)-finalrs(k,ntype)
             tmp2=tmp*tmp1
             tmp3=tmp2*normrs(k,ntype)
             effect(j)=dexp(tmp2*tmp3)*expon(k,ntype)
             dedr(:,j)=4d0*effect(j)*tmp3*(dis(num)-cenrs(k,ntype))*drdx(:,num)
             effect(j)=effect(j)+expalpha(k,ntype)
           end if
         end do
       end do
     return
end subroutine
