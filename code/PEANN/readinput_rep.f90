subroutine readinput_rep
     use representation
     use sharedmod
     implicit none
     integer(kind=intype) :: i,j,k,l,ntype,num
     real(kind=typenum) :: tmp,dier_rs
       open(100,file=trim(adjustl(parafile))//'/input_density')
       read(100,*) maxnumtype
       read(100,*) maxnumatom,maxneff
       read(100,*) ipsin
       read(100,*) nwave
       read(100,*) cutnum,ncell
       read(100,*) table_grid
       read(100,*) rc
       dier=rc+1d-2
       pcutnum=cutnum*nwave
       wavecutnum=pcutnum*(ipsin+1)
       totpara=0
       rcsq=rc*rc
       interaction=ceiling(rc/dier)
       do i=0,ipsin
         do j=0,i
           do k=0,i-j
             totpara=totpara+1
           end do
         end do
       end do
       atomwave=0
       do i=0,ipsin
         do k=0,nwave-1
           do ntype=1,maxnumtype
             atomwave=atomwave+1
           end do
         end do
       end do
       call allocate_all_rep
       factorial(0)=1d0
       do i=1,ipsin
         factorial(i)=factorial(i-1)*dble(i)
       end do
       npara=0 
       num=0
       do i=0,ipsin
         do j=0,i
           do k=0,i-j
             l=i-j-k
             num=num+1
             npara(i)=npara(i)+1
             factor_wave(num)=factorial(i)/(factorial(j)*factorial(k)*factorial(l)) 
           end do
         end do
       end do
       maxnpara=maxval(npara)
       norbit=nwave*(1+ipsin)
       if (table_grid==0) then
         do ntype=1,maxnumtype
           read(100,*) 
           do j=0,nwave-1
             read(100,*) initrs(j,ntype),cenrs(j,ntype),alpha
             finalrs(j,ntype)=cenrs(j,ntype)+cenrs(j,ntype)-initrs(j,ntype)
             normrs(j,ntype)=-alpha/(finalrs(j,ntype)-cenrs(j,ntype))**4
             expon(j,ntype)=1d0/(1d0-dexp(-alpha))
             expalpha(j,ntype)=-dexp(-alpha)*expon(j,ntype)
           end do
         end do
       else 
         do ntype=1,maxnumtype
           read(100,*) breadth,alpha,init_wf
           dier_rs=(rc-init_wf)/((nwave-1)*1d0+breadth)
           breadth=dier_rs*breadth
           initrs(0,ntype)=init_wf
           cenrs(0,ntype)=init_wf+breadth
           if(cenrs(0,ntype)>rc) cenrs(0,ntype)=rc
           finalrs(0,ntype)=cenrs(0,ntype)+cenrs(0,ntype)-initrs(0,ntype)
           tmp=alpha*breadth/3.5d0**initrs(0,ntype)
           normrs(0,ntype)=-tmp/(finalrs(0,ntype)-cenrs(0,ntype))**4
           expon(0,ntype)=1d0/(1d0-dexp(-tmp))
           expalpha(0,ntype)=-dexp(-tmp)*expon(0,ntype)
           do j=1,nwave-1
             initrs(j,ntype)=initrs(j-1,ntype)+dier_rs
             cenrs(j,ntype)=cenrs(j-1,ntype)+dier_rs
             if(initrs(j,ntype)<init_wf) initrs(j,ntype)=init_wf
             if(cenrs(j,ntype)>rc) cenrs(j,ntype)=rc
             finalrs(j,ntype)=cenrs(j,ntype)+cenrs(j,ntype)-initrs(j,ntype)
             tmp=alpha*breadth/3.5d0**initrs(j,ntype)
             normrs(j,ntype)=-tmp/(finalrs(j,ntype)-cenrs(j,ntype))**4
             expon(j,ntype)=1d0/(1d0-dexp(-tmp))
             expalpha(j,ntype)=-dexp(-tmp)*expon(j,ntype)
           end do
         end do
       end if
       close(100)
       do i=1,maxnumtype
         do j=1,atomwave
           call random_number(tmp)
           weight_wave(j,i)=2d0*tmp-1d0
         end do
       end do
     return
end subroutine
