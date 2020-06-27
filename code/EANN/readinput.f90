subroutine readinput
     use constant
     use nnmod 
     implicit none
     integer(kind=intype) :: i,j,k,l,ntype,num
     real(kind=typenum) :: hyperpara,dier_rs,tmp(3),vec1(3,3),vec2(3,3),mind(3),maxd(3)
     character*10 a
       open(100,file=trim(adjustl(parafile))//'/input')
       read(100,*) 
       read(100,*) 
       read(100,*) 
       read(100,*)
       read(100,*) table_grid
       read(100,*) 
       read(100,*) numtype
       read(100,*) numatom,neff
       read(100,*) 
       read(100,*) ipsin
       read(100,*) maxnwave
       call allocate_wave
       read(100,*) nwave
       read(100,*) cutnum,ncell
       read(100,*) 
       read(100,*) 
       read(100,*) 
       read(100,*) 
       read(100,*) 
       read(100,*) 
       read(100,*) rc
       read(100,*) mnl
       read(100,*) mnhid
       read(100,*) nkpoint
       read(100,*) outputneuron
!-----------read the wave_nn structure----------------
       open(300,file=trim(adjustl(parafile))//'/cell')
       read(300,*) scalmatrix
       close(300)
       maxrc=maxval(rc)
       dier=maxrc+1d-2
       rcsq=rc*rc
       interaction=ceiling(rc/dier)
       call inverse_matrix(scalmatrix,inv_matrix)
       tmp(1)=scalmatrix(1,1)
       tmp(2)=scalmatrix(2,2)
       tmp(3)=scalmatrix(3,3)
       nimage=ceiling(maxrc/tmp)
       length=(2*nimage(1)+1)*(2*nimage(2)+1)*(2*nimage(3)+1)
       allocate(vector(atomdim,length))
       mind=1d5
       maxd=-1d5
       do i=0,1
         do j=0,1
           do k=0,1
             tmp(1)=(k*1)*1d0
             tmp(2)=(j*1)*1d0
             tmp(3)=(i*1)*1d0
             tmp=matmul(scalmatrix,tmp)
             do l=1,atomdim
               if(tmp(l)<mind(l)) mind(l)=tmp(l)
               if(tmp(l)>maxd(l)) maxd(l)=tmp(l)
             end do
           end do
         end do
       end do
       rangecoor=maxd-mind+2d0*maxrc+hinc
       numrs=ceiling(rangecoor/dier)
       vec2(:,1)=-scalmatrix(:,1)*nimage(1)
       vec2(:,2)=-scalmatrix(:,2)*nimage(2)
       vec2(:,3)=-scalmatrix(:,3)*nimage(3)
       l=0
       vec1(:,3)=vec2(:,3)
       do i=-nimage(3),nimage(3)
         vec1(:,2)=vec2(:,2)
         do j=-nimage(2),nimage(2)
           vec1(:,1)=vec2(:,1)
           do k=-nimage(1),nimage(1)
             l=l+1
             vector(:,l)=vec1(:,1)+vec1(:,2)+vec1(:,3)
             vec1(:,1)=vec1(:,1)+scalmatrix(:,1)
           end do
           vec1(:,2)=vec1(:,2)+scalmatrix(:,2)
         end do
         vec1(:,3)=vec1(:,3)+scalmatrix(:,3)
       end do
       norbit=(maxnwave+1)*(ipsin+1)
       pcutnum=norbit*cutnum
       atomwave=0
       do i=0,ipsin
         do k=0,maxnwave
           do ntype=1,numtype
             atomwave=atomwave+1
           end do
         end do
       end do
       factorial(0)=1d0
       do i=1,ipsin
         factorial(i)=factorial(i-1)*dble(i)
       end do
       totpara=0
       do i=0,ipsin
         do j=0,i
           do k=0,i-j
             totpara=totpara+1
           end do
         end do
       end do
       call allocate_all
       npara=0 
       num=0
       do i=0,ipsin
         do j=0,i
           do k=0,i-j
             num=num+1
             npara(i)=npara(i)+1
             l=i-j-k
             factor_wave(num)=factorial(i)/(factorial(j)*factorial(k)*factorial(l)) 
           end do
         end do
       end do
       nl(0,:)=norbit
       if(table_grid==0) then
         do i=1,numtype
           read(100,*) atomtype(i),nhid(i),nl(1:nhid(i),i)
           nl(nhid(i)+1,i)=outputneuron
           do j=0,nwave(i) 
             read(100,*) l,inta(i,j),rs(i,j)
           end do
         end do
       else
         rs=0d0
         do i=1,numtype
           dier_rs=rc(i)/(nwave(i)+1d0/3d0)
           read(100,*) atomtype(i),nhid(i),nl(1:nhid(i),i),hyperpara
           nl(nhid(i)+1,i)=outputneuron
           inta(i,:)=hyperpara/dier_rs**2
           do j=1,nwave(i)
             rs(i,j)=rs(i,j-1)+dier_rs
           end do
         end do
       end if
       close(100)
       numforce=neff*atomdim
       open(100,file=trim(adjustl(parafile))//'/atom')
       do i=1,numatom
         read(100,*) atom(i),a
       end do
       close(100)
     return
end subroutine
