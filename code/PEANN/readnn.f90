subroutine readnn
     use constant
     use representation
     use sharedmod
     implicit none
     integer(kind=intype) :: i,j,k,l,ntype
     real(kind=typenum) :: tmp(3),vec1(3,3),vec2(3,3),mind(3),maxd(3)
     character*30 a,ckpoint,b
       open(100,file=trim(adjustl(parafile))//'/input_nn')
       read(100,*)
       read(100,*) 
       read(100,*) 
       read(100,*) 
       read(100,*) 
       read(100,*) 
       read(100,*) 
       read(100,*) 
       read(100,*) 
       read(100,*) 
       read(100,*) 
       read(100,*) 
       read(100,*) 
       read(100,*) mnl
       read(100,*) mnhid
       read(100,*) nkpoint
       read(100,*) outputneuron
!-----------read the wave_nn structure---------------
       open(111,file=trim(adjustl(parafile))//'/cell')
       read(111,*) scalmatrix(:,:)
       close(111) 
       call inverse_matrix(scalmatrix(:,:),inv_matrix)
       nimage(1)=ceiling(rc/scalmatrix(1,1))
       nimage(2)=ceiling(rc/scalmatrix(2,2))
       nimage(3)=ceiling(rc/scalmatrix(3,3))
       length=(2*nimage(1)+1)*(2*nimage(2)+1)*(2*nimage(3)+1)
       mind=1d5
       maxd=-1d5
       do i=0,1
         do j=0,1
           do k=0,1
             tmp(1)=k*1d0
             tmp(2)=j*1d0
             tmp(3)=i*1d0
             tmp=matmul(scalmatrix,tmp)
             do l=1,atomdim     
               if(tmp(l)<mind(l)) mind(l)=tmp(l)
               if(tmp(l)>maxd(l)) maxd(l)=tmp(l)
             end do
           end do
         end do
       end do
       rangecoor=(maxd-mind+2d0*rc)+hinc
       numrs=ceiling(rangecoor/dier)
       call allocate_shared
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
       z(0,:)=1d0
       nl(0,:)=norbit
       do i=1,maxnumtype
         read(100,*) atomtype(i),nhid(i),nl(1:nhid(i),i)
         nl(nhid(i)+1,i)=outputneuron
       end do
       close(100)
       do k=1,nkpoint
         write(ckpoint,*) k
         do ntype=1,maxnumtype
           open(123,file=trim(adjustl(parafile))//'/W_'//trim(adjustl(atomtype(ntype)))//trim(adjustl(ckpoint)))
           read(123,*)
           read(123,*)
           do i=1,nhid(ntype)+1
             do j=1,nl(i,ntype)
               do l=1,nl(i-1,ntype)
                 read(123,*) w(l,j,i,ntype,k)
               enddo
             enddo
             read(123,*) 
             do j=1,nl(i,ntype)
               read(123,*) w(0,j,i,ntype,k)
             enddo
             read(123,*)
           enddo
           close(123)
         end do
       end do
       do i=1,maxnumtype
         open(123,file=trim(adjustl(parafile))//'/weight_wave_'//trim(adjustl(atomtype(i))))
         read(123,*) weight_wave(:,i)
         close(123)
       end do
       close(123)
       numatom=0
       neff=0
       open(100,file=trim(adjustl(parafile))//'/atom') 
1011   read(100,*,end=1010) b,a
         numatom=numatom+1
         atom(numatom)=b
         if(a=='T') neff=neff+1
         goto 1011
1010   continue
       close(100)
       numforce=neff*3
       do i=1,numatom
         do j=1,maxnumtype
           if(atom(i)==atomtype(j)) index_ele(i)=j
         end do
       end do
     return
end subroutine
