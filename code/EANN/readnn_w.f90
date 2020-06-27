subroutine readnn_w
     use constant
     use nnmod 
     implicit none
     integer ilayer,jneuron,jp,i,ntype,k
     character(len=80) ckpoint
       do k=1,nkpoint
         write(ckpoint,*) k
         do ntype=1,numtype
           open(123,file=trim(adjustl(parafile))//'/W_'//trim(adjustl(atomtype(ntype)))//trim(adjustl(ckpoint)))
           read(123,*)
           read(123,*)
           do ilayer=1,nhid(ntype)+1
             do jneuron=1,nl(ilayer,ntype)
               do jp=1,nl(ilayer-1,ntype)
                 read(123,*) w(jp,jneuron,ilayer,ntype,k)
               enddo
             enddo
             read(123,*) 
             do jneuron=1,nl(ilayer,ntype)
               read(123,*) w(0,jneuron,ilayer,ntype,k)
             enddo
             read(123,*)
           enddo
           close(123)
         end do
       end do
       close(123)
       do i=1,numtype
         open(123,file=trim(adjustl(parafile))//'/weight_wave_'//trim(adjustl(atomtype(i))))
         read(123,*) weight_wave(:,i)
         open(124,file=trim(adjustl(parafile))//'/scalfactor_'//trim(adjustl(atomtype(i))))
         read(124,*) maxwf(:,i)
         close(123)
         close(124)
         do jneuron=1,atomwave
           jp=index_orbit(jneuron)
           weight_wave(jneuron,i)=weight_wave(jneuron,i)/maxwf(jp,i)
         end do
       end do
     return
end subroutine
