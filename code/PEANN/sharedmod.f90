module sharedmod
    use constant
    implicit none
      integer(kind=intype) :: nwave,norbit,atomwave,pcutnum,wavecutnum,ncell
      integer(kind=intype) :: numforce,maxnumtype,maxnumatom,maxneff,neff,numatom,length
      integer(kind=intype) :: numrs(atomdim),interaction,nimage(atomdim)
      integer(kind=intype) :: mnl,mnhid,nkpoint,outputneuron
      integer(kind=intype),allocatable :: index_ele(:)
      integer(kind=intype),allocatable :: nl(:,:),nhid(:)
      integer(kind=intype),allocatable :: index_rs(:,:,:),index_numrs(:,:,:,:,:)
      real(kind=typenum) :: rc,rcsq
      real(kind=typenum) :: scalmatrix(atomdim,atomdim),inv_matrix(atomdim,atomdim),rangecoor(atomdim)
      real(kind=typenum),allocatable :: vector(:,:)
      real(kind=typenum),allocatable :: z(:,:),w(:,:,:,:,:)
      real(kind=typenum),allocatable :: weight_wave(:,:)
      character(len=5),allocatable :: atom(:),atomtype(:)
!-------------------------for test time--------------------------
end module
