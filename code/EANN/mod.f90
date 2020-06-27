module constant
     integer(kind=4),parameter :: intype=4,typenum=8,atomdim=3
     real(kind=typenum),parameter :: pi=dacos(-1d0),hinc=1d-8
end module
module nnmod
     use constant
     integer(kind=intype) :: numatom,neff,numforce,gasatom,numtype,length,table_grid,ncell
     integer(kind=intype) :: maxnwave,norbit,ipsin,atomwave,totpara
     integer(kind=intype) :: mnl,mnhid,nkpoint,outputneuron,cutnum,pcutnum
     integer(kind=intype) :: nimage(atomdim),numrs(atomdim)
     real(kind=typenum) :: dier,maxrc
!----------------------------------real-------------------------------------------------
     real(kind=typenum) :: scalmatrix(atomdim,atomdim),inv_matrix(atomdim,atomdim),rangecoor(atomdim)
!---------------------------------integer allocatable array---------------------------------------------------
     integer(kind=intype),allocatable :: nl(:,:),nhid(:)
     integer(kind=intype),allocatable :: index_ele(:)
     integer(kind=intype),allocatable :: nwave(:),npara(:)
     integer(kind=intype),allocatable :: index_orbit(:),index_power(:,:),inv_power(:,:,:)
     integer(kind=intype),allocatable :: index_rs(:,:,:),index_numrs(:,:,:,:,:),interaction(:)
!----------------------------------real allocatable array----------------------------------------------------------------
     real(kind=typenum),allocatable :: inta(:,:),rs(:,:),rc(:),rcsq(:)
     real(kind=typenum),allocatable :: w(:,:,:,:,:),weight_wave(:,:),maxwf(:,:)
     real(kind=typenum),allocatable :: factorial(:),factor_wave(:)
     real(kind=typenum),allocatable :: vector(:,:)
!----------------------------------character allocatable array----------------------------------------------
     character(len=5),allocatable :: atom(:),atomtype(:)
     character(len=20) :: parafile='para'
end module
