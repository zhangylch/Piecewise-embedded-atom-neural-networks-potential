module representation
    use constant
    implicit none
      integer(kind=intype) :: totpara,table_grid
      integer(kind=intype) :: ipsin,maxnpara
      integer(kind=intype) :: cutnum
!---------------------------------------------------------------------------------------------
      integer(kind=intype),allocatable :: npara(:)
      integer(kind=intype),allocatable :: index_power(:,:),inv_power(:,:,:)
!-------------------------------------------------------------------------------------------
      real(kind=typenum),allocatable :: finalrs(:,:),initrs(:,:),cenrs(:,:),normrs(:,:)
      real(kind=typenum),allocatable :: factorial(:),factor_wave(:)
      real(kind=typenum),allocatable :: expon(:,:),expalpha(:,:)
      real(kind=typenum) :: init_wf,breadth,dier,alpha
      character(len=30) :: parafile='para'
end module
