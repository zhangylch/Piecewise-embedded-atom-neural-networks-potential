subroutine get_wave(natom,dire,cart,density)
     use constant
     use representation
     use sharedmod
     implicit none
     integer(kind=intype) natom,scutnum
     integer(kind=intype) index_dis(cutnum)
     real(kind=typenum) dire(atomdim,numatom,length),cart(atomdim),dis(cutnum)
     real(kind=typenum) effvec(atomdim,cutnum)
     real(kind=typenum) density(norbit)
!f2py integer(kind=intype),intent(aux) :: cutnum,norbit,numatom,length
       call period_dis(natom,dire,cart,scutnum,index_dis,dis,effvec)
       call get_density(natom,scutnum,index_dis,dis,effvec,density)
     return
end subroutine
