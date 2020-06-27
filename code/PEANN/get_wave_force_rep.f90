subroutine get_wave_force(natom,dire,cart,index_natom,index_tnwave,density,pdensity)
     use constant
     use representation
     use sharedmod
     implicit none
     integer(kind=intype) natom,scutnum,ntype
     integer(kind=intype) index_natom(2,wavecutnum),index_dis(cutnum),index_tnwave(0:nwave-1) 
     integer(kind=intype) index_nwave(2,0:nwave-1),index_numdis(pcutnum),index_numdis1(pcutnum)
     real(kind=typenum) dire(atomdim,numatom,length),cart(atomdim)
     real(kind=typenum) effvec(atomdim,cutnum),dis(cutnum),effect(pcutnum),dedr(atomdim,pcutnum)
     real(kind=typenum) density(norbit)
     real(kind=typenum) pdensity(atomdim,wavecutnum)
!f2py integer(kind=intype),intent(aux) :: cutnum,totpara,maxnumtype,nwave,norbit,pcutnum,wavecutnum,numatom,length
       ntype=index_ele(natom)
       call period_dis(natom,dire,cart,scutnum,index_dis,dis,effvec)
       call der_cutoff(scutnum,ntype,index_nwave,index_numdis,index_numdis1,effvec,dis,effect,dedr)
       call get_pdensity(natom,scutnum,index_dis,index_nwave,index_numdis,index_numdis1,effect, &
       dedr,effvec,index_natom,index_tnwave,density,pdensity)
     return
end subroutine
