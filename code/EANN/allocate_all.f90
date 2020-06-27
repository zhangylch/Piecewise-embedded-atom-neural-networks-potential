subroutine allocate_all
     use nnmod
     implicit none
     integer(kind=intype) :: numrs1,numrs2,numrs3
       numrs1=numrs(1)
       numrs2=numrs(2)
       numrs3=numrs(3)
       allocate(atom(numatom))
       allocate(atomtype(numtype))
       allocate(nl(0:mnhid+1,numtype))
       allocate(nhid(numtype))
       allocate(w(0:mnl,1:mnl,mnhid+1,numtype,nkpoint))
       allocate(index_ele(numatom))
       allocate(weight_wave(atomwave,numtype))
       allocate(index_orbit(atomwave))
       allocate(index_power(3,totpara))
       allocate(inv_power(-1:ipsin,-1:ipsin,-1:ipsin))
       allocate(factor_wave(totpara))
       allocate(maxwf(norbit,numtype))
       allocate(index_rs(numrs1,numrs2,numrs3))
       allocate(index_numrs(2,ncell,numrs1,numrs2,numrs3))
     return
end subroutine
