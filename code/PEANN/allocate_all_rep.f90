subroutine allocate_all_rep
     use representation
     use sharedmod
     implicit none
       allocate(cenrs(0:nwave-1,maxnumtype))
       allocate(normrs(0:nwave-1,maxnumtype))
       allocate(finalrs(0:nwave-1,maxnumtype))
       allocate(initrs(0:nwave-1,maxnumtype))
       allocate(npara(0:ipsin))
       allocate(factorial(0:ipsin))
       allocate(factor_wave(totpara))
       allocate(index_power(3,totpara))
       allocate(weight_wave(atomwave,maxnumtype))
       allocate(inv_power(-1:ipsin,-1:ipsin,-1:ipsin))
       allocate(expon(0:nwave-1,maxnumtype))
       allocate(expalpha(0:nwave-1,maxnumtype))
     return
end subroutine
