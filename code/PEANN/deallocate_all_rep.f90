subroutine deallocate_all
     use sharedmod
     use representation
     implicit none
       deallocate(cenrs)
       deallocate(normrs)
       deallocate(finalrs)
       deallocate(initrs)
       deallocate(npara)
       deallocate(factorial)
       deallocate(factor_wave)
       deallocate(index_power)
       deallocate(weight_wave)
       deallocate(inv_power)
       deallocate(index_ele)
       deallocate(atom)
       deallocate(atomtype)
       deallocate(nl)
       deallocate(nhid)
       deallocate(z)
       deallocate(w)
       deallocate(vector)
       deallocate(index_rs)
       deallocate(index_numrs)
       deallocate(expon)
       deallocate(expalpha)
     return
end subroutine
