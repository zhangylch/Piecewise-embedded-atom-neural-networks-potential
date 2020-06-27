subroutine deallocate_all
     use nnmod
     implicit none
       deallocate(nl)
       deallocate(nhid)
       deallocate(index_ele)
       deallocate(nwave)
       deallocate(npara)
       deallocate(factorial)
       deallocate(factor_wave)
       deallocate(index_power)
       deallocate(inv_power)
       deallocate(index_orbit)
       deallocate(inta)
       deallocate(rs)
       deallocate(w)
       deallocate(weight_wave)
       deallocate(maxwf)
       deallocate(vector)  ! allocate in readinput
       deallocate(index_rs)
       deallocate(index_numrs)
       deallocate(atom)
       deallocate(atomtype)
       deallocate(rc)
       deallocate(rcsq)
       deallocate(interaction)
     return
end subroutine
