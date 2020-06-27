subroutine allocate_wave
     use nnmod
     implicit none
       allocate(nwave(numtype))
       allocate(inta(numtype,0:maxnwave))
       allocate(rs(numtype,0:maxnwave))
       allocate(npara(0:ipsin))
       allocate(factorial(0:ipsin))
       allocate(rc(numtype))
       allocate(rcsq(numtype))
       allocate(interaction(numtype))
     return
end subroutine
