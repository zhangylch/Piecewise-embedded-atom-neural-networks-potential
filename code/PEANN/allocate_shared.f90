subroutine allocate_shared
     use sharedmod
     use representation
     implicit none
       allocate(index_ele(maxnumatom))
       allocate(atom(maxnumatom))
       allocate(atomtype(maxnumtype))
       allocate(nl(0:mnhid+1,maxnumtype))
       allocate(nhid(maxnumtype))
       allocate(z(0:mnl,0:mnhid+1))
       allocate(w(0:mnl,1:mnl,mnhid+1,maxnumtype,nkpoint))
       allocate(vector(3,length))
       allocate(index_rs(numrs(1),numrs(2),numrs(3)))
       allocate(index_numrs(2,ncell,numrs(1),numrs(2),numrs(3)))
     return
end subroutine 
