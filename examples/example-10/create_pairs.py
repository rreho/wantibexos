s1 = 'begin nnkpts\n'

s2 = 'end nnkpts\n'


# we canot create only half of the k-point pairs, since the in the nnkpts block all kpoints must have the same number of neighbours
with open('nnkpts_list.dat','w') as f:
  f.write(s1)
  for i in range(1,217):
    for j in range(1,217):
      f.write("%i %i 0 0 0\n" %(i,j))
  f.write(s2)
