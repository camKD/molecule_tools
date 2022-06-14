from md_to_pd_tools import dihedral_into_pd
from identify_rotatable_bonds import identify_main_rotatable_bonds

struc = 'epctcng.pdb'
traj = 'crest_rotamers.xyz'

dihedral_df = dihedral_into_pd(identify_main_rotatable_bonds(struc), struc, traj)

print(dihedral_df.tail())
