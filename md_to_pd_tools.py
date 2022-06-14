import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Dihedral
import pandas as pd

def dihedral_into_pd(dihedrals, struc, traj):
    # Load the universe with structure and trajectory
    u = mda.Universe(struc, traj)
    # Create an empty dataframe
    dihedral_df = pd.DataFrame()
    counter = 1
    # Iterate through the list of dihedral atom lists
    for d in dihedrals:
        # Make a selection of the 4 atoms involved in the current dihedral
        latest_d = conformers.atoms[d]
        # Create an empty list to store the each dihedrals values
        dis = []
        # Iterate through the trajectory
        for ts in u.trajectory:
            dis.append(latest_d.dihedral.value())
        d_label = 'dihedral_'+str(counter)
        dihedral_df[d_label] = dis
        counter += 1
    return dihedral_df

