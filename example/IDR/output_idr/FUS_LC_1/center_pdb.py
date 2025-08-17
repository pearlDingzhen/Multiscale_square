import MDAnalysis as mda 
u = mda.Universe ('FUS_LC_1_20251708_00h39m56s.pdb')

with mda.Writer('center_pdb.pdb', len(u.atoms)) as W:
	u.atoms.translate(-u.atoms.center_of_geometry()+0.5*u.dimensions[:3])
	W.write(u.atoms)