from openmm.unit import *
from openmm import *
from openmm.app import *
import PACE2openmm as PACE
from mdtraj.reporters import XTCReporter
from sys import stdout
import time


def add_ring_bond_to_system(system, topology, ring_distance_dict=None):
	"""
	为芳香族氨基酸的环结构添加键约束，以保持环的几何形状
	
	Parameters
	----------
	system : openmm.System
		OpenMM系统对象
	topology : openmm.app.Topology
		拓扑对象
	ring_distance_dict : dict, optional
		环键距离字典，如果为None则使用默认值
		
	Returns
	-------
	ring_constraint_force : openmm.CustomBondForce
		添加的环约束力对象
	"""
	if ring_distance_dict is None:
		ring_distance_dict = {'HIS': {'ring': {'CD2-CE1': 2.1841,
						'CD2-ND1': 2.1878,
						'CE1-CG': 2.191,
						'CG-NE2': 2.1907,
						'ND1-NE2': 2.1715}},
		'PHE': {'ring': {'CD1-CD2': 2.407,
						'CD1-CE2': 2.7805,
						'CD1-CZ': 2.408,
						'CD2-CE1': 2.7794,
						'CD2-CZ': 2.4073,
						'CE1-CE2': 2.4079,
						'CE1-CG': 2.4075,
						'CE2-CG': 2.4072,
						'CG-CZ': 2.78}},
		'TRP': {'five_membered': {'CD1-CD2': 2.2398,
								'CD1-CE2': 2.258,
								'CD2-NE1': 2.2095,
								'CE2-CG': 2.2965,
								'CG-NE1': 2.2316},
				'six_membered': {'CD2-CH2': 2.7998,
								'CD2-CZ2': 2.425,
								'CD2-CZ3': 2.425,
								'CE2-CE3': 2.3957,
								'CE2-CH2': 2.4244,
								'CE2-CZ3': 2.75,
								'CE3-CH2': 2.3675,
								'CE3-CZ2': 2.7504,
								'CZ2-CZ3': 2.3393}},
		'TYR': {'ring': {'CD1-CD2': 2.3852,
						'CD1-CE2': 2.742,
						'CD1-CZ': 2.4075,
						'CD2-CE1': 2.7418,
						'CD2-CZ': 2.3642,
						'CE1-CE2': 2.3416,
						'CE1-CG': 2.4078,
						'CE2-CG': 2.4077,
						'CG-CZ': 2.7802}}}
	
	# 创建环约束力
	ring_constraint_force = CustomBondForce('k_ring * (r - r0)^2')
	ring_constraint_force.addPerBondParameter('k_ring')
	ring_constraint_force.addPerBondParameter('r0')
	
	# 环约束的力常数 (kJ/mol/nm^2)
	k_ring = 1000.0  # 较强的约束力常数
	
	# 创建原子名称到索引的映射
	atom_name_to_index = {}
	for atom in topology.atoms():
		atom_name_to_index[atom.name] = atom.index
	
	# 为每个芳香族氨基酸添加环约束
	ring_constraints_added = 0
	
	for residue in topology.residues():
		residue_name = residue.name
		
		if residue_name in ring_distance_dict:
			print(f"为 {residue_name} 残基添加环约束...")
			
			# 获取该残基的环结构信息
			ring_info = ring_distance_dict[residue_name]
			
			# 处理不同类型的环结构
			for ring_type, bonds in ring_info.items():
				print(f"  添加 {ring_type} 环约束...")
				
				for bond_name, distance in bonds.items():
					atom1_name, atom2_name = bond_name.split('-')
					
					# 在当前残基中查找原子
					atom1_index = None
					atom2_index = None
					
					for atom in residue.atoms():
						if atom.name == atom1_name:
							atom1_index = atom.index
						elif atom.name == atom2_name:
							atom2_index = atom.index
					
					# 如果找到两个原子，添加约束
					if atom1_index is not None and atom2_index is not None:
						ring_constraint_force.addBond(
							atom1_index, 
							atom2_index, 
							[k_ring, distance * 0.1]  # 距离转换为nm
						)
						ring_constraints_added += 1
						print(f"    添加约束: {atom1_name}-{atom2_name} = {distance} Å")
					else:
						print(f"    警告: 未找到原子 {atom1_name} 或 {atom2_name} 在 {residue_name} 残基中")
	
	# 将环约束力添加到系统
	if ring_constraints_added > 0:
		system.addForce(ring_constraint_force)
		print(f"总共添加了 {ring_constraints_added} 个环约束")
	else:
		print("未找到需要添加环约束的芳香族氨基酸")
	
	return ring_constraint_force


def add_ring_dihedral_to_system(system, topology, k_force=-20.0, ring_dihedral_dict=None):
	"""
	为芳香族氨基酸的环结构添加二面角约束，以保持环的平面性
	
	Parameters
	----------
	system : openmm.System
		OpenMM系统对象
	topology : openmm.app.Topology
		拓扑对象
	k_force : float, default=20.0
		力常数 (kJ/mol)
	ring_dihedral_dict : dict, optional
		环二面角字典，如果为None则使用默认值
		
	Returns
	-------
	ring_dihedral_force : openmm.PeriodicTorsionForce
		添加的环二面角约束力对象
	"""
	added_dihedral_info = []
	if ring_dihedral_dict is None:
		ring_dihedral_dict = {'HIS': {'CD2-CG-ND1-CE1': 0,
		         'CE1-NE2-CD2-CG': 0,
		         'CG-ND1-CE1-NE2': 0,
		         'ND1-CE1-NE2-CD2': 0,
		         'NE2-CD2-CG-ND1': 0},
		 'PHE': {'CD1-CE1-CZ-CE2': 0,
		         'CD2-CG-CD1-CE1': 0,
		         'CE1-CZ-CE2-CD2': 0,
		         'CE2-CD2-CG-CD1': 0,
		         'CG-CD1-CE1-CZ': 0,
		         'CZ-CE2-CD2-CG': 0},
		 'TRP': {'CD1-CG-CD2-CE3': 180,
		         'CD1-NE1-CE2-CD2': 0,
		         'CD1-NE1-CE2-CZ2': 180,
		         'CD2-CE2-CZ2-CH2': 0,
		         'CD2-CG-CD1-NE1': 0,
		         'CE2-CD2-CG-CD1': 0,
		         'CE2-CZ2-CH2-CZ3': 0,
		         'CE3-CD2-CE2-CZ2': 0,
		         'CG-CD1-NE1-CE2': 0,
		         'CG-CD2-CE3-CZ3': 180,
		         'CH2-CZ3-CE3-CD2': 0,
		         'CZ2-CE2-CD2-CG': 180,
		         'CZ2-CH2-CZ3-CE3': 0,
		         'CZ3-CE3-CD2-CE2': 0,
		         'NE1-CE2-CD2-CE3': 180,
		         'NE1-CE2-CD2-CG': 0,
		         'NE1-CE2-CZ2-CH2': 180},
		 'TYR': {'CD1-CE1-CZ-CE2': 0,
		         'CD2-CG-CD1-CE1': 0,
		         'CE1-CZ-CE2-CD2': 0,
		         'CE2-CD2-CG-CD1': 0,
		         'CG-CD1-CE1-CZ': 0,
		         'CZ-CE2-CD2-CG': 0}}
	
	# 创建环二面角约束力，使用PeriodicTorsionForce
	ring_dihedral_force = PeriodicTorsionForce()
	
	# 角度转换常数
	degToRad = 3.14159265359 / 180.0
	
	# 为每个芳香族氨基酸添加环二面角约束
	ring_dihedrals_added = 0
	
	for residue in topology.residues():
		residue_name = residue.name
		
		if residue_name in ring_dihedral_dict:
			print(f"为 {residue_name} 残基添加环二面角约束...")
			
			# 获取该残基的环二面角信息
			dihedral_info = ring_dihedral_dict[residue_name]
			
			# 处理每个二面角约束
			for dihedral_name, target_angle in dihedral_info.items():
				atom_names = dihedral_name.split('-')
				
				if len(atom_names) != 4:
					print(f"    警告: 二面角 {dihedral_name} 格式不正确，需要4个原子")
					continue
				
				# 在当前残基中查找原子
				atom_indices = []
				atoms_found = True
				
				for atom_name in atom_names:
					atom_index = None
					for atom in residue.atoms():
						if atom.name == atom_name:
							atom_index = atom.index
							break
					
					if atom_index is not None:
						atom_indices.append(atom_index)
					else:
						print(f"    警告: 未找到原子 {atom_name} 在 {residue_name} 残基中")
						atoms_found = False
						break
				
				# 如果找到所有四个原子，添加约束
				if atoms_found and len(atom_indices) == 4:
					# 使用PeriodicTorsionForce的格式：periodicity, phase, k
					# periodicity=1, phase=target_angle*degToRad, k=k_force
					ring_dihedral_force.addTorsion(
						atom_indices[0], 
						atom_indices[1], 
						atom_indices[2], 
						atom_indices[3],
						1,  # 多重度
						target_angle * degToRad,  # 平衡位置（弧度）
						k_force  # 力常数
					)
					added_dihedral_info.append([atom_indices[0], atom_indices[1], atom_indices[2], atom_indices[3], target_angle, k_force])
					ring_dihedrals_added += 1
					print(f"    添加二面角约束: {dihedral_name} = {target_angle}°, k = {k_force} kJ/mol")
				else:
					print(f"    警告: 无法为 {dihedral_name} 添加约束，缺少原子")
	
	# 将环二面角约束力添加到系统
	if ring_dihedrals_added > 0:
		system.addForce(ring_dihedral_force)
		print(f"总共添加了 {ring_dihedrals_added} 个环二面角约束")
	else:
		print("未找到需要添加环二面角约束的芳香族氨基酸")
	
	return ring_dihedral_force



conf_name = 'conf.gro'
top_name = 'PACE.top'
time_step = 0.004
ref_t = 325
ref_p = 1

GMXLIB = '/mnt/hdd1/tianxj_out/gmx2022.6.cuda/share/gromacs/top/'
start_from_checkpoint = False


properties = {'Precision': 'double'}
platform = Platform.getPlatformByName("CUDA")


epsilon_r = 15
conf = GromacsGroFile(conf_name)
gro_postion = conf.getPositions()
box_vectors = conf.getPeriodicBoxVectors()


top = PACE.PACETopFile(
	top_name,
	periodicBoxVectors=box_vectors,
	defines={},
	epsilon_r=epsilon_r,
	#includeDir=GMXLIB
)
### minimization ###
### minimize the system without nonbonded interactions ###
system = top.create_system(nonbonded_cutoff=1.1* nanometer, add_nonbonded_force=False)
forces = system.getForces()

# 添加环约束
#print("添加芳香族氨基酸环约束...")
#ring_constraint_force = add_ring_bond_to_system(system, top.topology)

# 添加环二面角约束
print("添加芳香族氨基酸环二面角约束...")
ring_dihedral_force = add_ring_dihedral_to_system(system, top.topology, k_force=20.0)


integrator = LangevinIntegrator(ref_t * kelvin,
								10.0 / picosecond,
								time_step * picosecond)
integrator.setRandomNumberSeed(0)

### add CA restraint ###
restraint = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
system.addForce(restraint)
restraint.addGlobalParameter('k', 100.0*kilojoules_per_mole/nanometer)
restraint.addPerParticleParameter('x0')
restraint.addPerParticleParameter('y0')
restraint.addPerParticleParameter('z0')


for atom in top.topology.atoms():

	if atom.name == 'CA':
		restraint.addParticle(atom.index, gro_postion[atom.index])

simulation = Simulation(top.topology, system, integrator,
						platform)

simulation.context.setPositions(conf.getPositions())

################################################################################
	### Minimization ###

### energy before minimization ###
print ('starting minimize energy , add CA restraint 100 kcal/mol/nm^2')
print ('minimize energy without nonbonded interactions')
energies = simulation.context.getState(getEnergy=True).getPotentialEnergy()
print ('energy(without nonbond) before minimization', energies)





simulation.minimizeEnergy(maxIterations=5000,tolerance=100)

energies = simulation.context.getState(getEnergy=True).getPotentialEnergy()
state_bond = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
PDBFile.writeFile(simulation.topology, state_bond.getPositions(), open('optibond.pdb', 'w'))

print ('energy(without nonbond) after minimization', energies, '\n')
print ('minimize structure is saved in optibond.pdb')


system1 = top.create_system(nonbonded_cutoff=1.1* nanometer, add_nonbonded_force=True)

# 为第二个系统也添加环约束
#print("为包含非键相互作用的系统添加芳香族氨基酸环约束...")
#ring_constraint_force1 = add_ring_bond_to_system(system1, top.topology)

# 为第二个系统也添加环二面角约束
print("为包含非键相互作用的系统添加芳香族氨基酸环二面角约束...")
ring_dihedral_force1 = add_ring_dihedral_to_system(system1, top.topology, k_force=20.0)

restraint1 = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
system1.addForce(restraint1)
restraint1.addGlobalParameter('k', 100*kilojoules_per_mole/nanometer)
restraint1.addPerParticleParameter('x0')
restraint1.addPerParticleParameter('y0')
restraint1.addPerParticleParameter('z0')

for atom in top.topology.atoms():
	
	if atom.name == 'CA':
		restraint1.addParticle(atom.index, gro_postion[atom.index])

integrator = LangevinIntegrator(ref_t * kelvin,
								10.0 / picosecond,
								time_step * picosecond)
integrator.setRandomNumberSeed(0)

#nbforce = system1.getForces()[0]
#nbforce.getGlobalParameterName(0)
#nbforce.setGlobalParameterDefaultValue(0, 0.001)

simulation1 = Simulation(top.topology, system1, integrator,
						platform, properties)

simulation1.context.setPositions(state_bond.getPositions())

state = simulation1.context.getState(getEnergy=True, getForces=False)

print ('energy(with nonbond) before minimization', state.getPotentialEnergy(), '\n')

simulation1.minimizeEnergy(maxIterations=5000,tolerance=200)

#state = simulation1.context.getState(getEnergy=True, getForces=False)

#print ('energy(with nonbond) after minimization', state.getPotentialEnergy(), '\n')

state_nonbond = simulation1.context.getState(getPositions=True, enforcePeriodicBox=True)
PDBFile.writeFile(simulation1.topology, state_nonbond.getPositions(), open('optinonbond.pdb', 'w'))




state = simulation1.context.getState(getEnergy=True, getForces=False)

print ('energy(with bond) after minimization', state.getPotentialEnergy(), '\n')
# in green
print ('\033[32mminimize success, saved in optinonbond.pdb  \033[0m')
print ('')

################################################################################
	### Minimization ###
'''
simulation.reporters.append(StateDataReporter(stdout, nstlog,
												step=True,
												potentialEnergy=True,
												temperature=True,
												volume=True,
												speed=True)
							)
'''
'''
print("Minimizing energy...")
energies = simulation.context.getState(getEnergy=True).getPotentialEnergy()
print ('energy after bond minimization', energies, '\n')
state = simulation.context.getState(getEnergy=True, getForces=False)



simulation1.minimizeEnergy(maxIterations=5000,tolerance=100)
print ('energy after nonbond minimization', energies, '\n')

energies = simulation.context.getState(getEnergy=True).getPotentialEnergy()
state_bond = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)

PDBFile.writeFile(simulation.topology, state_bond.getPositions(), open('optinonbond.pdb', 'w'))
'''
'''
state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
PDBFile.writeFile(simulation.topology, state.getPositions(), open('md.pdb', 'w'))




#PDBFile.writeFile(simulation.topology, state.getPositions(), open('openmm.pdb', 'w'))

################################################################################
	### NVT equilibration ###

simulation.context.setVelocitiesToTemperature(ref_t * kelvin)
print('Running NVT equilibration...')
simulation.step(nsteps_pr/2) #1ns

################################################################################
	### NPT equilibration ###

system.addForce(MonteCarloBarostat(ref_p * bar, ref_t * kelvin))
# to update the simulation object to take in account the new system
simulation.context.reinitialize(True)
print('Running NPT equilibration...')
simulation.step(nsteps_pr/2) #1ns

# save the equilibration results to file
simulation.saveState('equi.state')
simulation.saveCheckpoint('equi.chk')

################################################################################
	### Production run ###

# Set up the reporters to report energies every 1000 steps.
time_start = time.time()
simulation.reporters.append(StateDataReporter("md.log", nstlog,
												step=True,
												potentialEnergy=True,
												totalEnergy=True,
												density=True,
												temperature=True,
												volume=True,
												progress=True,
												remainingTime=True,
												speed=True,
												totalSteps=nsteps + nsteps_pr)
							)
# save the trajectory in XTC format
xtc_reporter = XTCReporter('md.xtc', nstxtcout)
simulation.reporters.append(xtc_reporter)

# run simulation
simulation.reporters.append(CheckpointReporter('md.chk',10000))
print("Running simulation...")
simulation.step(nsteps) #
length_of_simulation = time_step * nsteps / 1000
time_end = time.time()
time_cost = time_end - time_start
perfermance = (length_of_simulation / time_cost ) * 3600 * 24
print ('time cost: {:.2f} second'.format(time_cost))
print ('perfermance: {:.2f} ns/day'.format(perfermance))	
state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
PDBFile.writeFile(simulation.topology, state.getPositions(), open('md.pdb', 'w'))
'''	

