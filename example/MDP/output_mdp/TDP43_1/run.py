from calvados import sim
from argparse import ArgumentParser

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--path',nargs='?', default='.', const='.', type=str)
    parser.add_argument('--config',nargs='?', default='config.yaml', const='config.yaml', type=str)
    parser.add_argument('--components',nargs='?', default='components.yaml', const='components.yaml', type=str)

    args = parser.parse_args()

    path = args.path
    fconfig = args.config
    fcomponents = args.components

    sim.run(path=path,fconfig=fconfig,fcomponents=fcomponents)

from calvados.analysis import SlabAnalysis, calc_com_traj, calc_contact_map

slab = SlabAnalysis(name="TDP43_1", input_path="/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square/example/MDP/output_mdp/TDP43_1",
  output_path="data", ref_name="TDP43_1", verbose=True)

slab.center(start=0, center_target='all')
slab.calc_profiles()
slab.calc_concentrations()
print(slab.df_results)
slab.plot_density_profiles()

chainid_dict = dict(TDP43_1 = (0,14))
calc_com_traj(path="/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square/example/MDP/output_mdp/TDP43_1",sysname="TDP43_1",output_path="data",residues_file="/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square/example/MDP/output_mdp/input/residues_CALVADOS3.csv",chainid_dict=chainid_dict)
calc_contact_map(path="/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square/example/MDP/output_mdp/TDP43_1",sysname="TDP43_1",output_path="data",chainid_dict=chainid_dict,is_slab=True)
