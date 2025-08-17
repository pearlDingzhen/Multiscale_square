#!/bin/bash
#SBATCH --job-name=TDP43_droplet
#SBATCH --partition=MD4G-cluster
#SBATCH --nodes=1
#SBATCH --nodelist=MD4G01
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2048
#SBATCH --gres=gpu:1



source activate calvados
conda activate calvados

replica=1

python prepare.py --name TDP43 --replica $replica --temp 240 --number 1 --box 15
python TDP43_$replica/run.py --path TDP43_$replica


