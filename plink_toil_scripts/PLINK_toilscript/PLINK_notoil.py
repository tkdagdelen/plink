from os import listdir
from os.path import isfile, join
import argparse
from toil.job import Job
import subprocess

def concordanceCheck():
	# read in the parameters
	input_args= readParameters('Run4_Parameters.txt')
	# make the commands
	# PLINK_IBS1 = ['--bfile', input_args['bfile'], 
	# 				'--Z-genome',
	# 				'--out', 'nodocker_pairs']
	# PLINK_IBS2 = ['--bfile', 'nodocker_filtered',
	# 	'--read-genome', 'nodocker_pairs',
	# 	'--cluster',
	# 	'--K', input_args['num_clusters'],
	# 	'--out', 'nodocker_clusterfiles']
	PLINK_IBS2 = ['--bfile', input_args['bfile'],
		'--cluster',
		'--K', input_args['num_clusters'],
		'--out', 'nodocker_clusterfiles']
	PLINK_regression = ['--bfile', input_args['bfile'], '--' + input_args['analysis_type'],
		'--within', 'nodocker_clusterfiles.cluster2',
		'--covar', 'fixed_' + input_args['covar_file'],
		'--covar-name', input_args['covar_name'],
		'--pheno', 'fixed_' + input_args['pheno'],
		'--pheno-name', input_args['pheno_name'],

		'--out', 'nodocker_output',
		'--adjust']
	commands = [PLINK_IBS2] + [PLINK_regression]
	print(commands)
	
	# call the commands using subprocess.check_call()
	for command in commands: 
		print(command)
		subprocess.check_call(['plink2'] + command)
	print ('Done.')

def readParameters(parameter_file): 
	param_dict = {}
	with open(parameter_file, 'r') as f: 
		lines = f.readlines()
		grabNextLine=False
		for line in lines: 
			if grabNextLine:
				param_name = line[1:-2].strip('\n')
				param_dict[param_name] = param_value
				grabNextLine = False
				param_value = None
			if line.startswith('@parameter'):
				parts = line.split(':')
				if len(parts)>1:
					param_value = parts[1].strip(' ').strip('\n')
					grabNextLine=True
	return param_dict


concordanceCheck()

"""to run this on the cluster: 
ssh tkdagdelen@amp-bdg-master.amplab.net
source activate py27
tmux new -s ConcordanceCheck 
python toilConcordanceCheck.py 
ctrl + d + b
"""

"""
simple run with only plink
plink2 --noweb --bfile SUCCESSA_Top_subject_level_filtered --linear --pheno fixed_phs000547.v1.pht003007.v1.p1.c1.GARNET_SUCCESS_Subject_Phenotypes.BCETC.txt --pheno-name IS_AE_CYCLE123 --within clusterFilez.cluster2 --out Output --adjust
"""

