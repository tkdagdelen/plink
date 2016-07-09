#!/usr/bin/env python2.7
# Taner Dagdelen
# Feb. 22, 2016
# Aknowledgements: 

from os import listdir
from os.path import isfile, join
import argparse
import os
from toil.job import Job
from subprocess import call
import subprocess
from tempfile import mkstemp
from shutil import move
from os import remove, close
import csv
import sys
import urllib
import pygame

"""
Genome Wide Association Study pipeline with MAF and missingness thresholding, IBS_clustering for stratification correction, association analysis using linear/logistic regression, bonferroni correction, and sorting by P-value, using PLINK. 

Any version of PLINK can be used, and dockerfiles for both PLINK and PLINK2 can be found in the git repository:
PLINK: https://github.com/tkdagdelen/plink/blob/master/plinkdockers/PlinkDockerfile
PLINK2: https://github.com/tkdagdelen/plink/blob/master/plinkdockers/Plink2Dockerfile 


		0 --> 1 --> 2 --> 3

 0  Data massaging (mostly cutting off headers and changing ID lables to accomodate PLINK's formatting requirements)
 1  MAF and Missingness Filtering
 2  Clustering
 3  Regression Analysis (including clusters as covariates) & Results Correction (bonferroni)

Dependencies: 
Docker  -   apt-get install docker.io
Toil    -   pip install toil
PLINK   -   

Minimum data requirements: 
	- the bed/fam/bim files 
	- a phenotype file that follows convention found here: http://pngu.mgh.harvard.edu/~purcell/plink/tutorial.shtml
	- a parameter file that follows the template found here: [INSERT GITHUB LINK HERE]
 """

################# Utility functions #####################

def build_parser():
	"""
	Builds the parser to specify the parameter file
	"""
    parser = argparse.ArgumentParser(description="main.__doc__", add_help=True)
    parser.add_argument('-f', '--file', required=True, default = 'Run3_parameters.txt', help='Parameter file. See README')
    return parser

def fixPhenoFile(filename): 
	"""
	Fixes the dbGaP-standard phenotype file to match the formatting requirements of PLINK. Specifically it 
	changes the labeling convention of the phenotype file to match that of the genotype files, and cuts off 
	the headers etc. 
	"""
	fh, output_file = mkstemp()
	with open(filename) as to_read:
		with open(output_file, "wb") as tmp_file:
			reader = csv.reader(to_read, delimiter = "\t")
			writer = csv.writer(tmp_file, delimiter = " ")
			isTrue = False
			isfirst = True
			for row in reader:     # read one row at a time
				if len(row) >= 2:
					if not row[0].startswith('#'):#'dbGaP SubjID': 
						isTrue = True
				if isTrue and not isfirst:     
					myColumn = [row[1], row[1]] + list(row[i] for i in range(2,len(row)))  
					writer.writerow(myColumn) # write it
				if isfirst and isTrue: 
					myColumn = ['FID', 'IID'] + list(row[i] for i in range(2,len(row)))
					writer.writerow(myColumn)
					isfirst = False	
	close(fh)
	move(output_file,"fixed_"+filename)

def readParameters(parameter_file): 
	"""
	Reads in the parameters from the parameter_file.
	"""
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
				if param_value is not None:
					grabNextLine=True
	return param_dict

def fetchContainerFiles(filelist):
	"""
	When you mount a local directory, all the files from the local directory are brought into the Docker containter, 
	but when you write to files or create new files, those changes are not brought back into the local directory by default. 
	Instead you have to manually copy the files back to the local directory, so that they are available for use in 
	downstream docker containers/plink calls. This function manually copies the files in filelist in the docker container back into the
	local directory so that any changes that occur in the container are captured. 

	Another way to do it would be to create and mount a Data Volumes container (as described here: https://docs.docker.com/engine/userguide/containers/dockervolumes/)
	but since you'll eventually have to fetch the files from that container, this works fine for now. 
	"""
	work_dir = os.getcwd()
	for fname in filelist: 
		docker_call = 'docker cp plinkContainer:{0} {1}'.format(fname, work_dir)
		subprocess.check_call(docker_call.split())
	subprocess.check_call('docker rm plinkContainer'.split())

def docker_call(work_dir, tool_parameters = [], tool = [], outfile=None, sudo=False):
    """
    Makes subprocess call to spawn docker container, mounting the current local directory.
    tool_parameters: list   An array of the parameters to be passed to the tool
    tool: str               Name of the Docker image to be used 
    outfile: file           Filehandle that stderr will be passed to
    sudo: bool              If the user wants the docker command executed as sudo
    """
    # subprocess.check_call('docker-machine env default'.split())
    # subprocess.check_call('eval',"$(docker-machine env default)")
    work_dir = os.getcwd()
    base_docker_call = 'docker run --name plinkContainer -v {0}:/data plink2 {1}'.format(work_dir,tool).split() + tool_parameters
    try:
        subprocess.check_call(base_docker_call)
    except subprocess.CalledProcessError:
        raise RuntimeError('docker command returned a non-zero exit status. Check error logs.')
    except OSError:
        raise RuntimeError('docker not found on system. Install on all nodes.')

######################## Jobs ###########################

# Start the pipeline 
def start_pipeline(job, input_args):
	"""
	Launches the TOIL pipeline. 
	"""   
	print("Pipeline launched...") 
	job.addFollowOnJobFn(apply_filters,input_args)


def apply_filters(job, input_args): 
	"""
	This module applies MINOR Allele Frequency and missingness filters to the data and outputs a new file. 

	Equivalent PLINK command line arguments:

		>> plink --bfile [input_filename] --maf [mafthreshold] --mind [SNPmissingnessThreshold] --geno [genoRate] --out [output_filename_prefix]
					 0                       1                  2                                 3                 4     

	Operations performed (see equiv. command line argument above): 
		0   Loads the binary MAP/PED/etc. PLINK files with the filename [input_filename]
		1   Filters out SNPs with MINOR Allele Frequency below [mafthresh].
		2   Filters out sampels in which the percentage of missing SNPs is greater than [SNPmiss]
		3   Filters out SNPs for which the genotyping rate is less than [genoRate]
		4   Outputs the files with the prefix [output_filename] with the filtered data.  

	Parameters:
			 --param--                           --type--        --description-- 
		:param job:                              Job instance
		:param MapPed_filename:                  string
		:param MAF_cutoff:                       float
		:param percent_SNPs_missing_cutoff:      float
		:param genoRate_cutoff:                  float
		:param output_prefix:                    string

	Filters:

		--maf    Filters out SNPs with MINOR Allele Frequency below the given threshold. Standard MAF threshold used in 
					GWAS is MAF<.01 
		--mind   Filters out samples in which the percentage of missing SNPs is greater than SNPmissingnessThreshold. Standard 
					according to PLINK website is 10% (.1), which is also the default value. 
		--geno   Filters out SNPs for which the genotyping rate is less than rateThresh. Standard genotyping rate cutoff according to 
					PLINK website is 90%. Default is to include all SNPs. 

	For more information: http://pngu.mgh.harvard.edu/~purcell/plink/thresh.shtml
	"""
	print("Filtering and cleaning samples...")
	# if the fixed version of the phenotype file is present then use that, otherwise fix the phenotype file: 
	fname = input_args['pheno']
	if not isfile('fixed_' + fname): 
		fixPhenoFile(fname)
		input_args['pheno'] = 'fixed_' + fname
	if 'filter' in input_args and input_args['filter'] == 'yes':
		work_dir = job.fileStore.getLocalTempDir()
		lst = ['maf', 'mind', 'geno']
		if 'maf' in input_args: 
			maf = input_args['maf']
		else: 
			maf = '0.01'
		if 'mind' in input_args: 
			mind = input_args['mind']
		else:
			mind = '.1'

		if 'geno' in input_args:
			geno = input_args['geno']
		else: 
			geno = '1'

		# Put together PLINK arguments
		PLINK_command = ['--bfile', 'data/' + input_args['bfile'], # going to need to do --make-bed first!
				   '--maf', maf,
				   '--mind', mind,
				   '--geno', geno,
				   '--out', input_args['filtered_out']]
		
		# Make the docker call with the given command line arguments 
		try: 
			docker_call(work_dir,
				tool_parameters = PLINK_command, 
				tool = 'plink',
				sudo = False)
		except: 
			sys.stderr.write('Running the GWAS pipeline with %s in %s failed.' %(" ".join(PLINK_command), work_dir))
			raise

		# Link to the next job in the pipeline
		output_prefix = input_args['filtered_out']
		fetchContainerFiles([output_prefix + '.bed',output_prefix + '.fam', output_prefix + '.bim', output_prefix + '.log'])
		input_args['bfile'] = input_args['filtered_out']
	job.addFollowOnJobFn(IBS_cluster, input_args)

def IBS_cluster(job, input_args):
	"""
		This module performs 2 primary functions: 
			- Computes the pairwise identity by state for all resulting samples.
			- Clusters Samples together based on their IBS similarity 

		Equivalent PLINK command line arguments: 

			If .genome file not found in current directory

				>> plink --bfile [input_filename] --Z-genome --out [IBS_pairsfile_prefix]
							0                       1            2                  
				>> plink --bfile [input_filename] --read-genome [output_filename].genome --cluster --K [num_clusters] --out [IBS_clusterfilename]
							3                        4                                      5       6                   7

			If .genome file found in current directory

				>> plink --bfile [input_filename] --read-genome [output_filename].genome --cluster --K [num_clusters] --out [IBS_clusterfilename]
							3                        4                                      5       6                   7

		Operations performed (see equiv. command line argument above): 
			0   loads the binary MAP/PED/etc. PLINK files with the filename [input_filename]
			1   calculates the pairwise identity by state of the samples in the dataset (Note: this can take a while)
			2   outputs a gzipped version of the pairwise IBS results as [output_filename].genome
			3   loads the binary MAP/PED/etc. PLINK files with the filename [input_filename]
			4   loads the (possibley g-zipped) [output_filename].genome file
			5   calls the clustering command to group the samples based on the falgs/parameters and outputs the four output files. 
			6   specifies the number of clusters [num_clusters] into which to group samples  
			7   outputs the clustering results in a file named [IBS_clusterfilename]

		Parameters:
				 --param--                           --type--        --description-- 
			:param job:                              job instance
			:param MapPed_filename:                  string
			:param IBS_pairs_output_filename:        string
			:param num_clusters:                     int
			:param clusterfile_prefix:               string  
		
		Note: 
			- For speed-up, this function can be split into batches, pushed to a cluster and executed in parallel, but this implementation does not do this. 
			- MAF filtering and missingness filtering included in this step, as it doesn't make sense to add a read/write step just for filtering.
			- the --cluster call generates four output files:
				 *.cluster0 - contains some information on the clustering process
				 *.cluster1 - contains information on the final solution, listed by cluster
				 *.cluster2 - contains the same information, but listed one line per individual
				 *.cluster3 - in the same format as cluster2 (one line per individual) but contains all solutions (i.e. every step of the clustering
								from moving from N clusters each of 1 individual (leftmost column after family and individual ID) to 1 cluster (labelled 0) 
								containing all N individuals (the final, rightmost column).

		For more information and examples of the clustering output formats: http://pngu.mgh.harvard.edu/~purcell/plink/strat.shtml
	"""
	bfile = input_args['bfile']
	if input_args['cluster'] == 'yes':
		work_dir = job.fileStore.getLocalTempDir()
		# PLINK_command1 = ['--bfile', 'data/' + bfile, 
		# 				'--Z-genome',
		# 				'--out', input_args['pairs_filename']]
		# PLINK_command2 = ['--bfile', 'data/' + bfile,
		# 			   '--read-genome', input_args['pairs_filename'],
		# 			   '--cluster',
		# 			   '--K', input_args['num_clusters'],
		# 			   '--out', input_args['clusterfiles_prefix']]
		PLINK_command = ['--bfile', 'data/' + bfile,
					   '--cluster',
					   '--K', input_args['num_clusters'],
					   '--out', input_args['clusterfiles_prefix']]
		if (not os.path.exists(input_args['pairs_filename'].join(".cluster2"))) or input_args['override_clusterfiles']== 'yes':
			if os.path.exists(input_args['pairs_filename'].join(".genome")) and input_args['override_clusterfiles']== 'yes': 
				# Make the docker call to cluster the IBS pairs found in the .genome file
				try: 
					print("Clustering samples based on previous pairings...")
					docker_call(work_dir,
						tool_parameters = PLINK_command, 
						tool = 'plink',
						sudo = False)
				except: 
					sys.stderr.write('Running the GWAS pipeline with %s in %s failed.' %(" ".join(PLINK_command), work_dir))
					raise
			else: 
				try: 
					print("Calculating IBS pairs...")
					docker_call(work_dir,
						tool_parameters = PLINK_command1, 
						tool = 'plink',
						sudo = False)
				except: 
					sys.stderr.write('Running the GWAS pipeline with %s in %s failed.' % (" ".join(PLINK_command1), work_dir))
				try: 
					print("Clustering based on pairings...")
					docker_call(work_dir,
						tool_parameters = PLINK_command2, 
						tool = 'plink',
						sudo = False)
				except: 
					sys.stderr.write('Running the GWAS pipeline with %s in %s failed.' % (" ".join(PLINK_command2), work_dir))
					raise
		flist = [input_args['clusterfiles_prefix'] + '.cluster{}'.format(i) for i in range(0,3)]
		fetchContainerFiles(flist)
	job.addFollowOnJobFn(regress, input_args)

def regress(job, input_args):
	"""
		Runs an association analysis (linear or logistic regression, based on the analysis_type parameter given), regressing the specified phenotypes on genotype including the specified covariates
		and taking into consideration the IBS clusters created in the previous step.

		Equivalent PLINK command line arguments: 
			For linear regression: 
				>>plink --bfile [input_filename] --linear --within [IBS_clusterfilename] --covar [covar_file.txt] --covar-name [covar_names] --pheno [pheno_filename.txt] --pheno-name [pheno_variablename] --out [output_filename] --adjust
							0                        1        2                             3                            4                      5                            6                                7                         8
			For logistic regression: 
				>>plink --bfile [input_filename] --logistic --within [IBS_clusterfilename] --covar [covar_file.txt] --covar-name [covar_names] --pheno [pheno_filename.txt] --pheno-name [pheno_variablename] --out [output_filename] --adjust
							0                        1        2                             3                            4                      5                            6                                7                         8
		Operations performed (see equiv. command line argument above): 
			0   loads the binary MAP/PED/etc. PLINK files with the filename [input_filename]
			1   runs a linear or logistic regression
			2   flag to take the IBS clustering into account (can also explicitly test for homogeneity of effect between clusters, but this is not done in this implementation)
			3   flag to take include covariates, taking from [covar_file.txt]
			4   flag to specify a subset of the covariates in [covar_file.txt] using [covar_names] where covar_names is a string composed of a comma separated list of covariates to include
			5   specifies the phenotype file to be used instead of the phenotype in the column 6 of the .ped file. 
			6   specifies the phenotype in the phenotype file to be regressed
			7   Outputs the results into a file with the prefix [output_filename] and suffix ".assoc.linear"
			8   flag to make an additional file with prefix [output_filename] and suffix ".adjust" that has the adjusted (bonferroni) p-value results, sorted by p-value.

		Parameters:  
				  --param--                    --type--      --description--      
			:param MapPed_filename:             string
			:param IBS_clusters_filename:       string      
			:param covariates_filename:         string
			:param covariates:                  string
			:param phenotypes_filename:         string
			:param phenotype:                   string
			:param results_prefix:              string
			:param regression_type:             string
		
		Note: 
			- --linear and --logistic generate files with the extension *.assoc.linear or *.assoc.logistic respectively, the only difference being whether the regression coefficients (--linear) or odds ratios (--logistic) is reported.
			- to include covariates, a specific covariates text file must be created (see the links below for more information)
			- by default, plink uses whatever pheontype is in the .ped file for the regression. However, this implementation assumes that the phenotypes to be used in the regression are specified in a separate file that contains the phenotypes for each sample. 
				- the phenotype file must be a file that contains 3 columns (one row per individual):
					- Family ID
					- Individual ID
					- Phenotype
				- the phenotype file must have a header row
					- the first two variables must be labelled FID and IID
					- all subsequent variable names cannot have any whitespace in them
			- It is possible to specify multiple phenotypes to independently regress using --mpheno, but this implementation performs only a single regression per call to regress().
			- The original PED file must still contain a phenotype in column 6 (even if this is a dummy phenotype, e.g. all missing), unless the --no-pheno flag is given
			- Note: for categorical covariates, a set of binary dummy variables must be created and used as the covariates
			- CAUTION: If an individual is not present in the covariate file, or if the individual has a missing phenotype value (i.e. -9 by default) for the covariate, then that individual is set to missing (i.e. will be excluded from association analysis)
			- if using --assoc instead of linear/logistic, and the phenotype is quantitative (i.e. contains values other than 1, 2, 0 or missing) then PLINK will automatically treat the analysis as a quantitative trait analysis

		For more information: http://pngu.mgh.harvard.edu/~purcell/plink/anal.shtml
		For more information on covariate file and flags: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#covar
	"""
	work_dir = job.fileStore.getLocalTempDir()
	regression_type = input_args['analysis_type']
	if regression_type == "linear": 
		reg = "--linear"
	elif regression_type == "logistic": 
		reg = "--logistic"

	if 'cluster' in input_args: 
		# include clusters
		include_clusters = True
	else:
		include_clusters = False
	if 'covar' in input_args:
		# include covariate
		include_covariates = True 
	else:
		include_covariates = False

	if include_clusters: 
		if include_covariates: 
			PLINK_command = ['--bfile', 'data/' + input_args['bfile'],
			   reg,
			   '--within', 'data/' + input_args['clusterfiles_prefix'] + '.cluster2',
			   '--covar', 'data/' + input_args['covar_file'],
			   '--covar-name', input_args['covar_file'],
			   '--pheno', 'data/' + input_args['pheno'],
			   '--pheno-name', input_args['pheno_name'],

			   '--out', input_args['final_out'],
			   '--adjust']
		else: 
			PLINK_command = ['--bfile', 'data/' + input_args['bfile'],
			   reg,
			   '--within', 'data/' + input_args['clusterfiles_prefix'],
			   '--pheno', 'data/' + input_args['pheno'],
			   '--pheno-name', input_args['pheno_name'],
			   '--out', input_args['final_out'],
			   '--adjust']
	else: 
		if include_covariates: 
			PLINK_command = ['--bfile', 'data/' + input_args['bfile'],
			   reg,
			   '--covar', 'data/' + input_args['covar_file'],
			   '--covar-name', input_args['covar_file'],
			   '--pheno', 'data/' + input_args['pheno'],
			   '--pheno-name', input_args['pheno_name'],
			   '--out', input_args['final_out'],
			   '--adjust']
		else: 
			PLINK_command = ['--bfile', 'data/' + input_args['bfile'],
			   reg,
			   '--pheno', 'data/' + input_args['pheno'],
			   '--pheno-name', input_args['pheno_name'],
			   '--out', input_args['final_out'],
			   '--adjust']

	# Make the docker call with the given command line arguments 
	try: 
		print("Performing regression...")
		docker_call(work_dir,
			tool_parameters = PLINK_command, 
			tool = 'plink',
			sudo = False)
		print("Regression finished.")
	except: 
		sys.stderr.write('Running the GWAS pipeline with %s in %s failed.' %(" ".join(PLINK_command), work_dir))
		raise
	current_directory = os.getcwd()
	files = [f for f in listdir(current_directory) if isfile(join(current_directory, f)) and input_args['final_out'] in f]
	fetchContainerFiles(files)

if __name__ == '__main__':
	"""
	This is a Toil pipeline used to perform GWA studies.	
	"""
	# Define Parser object and add to Toil
	parser = build_parser()
	# parameter_file = parser.parse_args()[0]
	if os.path.isdir('new'): 
		subprocess.check_call(['rm','-r','./new'])
	Job.Runner.addToilOptions(parser)
	args = parser.parse_args()
	parameter_file = args.file
	print("Reading in parameters...")
	inputs = readParameters(parameter_file)
	print("Launching pipeline...")
	# Launch Pipeline
	Job.Runner.startToil(Job.wrapJobFn(start_pipeline, inputs), args)







