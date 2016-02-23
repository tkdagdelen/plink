#!/usr/bin/env python2.7
# Taner Dagdelen
# Feb. 22, 2016

"""
Genome Wide Association Study pipeline with MAF and missingness thresholding, stratification correction, association analysis, bonferroni correction, and sorting by P-value. 
This pipeling was created using PLINK [version] 
[link to PLINK docker script]

Need Imputation funciontality, VCF to PED/MAP conversion

 		0 --> 1 --> 2 --> 3 --> 4

 0	Imputation [must write docker file for imputation...]
 1	VCF-->PED/MAP conversion [???]
 2 	MAF and missingness thresholding
 3 	
 4	
 5	
 """

def build_parser():
	parser = argparse.ArgumentParser(description=main.__doc__, add_help=True)
	parser.add_argument()
	return parser

# Utility functions
def docker_call(work_dir, tool_parameters, tool, java_opts=None, outfile=None, sudo=False):
    """
    Makes subprocess call of a command to a docker container.
    tool_parameters: list   An array of the parameters to be passed to the tool
    tool: str               Name of the Docker image to be used (e.g. quay.io/ucsc_cgl/samtools)
    java_opts: str          Optional commands to pass to a java jar execution. (e.g. '-Xmx15G')
    outfile: file           Filehandle that stderr will be passed to
    sudo: bool              If the user wants the docker command executed as sudo
    """
    base_docker_call = 'docker run --rm --log-driver=none -v {}:/data'.format(work_dir).split()
    if sudo:
        base_docker_call = ['sudo'] + base_docker_call
    if java_opts:
        base_docker_call = base_docker_call + ['-e', 'JAVA_OPTS={}'.format(java_opts)]
    try:
        if outfile:
            subprocess.check_call(base_docker_call + [tool] + tool_parameters, stdout=outfile)
        else:
            subprocess.check_call(base_docker_call + [tool] + tool_parameters)
    except subprocess.CalledProcessError:
        raise RuntimeError('docker command returned a non-zero exit status. Check error logs.')
    except OSError:
        raise RuntimeError('docker not found on system. Install on all nodes.')


# Job Functions
#def impute(job)
#def convertVCFtoPEDMAP(job)
def pruneData(job, MAFthresh, missingThresh): 
	"""
	This module filters out SNPs with MAF below the given threshold and filters out samples 
	with missingness above the given threshold. 

	MAFThresh: float			The cutoff for Minor Allele Frequency
	missingThresh: float		THe cutoff for SNP missingness in samples
	"""






