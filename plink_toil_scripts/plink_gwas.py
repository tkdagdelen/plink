#!/usr/bin/env python2.7
# Taner Dagdelen
# Feb. 22, 2016
# Aknowledgements: 
#       1) Utility functions copied from John Vivian's batch_align.py script

"""
Genome Wide Association Study pipeline with MAF and missingness thresholding, IBS_clustering for stratification correction, association analysis using linear/logistic regression, bonferroni correction, and sorting by P-value. 
This pipeline is intended for use with PLINK2 (1.9) 
[link to PLINK docker script]


 		0 --> 1 --> 2 --> 3 --> 4

 0	Imputation ? (might move this out into its own script later)
 1	VCF-->PED/MAP conversion
 2 	MAF and Missingness Filtering
 3  Clustering/Stratification
 4	Association Analysis (including covariates) & Results Correction (bonferroni)

Flags: 
- no flags. Instead, put input parameters into the GWAS_parameters_template, save under a new filename, and put the path to that file
    in the command line call 

Dependencies: 
Docker  -   apt-get install docker.io
Toil    -   pip install toil
PLINK   -   

 """
##################### Parser #######################

 def build_parser():
    parser = argparse.ArgumentParser(description=main.__doc__, add_help=True)
    parser.add_argument('-v', '--vcf_filename', required=True, help='VCF_filename (string, ex. myVCF.vcf)')
    parser.add_argument('-q', '--quality', required=True, help='Minimum genotype quality threshold (float, default is )')
    parser.add_argument('-gp', '--min_gp_score', required=True, help='Minimum gp_score threshold (float, default is )')
    parser.add_argument('-gq', '--min_gq_score', required=True, help='Minimum gq_score threshold (float, default is )')
    parser.add_argument('-f', '--filename', required=True, help='Desired prefix for original PED/MAP files (string, ex. SUPRESS_trial_data)')
    parser.add_argument('-maf', '--maf_cutoff', required=True, help='Minor Allele Frequency filter MAFcutoff (float, default is )')
    parser.add_argument('-m', '--mind', required=True, help='Cutoff for percentage of samples missing for a given SNP (float, default is )')
    parser.add_argument('-g', '--geno', required=True, help='Cutoff for percentage of SNPs missing for a given sample (float, default is )')
    parser.add_argument('-ff', '--filtered_filename', required=True, help='Desired prefix for post-filter PED/MAP files (string, ex. SUPRESS_trial_data_filtered)')
    parser.add_argument('-ibs', '--ibs_pairsfile_prefix', required=False, help='Desired prefix for the file containing non-clustered IBS pairs (*.genome) (string, ex. SUPRESS_trial_data_filtered)')
    parser.add_argument('-n', '--num_clusters', required=True, help='Number of clusters to sort IBS pairs into (float, ex. 3)')
    parser.add_argument('-c', '--clusterfiles_prefix', required=True, help='Desired prefix for file containing final IBS clusters (string, ex. SUPRESS_clusters)')
    parser.add_argument('-u', '--sudo', dest='sudo', action='store_true', help='Docker usually needs sudo to execute '
                                                                               'locally, but not''when running Mesos '
                                                                               'or when a member of a Docker group.')
    parser.add_argument('-3', '--s3_dir', default=None, help='S3 Directory, starting with bucket name. e.g.: '
                                                             'cgl-driver-projects/ckcc/rna-seq-samples/')
    parser.add_argument('-a', '--analysis_type', required=True, help='Type of analysis to run (linear or logistic) (ex. linear)')
    parser.add_argument('-se', '--file_size', default='100G', help='Approximate input file size. Should be given as %d[TGMK], e.g., for a 100 gigabyte file, use --file_size 100G')
    parser.add_argument('--covf', '--covariates_filename', required=True, help='Name of covariate file (string, ex. covar.txt)')
    parser.add_argument('--covars', '--covariate_names', required=True, help='Comma separated list of covariates to include from covariate file (string, ex. BMI,AGE,Treatment_Arm)')
    parser.add_argument('-phenof', '--phenotypes_filename', required=True, help='Name of phenotype file contianing phenotype(s) to regress on (string, ex. pheno.txt)')
    parser.add_argument('-pheno', '--phenotype_name', required=True, help='Name of phenotype variable from phenotype file to regress (string, ex. 10yr_survival)')
    parser.add_argument('-o', '--output_filename', required=True, help='Desired prefix for final association .assoc.linear/logistic and .adjust output files (string, ex. SUPRESS_results)')
    parser.set_defaults(sudo=False)
    return parser

################# Utility functions #####################

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

def download_from_url(job, url, fname):
    """
    Downloads a file from a URL and places it in the jobStore
    Input1: Toil job instance
    Input2: Input arguments
    Input3: jobstore id dictionary
    Input4: Name of key used to access url in input_args
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_path = os.path.join(work_dir, fname)
    if not os.path.exists(file_path):
        try:
            subprocess.check_call(['curl', '-fs', '--retry', '5', '--create-dir', url, '-o', file_path])
        except OSError:
            raise RuntimeError('Failed to find "curl". Install via "apt-get install curl"')
    assert os.path.exists(file_path)
    return job.fileStore.writeGlobalFile(file_path)

def copy_to_output_dir(work_dir, output_dir, uuid=None, files=()):
    """
    A list of files to move from work_dir to output_dir.
    work_dir: str       Current working directory
    output_dir: str     Desired output directory
    files: list         List of files (by filename) to be copied to the desired output directory
    """
    for fname in files:
        shutil.copy(os.path.join(work_dir, fname), os.path.join(output_dir, fname))

######################## Jobs ###########################

# Start the pipeline 
def start_pipeline(job, input_args):
    """
    Downloads the data that will be used in the analysis. 
    """
    files_to_downlaod = "blah.vcf"
    url = input_args[files_to_downlaod]
    data = job.addChildJobFn(download_from_url, url, files_to_downlaod)
    job.addFollowOnJobFn(call_impute2, data)

def call_impute2(job, data): 
    """
    NEED TO FIGURE OUT HOW TO MAKE A CALL TO INPUTE2

    IMPUTE2 output to VCF conversion: 
        https://www.biostars.org/p/153728/
        https://samtools.github.io/bcftools/bcftools.html#convert

    """
    job.addFollowOnJobFn(VCF_to_PEDMAP, VCFfile_path, flags)


def VCF_to_PEDMAP(job, VCFfilename, min_qual, min_gq, min_gp, output_prefix):
    """
        Uses PLINK's built-in function to read VCF and write to PED/MAP files. 
         
        Equivalent PLINK command line arguments: 

            >> plink --vcf [VCF_filename] --id-idspace-to _ --vcf-min-qual [min_qual_threshold] --vcf-require-gt --missing-genotype - --vcf-min-gq [min_gq_score] --vcf-min-gp [min_gp_score] --vcf-half-call missing --const-fid {"0"} --out [output_filename] 
                        0                   1                2                                   3                  4                     5                       6                            7                         8                9         
        
        Operations performed (see equiv. command line argument above): 
            0   loads a (possibly gzipped) VCF file, extracting information which can be represented by the PLINK 1 binary format and ignoring everything else.
            1   replaces all the spaces in the sampleIDs with "_"
            2   filters out variants with QUAL value smaller than [val], or with no QUAL value at all.
            3   filters out any SNPs without "GT" in the INFO column, indicating that the genotype of the sample at that site is not included in that row. (i.e. row without genotype information)
            4   set the character to be used to denote missing genotype calls to be "-"
            5   excludes all genotype calls with GQ below the given (nonnegative, decimal values permitted) threshold. Note: Missing GQ values are not treated as being below the threshold.
            6   excludes all genotype calls with GP value below the given threshold, assuming GP is 0-1 scaled rather than phred-scaled.
            7   treat half-calls ('0/.' and similar GT values) as missing.
            8   sets the family id's to "0" and the within-family id's to the sampleIDs of the VCF file.
            9   sets the name of the output file
        
        Parameters:
                 --param--                  --type--        --description--  
            :param job:                     Job instance
            :param VCFfilename:             string
            :param min_qual:                float
            :param min_gq:                  float
            :param min_gp:                  float
            :param output_prefix:           string

        Note: 
            - phase and dosage information is discarded in the conversion (a consequence of PLINK's built-in conversion function) 
            - VCF reference alleles are set to A2, even if minor (PLINK normally forces major alleles to A2 in loading step)
            - VCFs use sample IDs, but PLINK tracks family and within-family ID's. This implementation sets family ID to '0' and just 
                uses the VCF's sample IDs as the within-family IDs. 
            - Also note that sample IDs cannot contain spaces. To prevent this, --id-idspace-to ["_"] is used in this implementation to convert spaces to underscores in the sampleIDs. 
            - By default, PLINK only tracks the reference allele and the most common alternate allele found in the sample (or samples 
                if the VCF is the product of an aggregation of samples). PLINK codes all the rest of the alternate alleles as missing calls (missing calls denoted by "-" in this implementation)
        
        For more information on the VCF format: http://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it
        For more inromation about PLINK's VCF conversion peculiarities/options: https://www.cog-genomics.org/plink2/input#vcf
       
    """
    # Set load necessary input parameters    
    VCF_filename = input_args['vcf_filename']
    min_qual_threshold = input_args['quality']
    min_gq_score = input_args['min_gp_score']
    min_gp_score = input_args['min_gq_score']
    output_filename = input_args['filename']
    work_dir = job.fileStore.getLocalTempDir()
    
    # Put together PLINK arguments
    PLINK_command = ['--vcf', VCF_filename,
               '--id-idspace-to', '_',
               '--vcf-min-qual', min_qual_threshold,
               '--vcf-require-gt', 
               '--missing-genotype', '-',
               '--vcf-min-gq', min_gq_score,
               '--vcf-min-gp', min_gp_score,
               '--vcf-half-call', 'missing',
               '--const-fid', '{"0"}',
               '--out', output_filename]
    
    # Make the docker call with the given command line arguments 
    try: 
        docker_call(work_dir,
            tool_parameters = PLINK_command, 
            tool = 'plink',
            sudo = input_args['sudo'])
    except: 
        sys.stderr.write('Running the GWAS pipeline with %s in %s failed.' %(" ".join(PLINK_command), work_dir))
        raise

    # Link to the next job in the pipeline
    job.addFollowOnJobFn(apply_filters, input_args)


def apply_filters(job, MapPed_filename, MAF_cutoff, percent_SNPs_missing_cutoff, genoRate_cutoff, output_prefix): 
	"""
	This module applies MINOR Allele Frequency and missingness filters to the data and outputs a new file. 

    Equivalent PLINK command line arguments:

        >> plink --bfile [input_filename] --maf [mafthresh] --mind [SNPmiss] --geno [genoRate] --out [output_filename]


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
                    GWAS is MAF<1%. 
        --mind   Filters out sampels in which the percentage of missing SNPs is greater than missingnessThresh. Standard 
                    according to PLINK website is 10%, which is also the default value. 
        --geno   Filters out SNPs for which the genotyping rate is less than rateThresh. Standard genotyping rate cutoff according to 
                    PLINK website is 90%. Default is to include all SNPs. 

    Note: 
        - This filter step can be performed in the VCF to PED/MAP conversion step but in this implementation it is kept separate for simplicity. 

    For more information: http://pngu.mgh.harvard.edu/~purcell/plink/thresh.shtml
	"""
    # Set load necessary input parameters
    MapPed_filename = input_args['filename']
    maf_cutoff = input_args['maf_cutoff']
    percent_SNPs_missing_cutoff = input_args['mind']
    genoRate_cutoff = input_args['geno']
    output_prefix = input_args['filtered_filename']
    work_dir = job.fileStore.getLocalTempDir()
    
    # Put together PLINK arguments
    PLINK_command = ['--bfile', MapPed_filename,
               '--maf', maf_cutoff,
               '--mind', percent_SNPs_missing_cutoff,
               '--geno', genoRate_cutoff,
               '--out', output_prefix]
    
    # Make the docker call with the given command line arguments 
    try: 
        docker_call(work_dir,
            tool_parameters = PLINK_command, 
            tool = 'plink',
            sudo = input_args['sudo'])
    except: 
        sys.stderr.write('Running the GWAS pipeline with %s in %s failed.' %(" ".join(PLINK_command), work_dir))
        raise

    # Link to the next job in the pipeline
    job.addFollowOnJobFn(IBS_cluster, input_args)

# This function can be done over a cluster 
def IBS_cluster(job, MapPed_filename, IBS_pairsfile_prefix, num_clusters, clusterfile_prefix):
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

    # Set load necessary input parameters
    MapPed_filename = input_args['filtered_filename']
    IBS_pairsfile_prefix = input_args['ibs_pairsfile_prefix']
    num_clusters = input_args['num_clusters']
    clusterfiles_prefix = input_args['clusterfiles_prefix']
    work_dir = job.fileStore.getLocalTempDir()
    
    # Put together PLINK arguments
    PLINK_command1 = ['--bfile', MapPed_filename, 
                    '--Z-genome',
                    '--out', IBS_pairsfile_prefix
    PLINK_command2 = ['--bfile', MapPed_filename,
                   '--read-genome', clusterfiles_prefix,
                   '--cluster',
                   '--K', num_clusters,
                   '--out', clusterfiles_prefix]
    # If the pairwise IBS calculation was already done, skip that step and just read in the IBS pairs from the 
    # previously created file to perform clustering.
    if os.path.exists("{0}".format(IBS_pairs_output_filename).join(".genome")): 
        # Make the docker call to cluster the IBS pairs found in the .genome file
        try: 
            docker_call(work_dir,
                tool_parameters = PLINK_command, 
                tool = 'plink',
                sudo = input_args['sudo'])
        except: 
            sys.stderr.write('Running the GWAS pipeline with %s in %s failed.' %(" ".join(PLINK_command), work_dir))
            raise
    else: 
        try: 
            docker_call(work_dir,
                tool_parameters = PLINK_command1, 
                tool = 'plink',
                sudo = input_args['sudo'])
            docker_call(work_dir,
                tool_parameters = PLINK_command2, 
                tool = 'plink',
                sudo = input_args['sudo'])
        except: 
            sys.stderr.write('Running the GWAS pipeline with %s in %s failed.' % (" ".join(PLINK_command), work_dir))
            raise
        # Make docker call to create the .genome file of IBS pairs and make a docker call to cluster them.

    # Link to the next job in the pipeline
    job.addFollowOnJobFn(regress, input_args)




def regress(job, MapPed_filename, IBS_clusters_filename, covariates_filename, covariates, phenotypes_filename, phenotype, results_prefix, regression_type):
    """
        Runs an association analysis (linear or logistic regression, based on the reg_type parameter given), regressing the specified phenotypes on genotype including the specified covariates
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

    # Set load necessary input parameters
    MapPed_filename = input_args['filtered_filename']
    clusterfiles_prefix = input_args['clusterfiles_prefix']
    covariates_filename = input_args['covariates_filename']
    covariates = input_args['covariate_names']
    phenotypes_filename = input_args['phenotypes_filename']
    phenotype = input_args['phenotype_name']
    results_prefix = input_args['output_filename']
    regression_type = input_args['analysis_type']
    work_dir = job.fileStore.getLocalTempDir()
    
    if regression_type == "linear": 
        reg = "--linear"
    elif regression_type == "logistic": 
        reg = "--logistic"

    # Put together PLINK arguments
>>plink --bfile [input_filename] --linear --within [IBS_clusterfilename] --covar [covar_file.txt] --covar-name [covar_names] --pheno [pheno_filename.txt] --pheno-name [pheno_variablename] --out [output_filename] --adjust
    PLINK_command = ['--bfile', MapPed_filename,
               '--maf', maf_cutoff,
               '--mind', percent_SNPs_missing_cutoff,
               '--geno', genoRate_cutoff,
               '--out', output_prefix]


    
    # Make the docker call with the given command line arguments 
    try: 
        docker_call(work_dir,
            tool_parameters = PLINK_command, 
            tool = 'plink',
            sudo = input_args['sudo'])
    except: 
        sys.stderr.write('Running the GWAS pipeline with %s in %s failed.' %(" ".join(PLINK_command), work_dir))
        raise

    # Link to the next job in the pipeline
    job.addFollowOnJobFn(IBS_cluster, input_args)

def main():
    """
    This is a Toil pipeline used to perform GWA studies, starting with vcfs.
    It uses IMPUTE2 to impute and uses PLINK to convert from vcf to PLINK's filesystem and to run all the subsequent necessary operations.
    
    """
    # Define Parser object and add to Toil
    parser = build_parser()
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()

    # Store input parameters in a dictionary
    inputs = {'vcf_filename': args.vcf_filename,
              'quality': args.quality,
              'min_gp_score': args.min_gp_score,
              'min_gq_score': args.min_gq_score,
              'filename': args.filename,
              'maf_cutoff': args.maf_cutoff,
              'mind': args.mind,
              'geno': args.geno,
              'filtered_filename': args.filtered_filename,
              'ibs_pairsfile_prefix': args.ibs_pairsfile_prefix,
              'num_clusters': args.num_clusters,
              'sudo': args.sudo,
              's3_dir': args.s3_dir,
              'clusterfiles_prefix': args.clusterfiles_prefix,
              'analysis_type': args.analysis_type,
              'file_size': args.file_size,
              'covariates_filename': args.covariates_filename,
              'covariate_names': args.covariate_names,
              'phenotypes_filename': args.phenotypes_filename,
              'phenotype_name': args.phenotype_name,
              'output_filename': args.output_filename}

    # Launch Pipeline
    Job.Runner.startToil(Job.wrapJobFn(start_pipeline, inputs), args)






