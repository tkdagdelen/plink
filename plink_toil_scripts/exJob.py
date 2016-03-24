def haplotype_caller(job, shared_ids, input_args):
    """
    Uses GATK HaplotypeCaller to identify SNPs and Indels and writes a gVCF.
    Calls per-sample genotyper to genotype gVCF.
    :param job: Job instance
    :param shared_ids: dictionary of shared file promises
    :param input_args: dictionary of input arguments
    """
    work_dir = job.fileStore.getLocalTempDir()
    input_files = ['ref.fa', 'ref.fa.fai', 'ref.dict', 'toil.bam', 'toil.bam.bai']
    read_from_filestore_hc(job, work_dir, shared_ids, *input_files)
    output = '%s.raw.BOTH%s.gvcf' % (input_args['uuid'],
                                     input_args['suffix'])
    
    # Call GATK -- HaplotypeCaller
    command = ['-nct', input_args['cpu_count'],
               '-R', 'ref.fa',
               '-T', 'HaplotypeCaller',
               '--genotyping_mode', 'Discovery',
               '--emitRefConfidence', 'GVCF',
               '-I', 'toil.bam',
               '-o', output,
               '-variant_index_type', 'LINEAR',
               '-variant_index_parameter', '128000',
               '--annotation', 'QualByDepth',
               '--annotation', 'DepthPerSampleHC',
               '--annotation', 'FisherStrand',
               '--annotation', 'ReadPosRankSumTest']
    try:
        docker_call(work_dir = work_dir,
                    tool_parameters = command,
                    tool = 'quay.io/ucsc_cgl/gatk',
                    sudo = input_args['sudo'])
    except:
        sys.stderr.write("Running haplotype caller with %s in %s failed." % (
            " ".join(command), work_dir))
        raise

    # Update fileStore and spawn child job


    Plink Functions to Run: 

    Function: Imputation 
    Command Line: 
    Input Files: 
    Output Files: 

    Function: VCF to PED/MAP converion
    Command Line: plink --vcf [VCFfilename]
                                          --double-id
                                          --const-fid {FID}
                                          --id-delim {delimiter}
                                          --id-delim ' '
                                                                  --vcf-idspace-to [character]
                                                                                                --biallelic-only <strict> <list>
                                                                                                --vcf-min-qual [val]
                                                                                                                                  --vcf-filter {exception(s)...}
                                                                                                                                                                  --vcf-require-gt
                                                                                                                                                                                    --vcf-min-gq [val]
                                                                                                                                                                                                      --vcf-min-gp [val] 
                                                                                                                                                                                                                        --out [outfile]

      
      #--vcf loads a (possibly gzipped) VCF file, extracting information which can be represented by the PLINK 1 binary format and ignoring everything else.
      # --double-id causes both family and within-family IDs to be set to the sample ID.
      # --const-fid converts sample IDs to within-family IDs while setting all family IDs to a single value (default '0').
      # --id-delim causes sample IDs to be parsed as [FID][delimiter][IID]; the default delimiter is '_'. If any sample ID does not contain exactly one instance of the delimiter, an error is normally reported; however, if you have simultaneously specified --double-id or --const-fid, PLINK will fall back on that approach to handle zero-delimiter IDs.
      # PLINK sample IDs cannot contain spaces, an error is normally reported when there's a space in a VCF sample ID. To work around this, you can use --vcf-idspace-to to convert all spaces in sample IDs to another character.    
      # when more than one alternate allele is present, the reference allele and the most common alternate are tracked (ties broken in favor of the lower-numbered allele) and the rest are coded as missing calls. To simply skip all variants where at least two alternate alleles are present in the dataset, use --biallelic-only.    
      # --vcf-min-qual causes all variants with QUAL value smaller than the given number, or with no QUAL value at all, to be skipped. (--qual-scores has similar functionality.)
      # To skip variants which failed one or more filters tracked by the FILTER field, use --vcf-filter. This can be combined with one or more (space-delimited) filter names to ignore.
      # By default, when the GT field is absent, the variant is kept and all genotypes are set to missing. To skip the variant instead, use --vcf-require-gt.
      # --vcf-min-gq excludes all genotype calls with GQ below the given (nonnegative, decimal values permitted) threshold. Missing GQ values are not treated as being below the threshold.
      # Similarly, --vcf-min-gp excludes all genotype calls with GP value below the given threshold, assuming GP is 0-1 scaled rather than phred-scaled.
    Input Files: [VCFfilename]
    Output Files: [outfile]

    Function: MAF Filtering 
    Command Line: 
    Input Files: 
    Output Files:

    Function: Clustering/Stratification
    Command Line: 
    Input Files: 
    Output Files:

    Function: Association Analysis (with covariates)
    Command Line: 
    Input Files: 
    Output Files:

    Function: Results Correction (bonferroni)
    Command Line: 
    Input Files: 
    Output Files:

    Function: 
    Command Line: 
    Input Files: 
    Output Files:

    Function: 
    Command Line: 
    Input Files: 
    Output Files:

    Parameters to keep track of: 
    General Params: 
      - whether to do linear or logistic regression
    VCF_toPED: 
        - VCFfilename
        - min quality threshold (what is the default?)
        - min gq score 
        - min gp score
        - prefix for output files
    apply_filters:
        - the prefix for the PED/MAP file 
        - MAF cutoff
        - % samples missing for given SNP cutoff
        - % SMPs missing for ginen sampel cutoff
        - prefix for output files
    IBS_clustering:
        - the prefix for the PED/MAP files 
        - IBS_pairs filename (.genome)
        - number of clusters to make
        - perfix for IBS cluster files
    Assoc: 
        - selection of linear or logistic
        - the prefix for the PED/MAP files 
        - prefix for IBS cluster files 
 -      - covariate file name 
        - covariate names (comma separated list)
 -      - phenotype filename 
        - prefix for output files (.assoc.linear file and .adjust file)

Process: 
      0   Make phenotype file and covariate files (can be the same if only phenotypic covariates will be used, 
            but the implementation can take two different files if instead genotypic covariates are desired)
              -Phonotype file
                - the phenotype file must be a file that contains at least 3 columns (one row per individual):
                      - Family ID
                      - Individual ID
                      - Phenotype(s)
                - the phenotype file must have a header row
                      - the first two variables must be labelled FID and IID
                      - all subsequent variable names cannot have any whitespace in them
              -covariate file
                - same as the phenotype file, except the columns following FID and IID are covariates (not necessarily phenotypic)
                - Note: for categorical covariates (and genotypic covariates) a set of binary dummy variables must be created, the categories must be translated 
                    into those binary variables, and each of the dummy variables should then be included as a covariate
      1   Fill in analysis_parameters into the analysis_parameters text file
      2   Make sure you are in the same directory as the VCF file and parameters file
      3   Run the analysis, specifying the analysis_parameters text file, using the following command: 
            >> toil .... analysis_parameters






