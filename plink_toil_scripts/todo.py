TODO: 
    now: 

    # - Build parser (and type checking)
    - write the toil functionality
    - write ReadMe
    - fix docker for PLINK2
    - put descriptions into the function descriptions
    - send put on GitHub
    - imputation???


    General:
        # - Outline the work flow 0—>1—>2—>3 etc.
        - Outline the individual functions you want to have 
            -What is does it do? 
            -What does it need as inputs?
                call_impute2
                # VCF_to_ PEDMAP
                # Filters
                    # still need to set params
                # IBS_cluster
                    # still need to set params
                # lin_reg
                # log_reg
                # bonferroni
        - arrange all the flags in the order that they are executed in PLINK
        - Go backwards through the workflow and make the output of the previous job the input of the current job.
        - Write the addFollowOnJobFn(…) calls to string together the functions
        - Copy/outline the utility functions you will need for data input, docker calling, etc.
        - Parser for the analysis_parameters file
        - error-handling for invalid input parameters 
        - Pick a subset of the data to use as tests 
        - Write tests for each function (where possible), based on the data you picked
        - For each function
            -implement the function
            -test the function (make sure you are getting the expected output before you move to the next function)