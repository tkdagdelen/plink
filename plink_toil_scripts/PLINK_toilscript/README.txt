README for PLINK_toilscript

Note: in order for the PLINK_toilscript to work, you need TOIL, Docker, and Python installed. 

The the script was made to make running a GWAS on dbGaP data very easy. Analyses can be run using a simple command line call. Instead of using flags in the call, all options are read in from a parameters file. Hopefully this removes the need to copy and paste long command line calls. 

The file that specifices all the parameters should be in the format specified by the template. In it, the mandatory parameters are at the top, and the optional parameters are at the bottom. For the optional parameters, once the first parameter is "yes" then the rest of the parameters in that block are are required, unless a default value is given in the description. 

Note that there is a parameter in each optional block for whether or not to override past files with the same name. This allows you to override the files if you have changed parameters upstream of that step, or use the old files instead of recomputing them if you are only changing parameters in downstream steps. This is particularly important for the clustering operation. If you don't need to recompute the clusters, you shouldn't, since it takes a long time to run. 

To run the script on dbGaP data, follow the following procedure: 

	1) Ensure that the relevant dbGaP data is decrypted and in the current directory. Also ensure that the PLINK_toilscript.py and the relevant parameters file are in the current directory. 
	2) Open the parameters file, input the necessary information for all the mandatory parameters and whatever other relevant parameters, and save. 
	3) run the script using the online command: 

		>> python PLINK_toilscript.py -f myParameterfile.txt ./new

* If it breaks in the middle of a run, you may need to remove the old docker container before you can try running it again. To remove the containter, use Docker rm <container_name>

