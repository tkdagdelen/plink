Plinkdockerfile
===============

These dockerfiles create PLINK images that can be used on any computer, provided the computer has docker installed. 

Setup and run instructions:
=========================== 


To build PLINK1.07 image, run: 

>>docker build -t plink:latest -f PlinkDockerfile .
(note, the period at the end is necessary)

To build PLINK1.9 aka PLINK2 run: 

>>docker build -t plink2:latest -f Plink2Dockerfile .


To create the docker container, run: 

>>docker run plink

or 

>>docker run plink2 

When no arguments are passed into the containers, the default action for the PLINK1.07 version is to call 

>>plink —noweb 

which should print out something like this: 

—————————————————————————————————————————————————————————————————
|                                                               |

@----------------------------------------------------------@
|        PLINK!       |     v1.07      |   10/Aug/2009     |
|----------------------------------------------------------|
|  (C) 2009 Shaun Purcell, GNU General Public License, v2  |
|----------------------------------------------------------|
|  For documentation, citation & bug-report instructions:  |
|        http://pngu.mgh.harvard.edu/purcell/plink/        |
@----------------------------------------------------------@

Skipping web check... [ --noweb ] 
Writing this text to log file [ plink.log ]
Analysis started: Wed Jun 22 15:50:34 2016

Options in effect:
	--noweb

Before frequency and genotyping pruning, there are 0 SNPs
0 founders and 0 non-founders found
0 SNPs failed missingness test ( GENO > 1 )
0 SNPs failed frequency test ( MAF < 0 )
After frequency and genotyping pruning, there are 0 SNPs

ERROR: Stopping as there are no SNPs left for analysis

|                                                               |
—————————————————————————————————————————————————————————————————

And the default action for the PLINK2 version is:

>>plink

and should produce something like this: 

——————————————————————————————————————————————————————————————————————————————————
|                                                                                |

PLINK v1.90b3.38 64-bit (7 Jun 2016)       https://www.cog-genomics.org/plink2
(C) 2005-2016 Shaun Purcell, Christopher Chang   GNU General Public License v3

  plink [input flag(s)...] {command flag(s)...} {other flag(s)...}
  plink --help {flag name(s)...}

Commands include --make-bed, --recode, --flip-scan, --merge-list,
--write-snplist, --list-duplicate-vars, --freqx, --missing, --test-mishap,
--hardy, --mendel, --ibc, --impute-sex, --indep-pairphase, --r2, --show-tags,
--blocks, --distance, --genome, --homozyg, --make-rel, --make-grm-gz,
--rel-cutoff, --cluster, --pca, --neighbour, --ibs-test, --regress-distance,
--model, --bd, --gxe, --logistic, --dosage, --lasso, --test-missing,
--make-perm-pheno, --tdt, --qfam, --annotate, --clump, --gene-report,
--meta-analysis, --epistasis, --fast-epistasis, and --score.

'plink --help | more' describes all functions (warning: long).

|                                                                                |
——————————————————————————————————————————————————————————————————————————————————

Notes on docker: 
===============

Docker requires a 64-bit installation regardless of your Ubuntu version. Additionally, your kernel must be 3.10 at minimum. The latest 3.10 minor version or a newer maintained version are also acceptable.

Kernels older than 3.10 lack some of the features required to run Docker containers. These older versions are known to have bugs which cause data loss and frequently panic under certain conditions. 





