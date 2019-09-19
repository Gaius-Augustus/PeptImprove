# PeptImprove User Guide

Author and Contact Information
------------------------------

Author: 

Leonie J. Lorenz

Contact: 

Katharina J. Hoff

University of Greifswald,
Institute for Mathematics and Computer Science,
Walther-Rathenau-Str. 47,
17489 Greifswald

University of Greifswald,
Center for Functional Genomics of Microbes,
Felix-Hausdorff-Str. 8,
17489 Greifswald

katharina.hoff@uni-greifswald.de

Contents
========

-   [What is PeptImprove?](#what-is-peptimprove)
-   [Installation](#installation)
    -   [Quick start](#quick-start)
    -   [Dependencies](#dependencies)
    -   [PeptImprove](#peptimprove)
-   [Data preparation](#data-preparation)
-   [Running PeptImprove](#running-peptimprove)
    -   [Creating a new prediction](#creating-a-new-prediction)
    -   [Improving a given reference annotation](#improving-a-given-reference-annotation)
    -   [Options explained](#options-explained)
-   [Example data](#example-data)
-   [Output of PeptImprove](#output-of-peptimprove)
-   [Bug reporting](#bug-reporting)
-   [Citing PeptImprove](#citing-peptimprove)
-   [License](#license)


What is PeptImprove?
================

PeptImprove is a pipeline for making an *ab initio* genome annotations or improving existing annotations by using peptide evidence from MS/MS experiments. If an *ab initio* annotation is desired, the first gene prediction will be done by GeneMarkS-2 <b id="f1">[1]</b> which will then be used to train AUGUSTUS<b id="f2">[2]</b><b id="f3">[3]</b>. AUGUSTUS includes peptide hints into its statistical model to make a gene prediction. In the end both predictions are combined to a final prediction.

Alternatively, it can be chosen to use an existing annotation for training AUGUSTUS which will in the end be combined with the AUGUSTUS prediction.

PeptImprove is a fully automated pipeline for genome annotation and is implemented in Python3. It automatically executes a six-frame translation program called get_orfs_or_cdss.py, the ThermoRawFileParser, IdentiPy<b id="f4">[4]</b>, GeneMarkS-2 and AUGUSTUS.

PeptImprove can be used to either improve an existing reference annotation or to make an annotation from scratch in both cases using peptide evidence given by the user.


Installation
============

Quick Start
-----------

PeptImprove is a Python3 pipeline for Linux with x86-64 architecture. It requires Python3, Biopython<b id="f5">[5]</b>, sort, mono, perl and the following scripts and programs: 

* get_orfs_or_cdss.py
* ThermoRawFileParser
* IdentiPy
* GeneMarkS-2 (if ab initio predictions are desired)
* AUGUSTUS
* blastp


Dependencies
------------

In the following, we give instructions on where dependencies can be obtained, and how they may be installed on Ubuntu Linux.

Python3 is available at <https://www.python.org/downloads/>, or as a package for many Unix systems. Choose version 3.5 or newer (because otherwise, subprocess module is not fully functional).

For example, on Ubuntu, install Python3 with:

```
sudo apt install python3
```

We recommend to use pip for installing further python modules. pip is available at <https://pypi.org/project/pip/>. It is also available as package for many Unix systems.
 
For example, on ubuntu, install pip with:

```
sudo apt install python3-pip
```

The python3 script get_orfs_or_cdss.py uses Biopython. Install Biopython with pip as follows:

```
pip3 install biopython
```

The script get_orfs_or_cdss.py can be downloaded from <https://github.com/peterjc/pico_galaxy/blob/master/tools/get_orfs_or_cdss/get_orfs_or_cdss.py>. The pipeline is proven to work with get_orfs_or_cdss.py in version 0.2.3.

The ThermoRawFileParser can be obtained from <https://github.com/compomics/ThermoRawFileParser>. 

Perl is available at <https://www.perl.org/get.html> or as a package for many Unix systems.

For example, on Ubuntu, install perl with:

```
sudo apt-get install perl 
```
Mono can be dowloaded from <https://www.mono-project.com/download/stable/> where also a detailed description on the dowload procedure is given.

IdentiPy can be obtained from <https://bitbucket.org/levitsky/identipy/src/default/>.

GeneMarkS-2 can be downloaded from <http://topaz.gatech.edu/GeneMark/license_download.cgi>.

The BLAST+ package which includes blastp can be downloaded in the version needed for a certain computer platform from <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/>. A detailed description on how to set up BLAST+ can be found here: <https://www.ncbi.nlm.nih.gov/books/NBK52640/>.

PeptImprove runs AUGUSTUS, so its binaries are necessary, and uses the following tools provided by AUGUSTUS at <https://github.com/Gaius-Augustus/Augustus/tree/master/scripts>:


* aa2nonred.pl
* computeFlankingRegion.pl
* getAnnoFastaFromJoingenes.pl
* gff2gbSmallDNA.pl 
* new_species.pl
* optimize_augustus.pl
* randomSplit.pl


PeptImprove uses Unix sort. sort should be installed by default on all Unix systems.


PeptImprove
-------

PeptImprove is a python3 script named PeptImprove.py. It does not require a particular installation procedure after download.

It can be executed either with

```
python3 PeptImprove.py
```

If you add PeptImprove.py to your $PATH (i.e. by adding the location of PeptImprove.py at the bottom of your ~/.bashrc file similar to ```PATH=/path/to/PeptImprove:$PATH```, followed by loading the ~/.bashrc file in case you did not re-open a new bash session with ```source ~/.bashrc```)and make it executable (i.e. with ```chmod u+x PeptImprove.py```), it can be executed with

```
PeptImprove.py
```

from  any location on your computer.

Please also add Augustus/bin/, IdentiPy, and GeneMarkS-2 (if you desire to use it) to your $PATH in the same way. The directory containing the required python script get_orfs_or_cdss.py can either be stored in the $PATH or the script can simply be stored in PeptImprove's directory. The directory of the ThermoRawFileParser can either be given when calling PeptImprove using the option --path_to_rawfileparser string or it can be stored in PeptImprove's directory, too.

The following PeptImprove scripts are essential for running PeptImprove and can be downloaded from <https://github.com/Gaius-Augustus/PeptImprove>:

* AugOutputToGff.py
* collapseRedundantHints.py
* combineGMandAUGUSTUS2.py
* compareCombinedResultsToRefAnno.py
* GMOutputToGff.py
* HelpForIdentiPy.py
* HelpForThermoRawFileParserMgf.py
* IdentiPyToGffAndFilterAndstartHints.py
* PeptImprove.py

All PeptImprove scripts need to be stored in the same directory as PeptImprove.py.

Data Preparation
==================

PeptImprove accepts files in the following formats:

* genome file in FASTA format (simple FASTA headers without whitespaces or special characters)
* one or multiple directories containing peptide evidence from MS/MS experiments in RAW file format
* reference annotation in GTF format (a detailed description and examples can be found at <https://www.ensembl.org/info/website/upload/gff.html>)
* a (modified) version of the AUGUSTUS extrinsic.cfg file (the original as used by the pipeline can be found in the directory ```config```)

Running PeptImprove
===============

PeptImprove can be used to make a gene prediction from scratch using GeneMarkS-2 and AUGUSTUS or to improve a given reference annotation, then only AUGUSTUS will be run for making a gene prediction.

Making a new gene prediction
------------------

The essential arguments for making a gene prediction with PeptImprove are:

* ```-g GENOME_FILE, --genome_file GENOME_FILE```
  Genome input file in FASTA format
* ```-e EXTRINSIC_EVIDENCE [EXTRINSIC_EVIDENCE ...], --extrinsic_evidence EXTRINSIC_EVIDENCE [EXTRINSIC_EVIDENCE ...]```
 The location of the directory containing all MS data. You may name multiple directories. 
* ```-n NAME, --name NAME```
  Please give the name of your species (avoid white spaces as it will also function as the name of a directory).



Improve a given reference annotation
------------------

If you'd like to improve an existing an existing reference annotation instead of making a new gene prediction with GeneMarkS-2, then you will have to give the following arguments:

* ```-g GENOME_FILE, --genome_file GENOME_FILE```
  Genome input file in FASTA format
* ```-e EXTRINSIC_EVIDENCE, --extrinsic_evidence EXTRINSIC_EVIDENCE```
 The location of one or multiple directories containing all MS data. 
* ```-n NAME, --name NAME```
  Please give the name of your species (avoid white spaces as it will also function as the name of a directory).
* ```-r REFERENCE_ANNOTATION, --reference_annotation REFERENCE_ANNOTATION```
  Reference annotation is gtf format. Can be used for comparing the results of for an different usage of the pipeline.
* ```--use_ref```
  Option to use reference annotation for training AUGUSTUS. GeneMarkS-2 won't be run in this version.

Please have a look at the section  [Options Explained](#options_explained) for information about options to further adapt the pipeline to your wishes. You may, for example, state how many cores are available to the pipeline, into which directory the output files should be written and whether temporary files should be deleted.

Usage example 1:

```
PeptImprove.py -g genome.fa -e DirectoryWithExtrinsicDataFromMSExperiments -n name_of_your_species 
```
The pipeline will create hints from the given raw data from MS/MS experiments, run GeneMarkS-2, train and optimize AUGUSTUS based on the GeneMarkS-2 prediction, make a prediction with AUGUSTUS based on the given hints and join the GeneMarkS-2 and the AUGUSTUS prediction.

Usage example 2:

```
PeptImprove.py -g genome.fa -e Directory1WithExtrinsicDataFromMSExperiments Directory2WithExtrinsicDataFromMSExperiments Directory3WithExtrinsicDataFromMSExperiments -n name_of_your_species -r referenceAnnotation.gtf --use_ref
```
The pipeline will create hints from the given raw files which in this case are stored in three different directories, will train and optimize AUGUSTUS based on the given reference annotation, make a prediction with AUGUSTUS based on the given hints and then join the AUGUSTUS prediction with the reference annotation.

Options explained
-----------------

In the following, we explain all options of PeptImprove.py:

* ```-h, --help```
  show this help message and exit
* ```-g GENOME_FILE, --genome_file GENOME_FILE```
  Genome input file in FASTA format
* ```-e EXTRINSIC_EVIDENCE [EXTRINSIC_EVIDENCE ...], --extrinsic_evidence EXTRINSIC_EVIDENCE [EXTRINSIC_EVIDENCE ...]```
 The location of the directory containing all MS data. You may name multiple directories. 
* ```-n NAME, --name NAME```
  Please give the name of your species (avoid white spaces as it will also function as the name of a directory).
  * ```-r REFERENCE_ANNOTATION, --reference_annotation REFERENCE_ANNOTATION```
  Reference annotation is gff format. Can be used for comparing the results of for an different usage of the pipeline.
* ```--use_ref```
  Option to use reference annotation for training AUGUSTUS. GeneMarkS-2 won't be run in this version.
* ```--cores CORES```
  If you have information about the number of available cores you can enter it here. Default will be 4.
* ```--re_run```
  Option to leave out GeneMark prediction and AUGUSTUS training. Only makes sense when the pipeline was used before on the same species but with different extrinsic evidence.
* ```--no_optimization```
  Option to omit optimize_augustus. Only advisable if there is not enough time or computational capacity.
* ```--output_directory OUTPUT_DIRECTORY```
  Option to name the output directory and to define                    where it should be. If this argument is not given the                        default directory 'output' will be created where the                         pipeline was called.
* ```--augustus_config_path AUGUSTUS_CONFIG_PATH``` Option give the path to the AUGUSTUS config directory.
* ```--path_to_rawfileparser PATH_TO_RAWFILEPARSER``` Option give the path to the directory of the ThermoRawFileParser.
* ```--delete_intermediate```
  Option to delete all results except the final one.
* ```--archaea``` Use option if species is of type Archaea.
* ```--translation_table TRANSLATION TABLE``` If the translation table of the species' genome is not 11, please enter the number her. Default will be 11.
* ```--start START```
  You may give an integer (incl. 0) that indicates how many steps should be left out at the beginning of the pipeline.
* ```--stop STOP```
  You may give an integer <=9 that indicates how many steps should be left out at the end of the pipeline.
* ```-c, --compare_to_ref```
  Option to compare results to a reference annotation. Only possible when the reference annotation is given via -r.
* ```-f FILTER_IP, --filter_ip FILTER_IP```
  Value at which the hints by IdentiPy should be filtered. 0.01 stands for 1 percent. Default is 0.001.
* ```--extrinsic_cfg EXTRINSIC_CFG```
  Option to use a specific Cfg file that is derived from the original one. Elsewise the pipeline's default Cfg file is used.


Example data
============

Example data is located in the directory ```example/```. It consists of the following files:

* ```GCF_000953275.1_CD630DERM_genomic_example.fna```: the first 800,000 base pairs of sequence NZ_LN614756.1 of *Clostridium difficile 630Δerm*, assembly version GCF_000953275.1 from the NCBI.

* ```ReferenceAnnotation_example.gtf```: all reference annotation's genes corresponding to the first 800,000 base pairs of NZ_LN614756.1 of *Clostridium difficile 630Δerm*, of reference annotation version GCF_000953275.1 from the NCBI

Due to the size of raw format files, the example file is not contained in the named directory, but can be downloaded from <ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2016/08/PXD003659/150505_OV3_P5_SK_SSW_AO_BHI_I_01.raw>. It is a raw file from MS/MS experiments done with *Clostridium difficile 630Δerm* by Otto et al. (2016) <b id="f6">[6]</b>.

There are to options to run the example data, either creating a new gene prediction or improving the above named reference annotation. Both runs take about 4 minutes when run on four cores (default setting of PeptImprove) and optimize_augustus.pl is skipped.

Example 1:
python3 PeptImprove.py -g GCF_000953275.1_CD630DERM_genomic_example.fna -n clost_diff_example1 --no_optimization --path_to_rawfileparser PATH_TO_RAWFILEPARSER

Example 2:
python3 PeptImprove.py --genome GCF_000953275.1_CD630DERM_genomic_example.fna --reference_annotation ReferenceAnnotation_example.gtf --use_ref --name clost_diff_example2 --no_optimization --path_to_rawfileparser PATH_TO_RAWFILEPARSER --output_directory PATH_TO_OUTPUT

Output of MakeHub
=================

PeptImprove.py creates a directory that is either named and located as stated via the option --output_directory PATH_TO_OUTPUT or, if no output directory is given, it will create the default output directory "output" as a subdirectory of PeptImprove's directory.

```output``` contains the following files:

* ```PeptImprove_final_prediction.gff``` - this file contains PeptImproves final gene prediction. It is either a combination of the AUGUSTUS prediction with the GeneMarkS-2 prediction of with the given reference annotation, if --use_ref was chosen.

All other files and directories in the ```output``` directory are temporary results which are deleted by the pipeline if option --delete_intermediate is chosen. The most important and possibly interesting intermediate results are the following:

* ```GffFromGMS2.gff``` - this file contains the GeneMarkS-2 gene prediction in GFF format is option --use_ref is not chosen, otherwise this file will not exist.

* ```augustus.gff``` - this file contains the AUGUSTUS gene prediction in GFFF format.

*  ```HintsByIP.gff``` - this file contains all hints created from peptive evidence by IdentiPy, namely CDSpart and start codon hints with different multiplicities.
  
* if option ```-c, --compare_to_ref``` and a reference annotation are given, the file ```ComparisonToReferenceAnnotation.txt``` exists and contains information on the sensitivity and specificity levels of PeptImproves final prediction as when compared to the reference annotation.


Bug reporting
=============

Before reporting bugs, please check that you are using the most recent versions of PeptImprove. Also, check the open and closed issues on github at <https://github.com/Gaius-Augustus/PeptImprove/issues> for possible solutions to your problem.

Reporting bugs on github
------------------------

If you found a bug, please open an issue at <https://github.com/Gaius-Augustus/PeptImprove/issues>  (or contact katharina.hoff@uni-greifswald.de).

Information worth mentioning in your bug report:

PeptImprove.py prints information about separate steps on STDOUT. Please let us know at which step and with what error message PeptImprove.py caused problems.


Citing PeptImprove
==============

Lorenz LJ, “A pipeline for prokaryotic genome annotation with peptide data from MS/MS experiments”. Bachelor thesis. 

License
=======

All source code of PeptImprove is under GNU public license 3.0 (see
<https://www.gnu.org/licenses/gpl-3.0.de.html>).

Please note that this is not true for all external programs that are used in the pipeline. Always refer to the license notice shown at each program's download location that are named above.

References
==========

<b id="f1">[1]</b> Lomsadze, A., Gemayel, K., Tang, S., & Borodovsky, M. (2018). Modeling leaderless transcription and atypical genes results in more accurate gene prediction in prokaryotes. *Genome Reseach*. doi: 10.1101/gr.230615.117

<b id="f2">[2]</b> Stanke, M., Steinkamp, R., Waack, S., & Morgenstern, B. (2004). AUGUSTUS: a web server for gene finding in eukaryotes. *Nucleic Acids Research*, Volume 32, Issue suppl_2, Pages W309-W312. https://doi.org/10.1093/nar/gkh379

<b id="f3">[3]</b> Hoff, K. J. & Stanke, M. (2019). Predicting genes in single genomes with AUGUSTUS. *Current Protocols in Bioinformatics*, 65, e57. doi: 10.1002/cpbi.57

<b id="f4">[4]</b> Levitsky, L. I., Ivanov, M. V., Lobas, A. A., Bubis, J. A.,Tarasova, I. A., Solovyeva, E. M., Pridatchenko, M. L., & Gorshkov, M. V. (2018). IdentiPy: An Extensible Search Engine for Protein Identication in Shotgun Proteomics. *Journal of Proteome Research*. doi: https://doi.org/10.1021/acs.jproteome.7b00640

<b id="f5">[5]</b> Cock, P. J. A., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., Friedberg, I., Hamelryck, T., Kauff, F., Wilczynski, B., & de Hoon, M. J. L. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics.
*Bioinformatics*, 25(11) 1422-3. https://doi.org/10.1093/bioinformatics/btp163 pmid:19304878

<b id="f6">[6]</b> Otto, A., Maaß, S., Lassek, C., Becher, D., Hecker, M., Riedel, K., & Sievers, S. (2016). The protein inventory of Clostridium diffcile grown in complex and minimal medium.
*Proteomics Clinical Applications*. https://doi.org/10.1002/prca.201600069

