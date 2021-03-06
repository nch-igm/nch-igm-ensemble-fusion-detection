# Fusion Detection Overlap

## Sections

- [Introduction](#introduction)
- [Getting Started](#getting-started)
- [Run Overlap](#running-overlap-from-kickoff_overlap.sh-using-bash)

## Introduction

Fusion detection overlap merges fusion output data from [Arriba](https://github.com/suhrig/arriba), [CICERO](https://github.com/stjude/CICERO), [FusionMap](http://www.arrayserver.com/wiki/index.php?title=Oshell#OmicScript_for_FusionMap), [FusionCatcher](https://github.com/ndaniel/fusioncatcher), [JAFFA](https://github.com/Oshlack/JAFFA/wiki), [MapSplice](http://www.netlab.uky.edu/p/bioinfo/MapSplice2), [SOAPfuse](https://sourceforge.net/projects/soapfuse/), [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki), [TopHat-Fusion](https://ccb.jhu.edu/software/tophat/fusion_index.shtml), and [DRAGEN RNA-Seq](https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html).

## Getting Started

### Download code

Download the overlap code from GitHub. Use `git` to clone the repo or navigate to the [nch-igm-ensemble-fusion-detection](https://github.com/nch-igm/nch-igm-ensemble-fusion-detection) GitHub page, download, and unzip the code.

```bash
git clone https://github.com/nch-igm/nch-igm-ensemble-fusion-detection.git
```

### Install dependencies

The main dependencies are R [(R Setup)](#r-setup) and a Postgres database [(Postgres Database Setup)](#postgres-database-setup).

#### R Setup

Install R (version 3.5 is preffered since the code is tested on this version). Visit the [R website](https://www.r-project.org/) for installation instructions or use your preferred install tool.

##### R Packages

The following packages are required for fusion detection overlap.

Start an R session.

```R
R
```

Run these commands in the R session.

```R
install.packages('devtools', repos='http://cran.us.r-project.org/')
install.packages('optparse', repos='http://cran.us.r-project.org/')
install.packages('dplyr', repos='http://cran.us.r-project.org/')
install.packages('readr', repos='http://cran.us.r-project.org/')
install.packages('DBI', repos='http://cran.us.r-project.org/')
install.packages('RPostgres', repos='http://cran.us.r-project.org/')
install.packages('dbplyr', repos='http://cran.us.r-project.org/')
install.packages('magrittr', repos='http://cran.us.r-project.org/')
install.packages('purrr', repos='http://cran.us.r-project.org/')
```

#### Postgres Database Setup

Overlap interacts with a database of fusions. The database is used to store previously called fusions. The overlap code queries the database for the frequency of reported fusions. Fusion calls with a high frequency are regarded as artifacts or unlikely to cause disease.

##### Step 1: Install Postgres

On macOS, you can use [Homebrew](https://brew.sh/) to install Postgres.

```bash
brew install postgres
```

On Linux or Windows, you will need to find alternative instructions.

##### Step 2: (Initialize and) Start the database

On macOS, Homebrew automatically initializes the space where the databases will be written on your filesystem at `/usr/local/var/postgres`. Users may experience different results on different systems. You may need to initialize the database at `/usr/local/var/postgres` or another location or you may not need to initialize the database at all.

```bash
initdb /usr/local/var/postgres
```

Start the database.

```bash
pg_ctl -D /usr/local/var/postgres start
```

##### Step 3: Setup the fusions database with test data

Create the fusions database.

```bash
createdb fusions
```

Create the postgres user if it doesn't already exist.

```bash
createuser -s postgres
```

Import the schema and test data for the fusions database. You may need to modify the host.

```bash
cd nch-igm-ensemble-fusion-detection
psql -h localhost fusions < sql/fusions_example_db_dump.sql
```

Modify the read-only and read-write users' passwords in the create_users.sql file.

```bash
vim sql/create_users.sql  # set passwords
```

Import the read-only and read-write users.

```bash
psql -h localhost fusions < sql/create_users.sql
```

##### Step 4: Setup the `.dbconfig.R` file

`R/.dbconfig.R` sets the configuration variables that the overlap code uses to connect to the database. R variables `dbname`, `host`, `port`, `user`, and `passwd` must be set in this file so that the overlap code can connect to the database.

```bash
echo 'dbname <- "fusions"
host <- "localhost"
port <- 5432
user <- "fusion_user_rw"
passwd <- "PASSWORD_RW_HERE"' > R/.dbconfig.R
```

## Running Overlap from kickoff_overlap.sh using bash
#### In order for the overlap script to run properly, there is an expected file structure hierarchy, please see
#### the Example below for assemble_results or view example data in "test_data" directory

### kickoff_overlap.sh
Run 'kickoff_overlap.sh' to run the necessary R scripts (described in detail below)

```bash
overlap/kickoff_overlap.sh -h
```
The result:
```txt
[USAGE]: This script kicks off a script that merges fusion detection results from many tools.

-h * [None] Print this help message.
-p * [path1,path2,...] paths  where fusion results are where the script will be kicked off.
-s * [samples_file1,samples_file2,...] Files paths that contain the lists of samples to kickoff scripts for. The number of samples files provided must be equal to the number of paths, and they must be listed in the same respective order. If not specifically added, the script will search through sample names that are listed in the samples file within the path (-p) given[Default: samples,samples,...].
```
Note that -s is not required, and if not added the script will automatically run on all samples listed in your samples file


#### Try kickoff_overlap.sh with test data:

```bash
overlap/kickoff_overlap.sh -p /test_data
```

Results (will show you which callers it found output for):
```txt
Test
WARNING: ignoring environment value of R_HOME
Assembling a list of files to merge...
[1] "starFusion"    "fusionMap"     "fusionCatcher" "jaffa"        
[5] "mapSplice"     "dragen"       
```

#### Explanation of and option to run R code directly:

Run `assemble_results.R` to merge results and generate the overlap report. See the help message for `assemble_results.R` by running

```bash
R/assemble_results.R -h
```

or

```bash
Rscript R/assemble_results.R -h
```

The result

```txt
Usage: assemble_results.R --sample s1

assemble_results.R reads a set of fusion detection result files and produces an aggregated report with the overlapping fusion events predicted in the input files. The program searches recursively from the current directory to find known files. The input files can also be specified manually.

Options:
	--sample=SAMPLE
		Name of the sample. (required)

	--baseDir=BASEDIR
		Base directory to search for input files. [default = .]

	--outReport=OUTREPORT
		Location of the output report. [default = overlap_$sample.tsv]

	--collapseoutReport=COLLAPSEOUTREPORT
		Location of the output report. [default = collapsed_3callers_$sample.tsv]

	--foutReport=FOUTREPORT
		Location of the output report. [default = filtered_overlap_2callers_$sample.tsv]

	--foutReport3=FOUTREPORT3
		Location of the output report. [default = filtered_overlap_3callers_$sample.tsv]

	--outSingleton=OUTSINGLETON
		Location of the Singleton output report. [default = Singleton_KnownFusions_$sample.tsv]

	--outBreakpoints=OUTBREAKPOINTS
		Location of the breakpoint report. [default = breakpoints_$sample.tsv]

	--dragen=DRAGEN
		Path to the dragen results file. [default = search for file named like 'DRAGEN.fusion_candidates.final']

	--fusionCatcher=FUSIONCATCHER
		Path to the fusionCatcher results file. [default = search for file named like 'final-list_candidate-fusion-genes.txt']

	--fusionMap=FUSIONMAP
		Path to the fusionMap results file. [default = search for file named like 'FusionDetection.FusionReport.Table.txt']

	--jaffa=JAFFA
		Path to the jaffa results file. [default = search for file named like 'jaffa_results.csv']

	--mapSplice=MAPSPLICE
		Path to the mapSplice results file. [default = search for file named like 'fusions_well_annotated.txt']

	--soapFuse=SOAPFUSE
		Path to the soapFuse results file. [default = search for file named like '.final.Fusion.specific.for.genes']

	--starFusion=STARFUSION
		Path to the starFusion results file. [default = search for file named like 'star-fusion.fusion_predictions.abridged.tsv']

	--tophatFusion=TOPHATFUSION
		Path to the tophatFusion results file. [default = search for file named like 'result.txt']. potential_fusion.txt must be in the same directory as result.txt for the import to work.

	--arriba=ARRIBA
		Specific location of the arriba results file (fusions.tsv).

	--cicero=CICERO
		Specific location of the cicero results file (annotated.fusion.txt).

	-h, --help
		Show this help message and exit```

### `--sample`

`--sample` is the only required argument for `assemble_results.R`.

```bash
R/assemble_results.R --samples s1
```

`--sample` is used to insert a comment into the top of the output files

```txt
# Sample: s1
# NumToolsAggregated: 5
# - starFusionCalls = 20
# - fusionCatcherCalls = 552
# - jaffaCalls = 948
# - mapSpliceCalls = 23
# - dragenCalls = 9
```

and name the output files

```txt
--outReport      	= overlap_$sample.tsv
--collapsedoutReport	= collapsed_3callers_$sample.tsv
--foutReport     	= filtered_overlap_2callers_$sample.tsv
--foutReport3    	= filtered_overlap_3callers_$sample.tsv
--outSingleton		= Singleton_KnownFusions_$sample.tsv
--outBreakpoints 	= breakpoints_$sample.tsv
```

You can override the default output file names with the `--*out*` parameters.

### `--baseDir`

The `--baseDir` argument sets the location of the input data. By default, `--baseDir` is set to the current directory `.`, but you can override the default location

```bash
R/assemble_results.R --sample s1 --baseDir fusion_results
```

`kickoff_overlap.R` recursively searches all directories under the `baseDir` for fusion detection results files. The default search looks for

```txt
--arriba        = fusions.tsv
--cicero        = annotated.fusion.txt
--dragen        = DRAGEN.fusion_candidates.final
--fusionCatcher = final-list_candidate-fusion-genes.txt
--fusionMap     = FusionDetection.FusionReport.Table.txt
--jaffa         = jaffa_results.csv
--mapSplice     = fusions_well_annotated.txt
--soapFuse      = .final.Fusion.specific.for.genes
--starFusion    = star-fusion.fusion_predictions.abridged.tsv
--tophatFusion  = result.txt (potential_fusion.txt must be in the same directory as result.txt)
```

For each tool, overlap uses the first occurence of a file that matches the pattern associated with that tool. If there are no files that match for the tool, then overlap assumes that the tool was not run.

_Further info: `kickoff_overlap.R` uses the R function `list.files` with the `pattern` parameter set to the pattern. The pattern is used as a regular expression used to match file names. The `soapFuse` pattern `.final.Fusion.specific.for.genes` starts with a period because the actual output file is named `[sample].final.Fusion.specific.for.genes`. Our pattern will match `soapFuse` output files with different samples names, like `s1.final.Fusion.specific.for.genes` or `random_sample_name.final.Fusion.specific.for.genes`._

### Example

Given the following directory structure

```txt
/data/fusion_output_1/
    final-list_candidate-fusion-genes.txt
    FusionDetection.FusionReport.Table.txt
    star-fusion.fusion_predictions.abridged.tsv
    DRAGEN.fusion_candidates.final
```

run the command

```bash
cd /data/fusion_output_1
~/igm-ensemble-fusion-detection/R/assemble_results.R --sample s1
```

or

```bash
cd ~/igm-ensemble-fusion-detection
R/assemble_results.R --sample s1 --baseDir /data/fusion_output_1
```

and overlap will be able to find fusion output files for `fusionCatcher`, `fusionMap`, `starFusion`, and `dragen`. Overlap will be able to find the correct files even if they are under subdirectories since the search is recursive.

```txt
/data/fusion_output_1/
    fusioncatcher/
        final-list_candidate-fusion-genes.txt
    fusionmap/
        results/
            FusionDetection.FusionReport.Table.txt
    starfusion/
        star-fusion.fusion_predictions.abridged.tsv
    dragen/
        DRAGEN.fusion_candidates.final
```

### `--starFusion` override example

You can override the default fusion output file search with direct fusion output file names.

```bash
R/assemble_results.R --samples s1 --starFusion path/to/custom-star-fusion-output.tsv
```

This command will recursively search under the current directory `.` for the `dragen`, `fusionCatcher`, `fusionMap`, etc., default files, but it will not perform the same search for the `starFusion` file. It will use the exact location to the `starFusion` file as provided by the user `path/to/custom-star-fusion-output.tsv`.

The overlap code only _knows_ how to parse the data from the output files in the default list. Make sure your custom locations point to files that hold the same data. There may be multiple output files from each fusion caller. Our overlap code is only parsing certain files for certain data to generate the overlap report.

## Upload to the database

Run `kickoff_upload.sh` to upload fusion results to the database.

```bash
overlap/kickoff_upload.sh -h
```

The result

```txt
[USAGE]: Use this script to upload fusion detection results to the fusion detection database.

-h * [None] Print this help message.
-d * [results_dir1,results_dir2,...] List of directories (comma-separated) where the fusion results are located. Directories must be full paths. Ex: /path/to/fusion_results1,/path/to/fusion_results2.
-p * [patient1,patient2,...] List of patient IDs (comma-separated). The number of patient IDs must be equal to the number of results directories and they must be listed in the same respective order. Ex: patient1,patient2. (Also referred to as the 'subject ID'.)
-s * [sample1,sample2,...] List of sample IDs (comma-separated). The number of sample IDs must be equal to the number of results directories and they must be listed in the same respective order. Ex: sample1,sample2.
-r * [reads1,reads2,...] List of read pairs numbers (comma-separated). The number of read pairs numbers must be equal to the number of results directories and they must be listed in the same respective order. Optional. Ex: 40000000,52000000.

[Example 1]: kickoff_upload.sh -d /path/to/fusion_results1 -p p1 -s s1
[Example 2]: kickoff_upload.sh -d /path/to/fusion_results1,/path/to/fusion_results2 -p p1,p2 -s s1,s2
[Example 3]: kickoff_upload.sh -d /path/to/fusion_results1 -p p1 -s s1 -r 42543879
```

The `-d` results directory must contain fusion detection results files. The results directory usually is the top-level directory for a sample where multiple fusion caller results are stored in sub-directories. The upload code recursively searches all directories under the results directory for fusion detection results files. See the [example](#example) above for an example results directory layout. `kickoff_upload.sh` can take a single results directory as input (for results from multiple tools on a single sample) or multiple results directories as input (to make the upload process faster for multiple samples).

The `-p` patient ID and `-s` sample ID are necessary identifiers to upload the data to the database. The `-r` number of read pairs is optional, depending on whether you want to store that information.

### Analysis

The database uses the term 'analysis' to refer to a fusion caller execution on a single sample. A `STAR-Fusion` run on `sample1` is considered an analysis, say `analysis1`. A `JAFFA` run on `sample1` is considered a different analysis, say `analysis2`.

#### Analysis upload process

When uploading data to the database, the upload script uses the concept of an analysis to decide whether the data already exists in the database. If the analysis already exists in the database, then the upload script will skip the upload for that analysis. If there is already `STAR-Fusion` data for `sample1` in the database (`analysis1`), then any additional attempt to upload or replace `STAR-Fusion` data for `sample1` will fail.

- Upside - Since an analysis will not be uploaded twice, you can re-run the upload script on the same results directory and not have to worry about overwriting data. This can be useful if, for example, you run `STAR-Fusion` and `JAFFA` on `sample1` and upload the results. At a later time, you decide to run `FusionMap` on the same sample. When you want to add the `FusionMap` results to the database, you can run the upload script on the `sample1` results directory and it will only add the new `FusionMap` analysis to the database.
- Downside - You cannot replace analysis data using the upload script. If, for example, you run a new version of `STAR-Fusion` on `sample1`, you might want to replace the old `STAR-Fusion` results in the database. The upload script does not handle this scenario and you will have to perform a manual replacement.

### Upload with R instead of bash

Use `upload_fusion_results.R` for a greater level of control over the database upload (or if you don't have `bash` on your system to run `kickoff_upload.sh`).

```bash
R/upload_fusion_results.R -h
```

The result

```txt
Usage: R/upload_fusion_results.R [options]

This script acts as an importer for fusion detection results to the fusion detection database.

Options:
        --subject=SUBJECT
                Subject ID, e.g. subject1. (required)
        --sample=SAMPLE
                Sample ID, e.g. sample1. (required)
        --phenotype=PHENOTYPE
                Phenotype of the subject. (optional)
        --storage=STORAGE
                Storage method for the sample, such as 'Frozen' or 'FFPE'. (optional)
        --reads=READS
                Number of read pairs for the sample. (optional)
        --tool=TOOL
                Custom tool for which the results are being entered, possible values are 'dragen', 'FusionCatcher', 'FusionMap', 'Jaffa', 'MapSplice', 'SoapFuse', 'StarFusion', and 'TophatFusion'. (optional)
        --results=RESULTS
                Custom location of the results file for the associated tool, e.g. 'star-fusion.fusion_predictions.abridged.tsv'. (optional)
        -h, --help
                Show this help message and exit
```

`upload_fusion_results.R` must be run from the results directory. It searches for results files under `.`. Alternatively, you can upload results one-by-one by providing the `--tool` and `--results` parameters.
