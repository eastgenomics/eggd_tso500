# TSO500 v1.2

## What does this app do?
Runs the Illumina TSO500 local analysis app.

## What inputs are required for this app to run?
* TSO500_ruo - a zip file (originally provided by Illumina) containing the TSO500 local analysis app.
* input_files - Either a single tar.gz runfolder, a set of tars from dx-streaming-upload or individual FastQs
* (Optional) Samplesheet for the run if not present in the tarred inputs
* isFastQ - flag to allow the app to start without demultiplexing (Must be used with a SampleSheet input)
* demultiplexOnly - Flag to only demultiplex the runfolder (cannot be used in conjunction with --isFastQ)
* analysis_options -  a string which can be passed to the ruo command line

## How does this app work?
* Downloads and unzips/untars the TSO500 local analysis app and runfolder
* Runs the TruSight_Oncology_500_RUO.sh (within the TSO500 local app zip file) providing arguments for analysis folder, runfolder/fastqfolder, samplesheet, resourcesFolder and any other analysis options given as an input ($analysis_options)

## What does this app output
* RUO_stdout.txt - STDout from RUO. Saved into /logs
* The analysis folder. Saved into /analysis_folder - this includes the Results folder containing the final results and the Log Intermediates folder containing all files generated during analysis.
* More information please refer to the local app documentation: https://emea.support.illumina.com/downloads/trusight-oncology-500-local-app-user-guide-1000000067616.html

## Notes
* Samplesheet input is optional, if not specified the analysis app looks for SampleSheet.csv in top level of runfolder, the app checks to find the SampleSheet and rename it so it can be automatically found but it won't find it if it doesn't end in *heet.csv and not given as input.
* If the app fails it will print the SampleSheet Validation to identify if the error for failure was due to a invalid SampleSheet

## This app was created by East Genomics GLH
