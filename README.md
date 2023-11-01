# TSO500 v2.0

## What does this app do?
Runs the Illumina TSO500 local analysis app.

## What inputs are required for this app to run?

### Required
- `TSO500_ruo` (`file`) - a zip file (originally provided by Illumina) containing the TSO500 local analysis app.
- `input_files` (`array:files`) - raw sequencing data in tar files

### Optional
- `samplesheet` (`file`) - samplesheet to use for demultiplexing and analysis, will use the one in the run data if not provided
- `demultiplex_only` (`bool`) - flag to demultiplex a runfolder and upload the fastq and log files
- `upload_demultiplex_output` (`bool`; default: true): controls if to upload the demultiplexing output files
- `analysis_options` (`str`) -  a string which can be passed to the ruo command line
- `isNovaSeq` (`bool`; default: true) - passes the `-isNovaSeq` flag to the TSO500 local app for running on NovaSeq data
- `scatter_instance` (`str`): DNAnexus instance type to use for the per sample analysis (default: `mem1_ssd1_v2_x36`)
- `include_samples` (`str`) - comma separated string of samples to run analyses for (mutually exclusive with `exclude_samples`)
- `exclude_samples` (`str`) - comma separated string of samples to NOT run analyses for (mutually exclusive with `include_samples`)
- `n_samples` (`int`) - maximum number of samples from samplesheet to run analysis on (this will take the first n sample rows from the samplesheet)

## How does this app work?
The app runs the TSO500 local app in the 'scatter / gather' mode (explained on [page 9 here][user-guide]), this works by splitting off the per sample analysis into separate sub jobs in parallel and then combining the output in the parent job to produce the final results. This greatly speeds up analysis vs running the local app sequentially on all samples. The general outline is as follows:

- download and unpack the TSO500 resources zip file
- download and unpack the run data tar files
- parse samplesheet and (optionally) include/exclude sample rows if `-iinclude_samples` or `-iexclude_samples` specified
- run demultiplexing
- launch sub jobs per sample to run the initial analysis
- hold parent job until all scatter jobs complete
- gather up scatter job outputs and run final gather analysis step
- upload all final output files and set app output spec


## What does this app output

This outputs the full TSO500 local app, including all analysis files and intermediate logs. Some of the sample analysis files are first uploaded to specific fields in the output spec for the app, these include:

| Field           	    | Type  	| Path                                                    	| Notes                                                        	|
|---------------------  |-------	|---------------------------------------------------------	|--------------------------------------------------------------	|
| `fastqs`        	    | FASTQ 	| `/DemultiplexOutput/Logs_Intermediates/FastqGeneration` 	| The FASTQs as output from bclconvert                         	|
| `dna_bams`      	    | BAM   	| `/Analysis/*/Logs_Intermediates/StitchedRealigned/`     	| The 'final' BAM file upon which variant calling is performed 	|
| `dna_bam_index` 	    | BAI   	| `/Analysis/*/Logs_Intermediates/StitchedRealigned/`     	| Index of the above BAM                                       	|
| `rna_bam`       	    | BAM   	| `/Analysis/*/Logs_Intermediates/RnaAlignment/*`         	| Downsampled and trimmed BAM aligned by STAR aligner          	|
| `rna_bam_index` 	    | BAI   	| `/Analysis/*/Logs_Intermediates/RnaAlignment/*`         	| Index of the above BAM                                       	|
| `msi_metrics`   	    | JSON  	| `/Analysis/*/Logs_Intermediates/Msi/*`                  	| MSI metrics in JSON format                                   	|
| `tmb_metrics`   	    | JSON  	| `/Analysis/*/Logs_Intermediates/Tmb/*`                  	| TMB metrics in JSON format                                   	|
| `fusions`             | CSV       | `/Results/Results/*/*`                                    | CSV of all identified fusions from RNA analysis               |
| `small_variant_annotation`    	    | JSON  	| `/Results/Results/*/*`                                  	| Nirvana annotated variant JSON file                          	|
| `splice_variant_annotation`   |   JSON    | `/Results/Results/*/*`    |   Annotated splice variants from RNA analysis in JSON     |
| `cvo`           	    | TSV   	| `/Results/Results/*/*`                                  	| CombinedVariantOutput TSV file                               	|
| `metricsOutput` 	    | TSV   	| `/Results/Results/`                                     	| Final MetricsOutput TSV file for the run                     	|
| `cnv_vcfs`      	    | VCF   	| `/Results/Results/*/*`                                  	| Final CNV VCF from Results directory                         	|
| `gvcfs`         	    | VCF   	| `/Results/Results/*/*`                                  	| Genome VCF from the merging step                             	|
| `splice_variants_vcf` | VCF       | `/Results/Results/*/*`                                    | VCF of splice variants from RNA analysis                      |



Other demultiplexing files (i.e logs etc) are uploaded and attributed to the `demultiplex_logs` output field, and all other analysis files are uploaded and attributed to the `analysis_folder` output field.

The general output directory structure is (n.b. per sample directories / files have been limited to one for brevity):
```
output_folder/
├── scatter
│   ├── Annotation
│   │   └── sample1
│   │       ├── sample1_SmallVariants_Annotated.json.gz
│   │       └── sample1_SmallVariants_Annotated.json.gz.jsi
│   ├── DnaQCMetrics
│   │   └── sample1
│   │       ├── sample1.aligned.metrics.json
│   │       ├── sample1.cnv.metrics.json
│   │       ├── sample1.collapsed.metrics.json
│   │       └── sample1.stitched.metrics.json
│   ├── Contamination
│   │   └── sample1
│   │       └── sample1.contamination.json
│   ├── StitchedRealigned
│   │   └── sample1
│   │       ├── GeminiMultiLogs
│   │       │   └── GeminiMultiOptions.used.json
│   │       ├── sample1.bam
│   │       └── sample1.bam.bai
│   ├── VariantMatching
│   │   └── sample1
│   │       └── sample1_MergedSmallVariants.genome.vcf.gz
│   ├── PhasedVariants
│   │   └── sample1
│   │       ├── PsaraLogs
│   │       │   └── PsaraOptions.used.json
│   │       ├── ScyllaLogs
│   │       │   └── ScyllaOptions.used.json
│   │       ├── sample1.Complex.vcf.gz
│   │       └── sample1_SmallVariants.filtered.genome.vcf.gz
│   ├── Msi
│   │   └── sample1
│   │       ├── sample1.msi.json
│   │       ├── diffs.txt
│   │       └── tumor.dist
│   ├── DnaRealignment
│   │   └── sample1
│   │       ├── sample1.bam
│   │       └── sample1.bam.bai
│   ├── VariantCaller
│   │   ├── sample1
│   │   │   ├── sample1.filtered.genome.vcf.gz
│   │   │   └── sample1.genome.vcf.gz
│   │   │   └── PiscesOptions.used.json
│   │       └── PsaraOptions.used.json
│   ├── SmallVariantFilter
│   │   └── sample1
│   │       ├── sample1_SmallVariants.genome.vcf.errorRates
│   │       └── sample1_SmallVariants.genome.vcf.gz
│   ├── CollapsedReads
│   │   └── sample1
│   │       ├── sample1_metrics.json
│   │       ├── sample1_R1.fastq.gz
│   │       └── sample1_R2.fastq.gz
│   ├── CnvCaller
│   │   └── sample1
│   │       ├── sample1_CopyNumberVariants.vcf.gz
│   │       ├── sample1_foldChange.tsv
│   │       ├── sample1_normalizedBinCount.tsv
│   │       ├── sample1_rawBinCount.tsv
│   │       └── Craft_6735_20231026_113343.txt
│   ├── SamplesheetValidation
|   |   └── sample1
│   │       └── 20231026_224827_SampleSheet.csv
│   ├── DnaAlignment
│   │   └── sample1
│   │       ├── sample1.bam
│   │       └── sample1.bam.bai
│   ├── Tmb
│   │   └── sample1
│   │       ├── sample1.tmb.json
│   │       └── sample1_TMB_Trace.tsv
│   ├── CombinedVariantOutput
│   │   └── sample1
│   │       └── sample1_CombinedVariantOutput.tsv
│   ├── SampleAnalysisResults
|   |   └── sample1
│   │       └── sample1_SampleAnalysisResults.json
│   ├── RnaAlignment
|   |   └── sample1
│   │       ├── starLoad.Aligned.out.sam
│   │       └── starUnload.Aligned.out.sam
│   ├── MetricsOutput
|   |   └── sample1
│   │       └── MetricsOutput.tsv
│   |── MergedAnnotation
│   │   └── sample1
│   |       ├── sample1_MergedVariants_Annotated.json.gz
│   |       └── sample1_MergedVariants_Annotated.json.gz.jsi
│   |── inputs.json
│   └── SampleSheet.csv
|
├── gather
│   ├── Logs_Intermediates
│   │   ├── CombinedVariantOutput
│   │   │   ├── sample1
│   │   │   │   └── sample1_CombinedVariantOutput.tsv
│   │   │   ├── sample1
│   │   │   │   └── sample1_CombinedVariantOutput.tsv
...
│   │   ├── Msi
│   │   │   ├── 125501523-23269S0032-23TSOD58-8471
│   │   │   │   └── 125501523-23269S0032-23TSOD58-8471.msi.json
│   │   │   ├── 125510260-23269S0002-23TSOD58-8471
│   │   │   │   └── 125510260-23269S0002-23TSOD58-8471.msi.json
...
│   │   ├── Contamination
│   │   │   ├── sample1
│   │   │   │   └── sample1.contamination.json
│   │   │   ├── sample1
...
│   │   ├── DnaQCMetrics
│   │   │   ├── sample1
│   │   │   │   ├── sample1.aligned.metrics.json
│   │   │   │   ├── sample1.cnv.metrics.json
│   │   │   │   ├── sample1.collapsed.metrics.json
│   │   │   │   └── sample1.stitched.metrics.json
...
│   │   ├── CollapsedReads
│   │   │   ├── sample1
│   │   │   │   └── sample1_metrics.json
...
│   │   ├── SampleAnalysisResults
│   │   │   ├── sample1_SampleAnalysisResults.json
...
│   │   ├── MetricsOutput
│   │   │   └── MetricsOutput.tsv
│   │   ├── Annotation
│   │   │   ├── sample1
│   │   │   │   └── sample1_SmallVariants_Annotated.json.gz
...
│   │   ├── MergedAnnotation
│   │   │   ├── sample1
│   │   │   │   └── sample1_MergedVariants_Annotated.json.gz
...
│   │   ├── RunQc
│   │   │   └── RunQCMetrics.json
│   │   └── Tmb
│   │       ├── sample1
│   │       │   ├── sample1.tmb.json
│   │       │   └── sample1_TMB_Trace.tsv
...
│   ├── Results
│   │   ├── sample1
│   │   │   ├── sample1_CombinedVariantOutput.tsv
│   │   │   ├── sample1_CopyNumberVariants.vcf
│   │   │   ├── sample1_MergedSmallVariants.genome.vcf
│   │   │   ├── sample1_MergedVariants_Annotated.json.gz
│   │   │   └── sample1_TMB_Trace.tsv
...
│   ├── inputs.json
│   └── SampleSheet.csv
├── demultiplexOutput
│   ├── Logs_Intermediates
│   │   ├── FastqGeneration
│   │   │   ├── Reports
│   │   │   │   ├── Lane_1
│   │   │   │   │   ├── Adapter_Metrics.csv
│   │   │   │   │   ├── Demultiplex_Stats.csv
│   │   │   │   │   ├── fastq_list.csv
│   │   │   │   │   ├── Index_Hopping_Counts.csv
│   │   │   │   │   ├── IndexMetricsOut.bin
│   │   │   │   │   ├── RunInfo.xml
│   │   │   │   │   ├── SampleSheet.csv
│   │   │   │   │   └── Top_Unknown_Barcodes.csv
│   │   │   │   └── Lane_2
│   │   │   │       ├── Adapter_Metrics.csv
│   │   │   │       ├── Demultiplex_Stats.csv
│   │   │   │       ├── fastq_list.csv
│   │   │   │       ├── Index_Hopping_Counts.csv
│   │   │   │       ├── IndexMetricsOut.bin
│   │   │   │       ├── RunInfo.xml
│   │   │   │       ├── SampleSheet.csv
│   │   │   │       └── Top_Unknown_Barcodes.csv
│   │   │   ├── Logs
│   │   │   │   ├── Lane_2
│   │   │   │   │   ├── Errors.log
│   │   │   │   │   ├── FastqComplete.txt
│   │   │   │   │   ├── Info.log
│   │   │   │   │   └── Warnings.log
│   │   │   │   └── Lane_1
│   │   │   │       ├── Errors.log
│   │   │   │       ├── FastqComplete.txt
│   │   │   │       ├── Info.log
│   │   │   │       └── Warnings.log
│   │   │   ├── sample1
│   │   │   │   ├── sample1_S2_L001_R1_001.fastq.gz
│   │   │   │   ├── sample1_S2_L001_R2_001.fastq.gz
│   │   │   │   ├── sample1_S2_L002_R1_001.fastq.gz
│   │   │   │   └── sample1_S2_L002_R2_001.fastq.gz
...
│   │   │   ├── dsdm.json
│   │   │   ├── FastqGeneration-150369-20231026-223016.stderr
│   │   │   ├── FastqGeneration-150369-20231026-223016.stdout
│   │   │   ├── FastqGeneration-20231026222542.log
│   │   │   ├── FastqGeneration-417-20231026-222542.stderr
│   │   │   ├── FastqGeneration-417-20231026-222542.stdout
│   │   │   └── SampleSheet_combined.csv
│   │   ├── SamplesheetValidation
│   │   │   ├── 20231026_222415_SampleSheet.csv
│   │   │   ├── dsdm.json
│   │   │   └── SamplesheetValidation-20231026222415.log
│   │   ├── RunQc
│   │   │   ├── dsdm.json
│   │   │   ├── RunQc-20231026222537.log
│   │   │   └── RunQCMetrics.json
│   │   └── ResourceVerification
│   │       ├── dsdm.json
│   │       └── ResourceVerification-20231026222421.log
│   ├── demultiplex_cromwell_executions.tar.gz
│   ├── inputs.json
│   ├── receipt
│   ├── SampleSheet.csv
│   └── trusight-oncology-500-ruo_ruo-2.2.0.12_20231026-232354.log
├── logs
│   ├── sample1_cromwell_executions.tar.gz
│   ├── sample1_scatter_logs.tar.gz
...
│   ├── analysis_results_logs.tar.gz
│   └── gather_cromwell_executions.tar.gz
└── MetricsOutput.tsv
```


## Notes
- Samplesheet input is optional, if not specified the analysis app looks for the samplesheet in top level of runfolder
- When running in scatter / gather mode, demultiplexing via the local app must always be first performed (i.e it can't be started from previously demultiplexed data and reusing fastqs) To achieve the equivalent, `-iupload_demultiplex_output=false` may be specified to not upload the output of demultiplexing from the job. When used in conjunction with `-iinclude_samples="sample1"`, this would effectively just run and output data for the given sample(s) as if the local app were just run from fastqs for a single sample
- Jobs are launched per sample parsed from the 'Pair_ID' column, this means if the samplesheet is formatted for running in paired analysis mode, the fastqs for all samples for the given pair ID will be used for the analysis
- samples specified to `-iinclude_samples` or `-iexclude_samples` should be specified as given in the Pair_ID column of the samplesheet
- Intermediary genome vcfs found in `scatter/` are compressed before uploading to save on storage
- All log files are gathered up before uploading and combined into tar files, there is one tar file per file from the scatter step and one from the gather step. This is to reduce the total number of files uploaded at the end.


## This app was created by East Genomics GLH

[user-guide]: https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/trusight/trusight-oncology-500/trusight-oncology-500-local-app-v2.2-user-guide-1000000137777-01.pdf