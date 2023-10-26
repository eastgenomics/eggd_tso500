#!/bin/bash

# prefixes all lines of commands written to stdout with datetime
PS4='\000[$(date)]\011'
export TZ=Europe/London
set -exo pipefail

# set frequency of instance usage in logs to 60 seconds
kill $(ps aux | grep pcp-dstat | head -n1 | awk '{print $2}')
/usr/bin/dx-dstat 60


_get_tso_resources() {
    : '''
    Download and unpack the TSO500 resources zip file
    '''
    sudo dpkg -i ripunzip_0.4.0_amd64.deb

    local zip_file=$(dx describe --json "$TSO500_ruo" | jq -r '.name')
    echo "Downloading and unpacking TSO500 resources (${zip_file})"

    SECONDS=0
    dx download "$TSO500_ruo"
    ripunzip -d TSO500_ruo file "$zip_file"
    chmod -R 777 TSO500_ruo/TSO500_RUO_LocalApp/
    sudo docker load --input TSO500_ruo/TSO500_RUO_LocalApp/trusight-oncology-500-ruo-dockerimage-ruo-*.tar

    duration=$SECONDS
    echo "Downloading and unpacking ${zip_file} took $(($duration / 60))m$(($duration % 60))s"
}


_get_input_files() {
    : '''
    Download all input tar files

    Expect input files to be  tars of the run data (i.e. bcl files) =>
    download the runfolder input, decompress and save in directory 'runfolder'
    '''
    SECONDS=0
    echo "Downloading tar files"
    
    # drop the $dnanexus_link from the file IDs
    file_ids=$(grep -Po  "file-[\d\w]+" <<< "${input_files[@]}")

    echo "$file_ids" | xargs -P ${THREADS} -n1 -I{} sh -c \
        "dx cat {} | tar -I pigz -xf - --no-same-owner --absolute-names -C /home/dnanexus/runfolder"

    total=$(du -sh /home/dnanexus/runfolder | cut -f1)

    duration=$SECONDS
    echo "Downloaded $(wc -w <<< ${file_ids}) files (${total}) in $(($duration / 60))m$(($duration % 60))s"
}


_get_scatter_job_outputs() {
    : '''
    Download output from all scatter jobs to run final gather step

    All scatter job IDs were logged to job_ids file, so we will read
    from here and use the ids to filter down the files in the output
    folder to ensure we keep the same directory structure
    '''
    SECONDS=0
    echo "Downloading scatter job output"

    set +x  # suppress this going to the logs as its long

    # files from sub jobs will be in the container- project context of the
    # current job ($DX_WORKSPACE-id) => search here for  all the files
    scatter_files=$(dx find data --json --verbose --path "$DX_WORKSPACE_ID:/Analysis")
    scatter_files+=" $(dx find data --json --verbose --path "$DX_WORKSPACE_ID:/Logs")"

    # turn describe output into id:/path/to/file to download with dir structure
    files=$(jq -r '.[] | .id + ":" + .describe.folder + "/" + .describe.name'  <<< $scatter_files)

    # build aggregated directory structure and download all files
    cmds=$(for f in  $files; do \
        id=${f%:*}; path=${f##*:}; dir=$(dirname "$path"); \
        echo "'mkdir -p out/$dir && dx download --no-progress $id -o out/$path'"; done)

    echo $cmds | xargs -n1 -P${THREADS} bash -c

    set -x

    total=$(du -sh /home/dnanexus/out/Analysis/ | cut -f1)
    duration=$SECONDS
    echo "Downloaded $(wc -w <<< ${files}) files (${total}) in $(($duration / 60))m$(($duration % 60))s"
}


_get_samplesheet() {
    : '''
    Download samplesheet either from input or find it from the run data
    '''
    if [[ "$samplesheet" ]]; then
        # using provided samplesheet
        echo "Using provided sample sheet: ${samplesheet_name}"
        dx download "$samplesheet" -o runfolder/SampleSheet.csv
    elif [[ $(find ./ -regextype posix-extended  -iregex '.*sample[-_ ]?sheet.csv$') ]]; then
        # Sample sheet not given, try finding it in the run folder
        # Use regex to account for anything named differently
        # e.g. run-id_SampleSheet.csv, sample_sheet.csv, Sample Sheet.csv, sampleSheet.csv etc.
        local_samplesheet=$(find ./ -regextype posix-extended  -iregex '.*sample[-_ ]?sheet.csv$')
        echo "Using sample sheet in run directory: $local_samplesheet"
        mv "$local_samplesheet" /home/dnanexus/runfolder/SampleSheet.csv

        # upload to get a dx file ID to pass to downstream scatter jobs
        samplesheet=$(dx upload --brief "$local_samplesheet")
    else
        dx-jobutil-report-error "No SampleSheet could be found."
        exit 1
    fi
}


_parse_samplesheet() {
    : '''
    Parse the sample names from the samplesheet into variable "sample_list"

    If specified, can also optionally:
        - exclude specified samples (-iexclude_samples)
        - filter samplesheet to specific samples (-iinclude_samples)
        - limit samplesheet to first n samples (-in_samples)
    '''
    # parse list of sample names from samplesheet
    sample_rows=$(sed -e '1,/Sample_ID/ d' runfolder/SampleSheet.csv)
    sample_list=$(cut -d, -f1 <<< "$sample_rows")
    echo "Samples parsed from samplesheet: ${sample_list}"
}


_modify_samplesheet() {
    : '''
    If specified, will modify the samplesheet in the following way to change
    what samples analysis is run for:

        - exclude specified samples (-iexclude_samples)
        - filter samplesheet to specific samples (-iinclude_samples)
        - limit samplesheet to first n samples (-in_samples)
    
    The modified samplesheet will be uploaded to the output
    folder as "modified_SampleSheet.csv" to be able to pass
    the file ID to the scatter job.

    This modifies the following global variables:

        - samplesheet : str of file ID of samplesheet
        - sample_list : str of all sample names parsed from samplesheet
    '''
    echo "Modifying samplesheet"
    # parse original rows from samplesheet
    samplesheet_header=$(sed -n '1,/Sample_ID/ p' runfolder/SampleSheet.csv)
    sample_rows=$(sed -e '1,/Sample_ID/ d' runfolder/SampleSheet.csv)
    sample_list=$(cut -d, -f1 <<< "$sample_rows")
    echo "Original samples parsed from samplesheet: ${sample_list}"

    if [[ "$n_samples" ]]; then
        echo "Limiting analysis to ${n_samples} samples"
        sample_rows=$(sed "2,${n_samples}p;d" <<< "${sample_rows}")
        echo "Samples to run analysis for: ${sample_list}"
    fi

    if [[ "$include_samples" ]] || [[ "$exclude_samples" ]]; then
        # first ensure any samples specified to include or exclude are
        # a valid sample names from our samplesheet
        samples=$(sed -e 's/^,\|,$//' <<< "${include_samples},${exclude_samples}")

        invalid=$(while read -d',' -r line; do \
            [[ "$sample_list" =~ $line ]] || echo "${line} "; done <<< "${samples},")

        if [[ "${invalid}" ]]; then
            echo "One or more samplenames provided to include/exclude are invalid: ${invalid}"
            dx-jobutil-report-error "Invalid samplename specified"
        fi
    fi

    if [[ "$include_samples" ]]; then
        # retaining rows containing only those specified for given samples
        echo "-iinclude_samples specified: ${include_samples}"
        include=$(sed 's/,/|/g' <<< "$include_samples")
        sample_rows=$(awk '/'"$include"'/ {print $1}' <<< "$sample_rows")
        sample_list=$(cut -d, -f1 <<< "$sample_rows")
    fi

    if [[ "$exclude_samples" ]]; then
        # exclude rows containing only those specified for given samples
        echo -e "-iexclude_samples specified: ${exclude_samples}"
        exclude=$(sed 's/,/|/g' <<< "$exclude_samples")
        sample_rows=$(awk '!/'"$exclude"'/ {print $1}' <<< "$sample_rows")
        sample_list=$(cut -d, -f1 <<< "$sample_rows")
    fi

    # write out new samplesheet with specified rows
    echo "$samplesheet_header" >> runfolder/modified_SampleSheet.csv
    echo "$sample_rows" >> runfolder/modified_SampleSheet.csv

    # overwrite global variable of samplesheet file ID with new modified
    # file to be passed to scatter jobs
    samplesheet=$(dx upload --brief runfolder/modified_SampleSheet.csv)

    echo -e "New samplesheet:\n$(cat runfolder/modified_SampleSheet.csv)"
}


_calculate_total_fq_size() {
    : '''
    Calculate total fastq size per sample, used to dump into logs
    to give an indication of how long sub jobs will take

    Will output something like:

        125635958-23277S0009-23TSOD62-8471 1.14GB
        125636529-23277S0023-23TSOD62-8471 1.22GB
        ...
        125636099-23277S0012-23TSOD62-8471 3.22GB
	    125635945-23271S0014-23TSOD62-8471 3.27GB
    '''
    echo "Calculating total fastq sizes"
    set +x  # suppress this going to the logs as its long
    local fastqs=$1
    describe=$(xargs -P16 -n1 dx describe --json <<< $fastqs)
    sizes=$(echo $describe \
        | jq -r '[.name,.size] | @tsv' \
        | sed  -E 's/_S[0-9]+_L00[0-9]{1}_R[12]_[0-9]+\.fastq\.gz//g' \
        | awk '{ sum[$1] += $2 } END { for (i in sum) printf "\t%s %.2fGB\n", i, sum[i]/1024/1024/1024 }' \
        | sort -k2 -n)
    set -x
    
    echo -e "Total fastq sizes per sample:\n${sizes}"
}


_upload_single_file() {
  : '''
  Uploads single file with dx upload and associates uploaded
  file ID to specified output field

  Arguments
  ---------
    1 : str
        path and file to upload
    2 : str
        app output field to link the uploaded file to
    3 : bool
        (optional) controls if to link output file to job output spec
  '''
  local file=$1
  local field=$2
  local link=$3

  local remote_path=$(sed s'/\/home\/dnanexus\/out\///' <<< "$file")

  file_id=$(dx upload "$file" --path "$remote_path" --parents --brief)

  if [[ "$link" == true ]]; then 
    dx-jobutil-add-output "$field" "$file_id" --array
  fi

  echo "Uploaded ${remote_path}"
}


_upload_scatter_output() {
    : '''
    Upload the per sample output files in a scatter job

    To speed up the upload, logs and stdout/stderr files and  all cromwell
    logs in cromwell-executions/ are also combined into a tar files.
    '''
    SECONDS=0
    echo "Uploading sample output"
    export -f _upload_single_file

    # tar up all cromwell logs for faster upload
    tar -I pigz -cf /home/dnanexus/out/Logs/${sample}_cromwell_executions.tar.gz \
        /home/dnanexus/out/Analysis/${sample}_output/cromwell-executions
    rm -rf /home/dnanexus/out/Analysis/${sample}_output/cromwell-executions/

    # tar up all the logs, stdout and stderr for faster upload
    logs=$(find out/Analysis/${sample}_output/ \
        -name "*.log" -o \
        -name "*.out" -o \
        -name "*.stderr" -o \
        -name "*.stdout" -o\
        -name "receipt")

    mkdir all_logs && xargs -I{} mv {} all_logs/ <<< "$logs"
    tar -czf ${sample}_scatter_logs.tar.gz all_logs
    mv ${sample}_scatter_logs.tar.gz /home/dnanexus/out/Logs/
    mv /home/dnanexus/${sample}_scatter_stdout.txt /home/dnanexus/out/Logs/

    # compress intermediate genome VCFs since we don't use these routinely
    # and they go from >300mb to < 10mb (plus its a vcf, it should be compressed)
    find "/home/dnanexus/out/Analysis/${sample}_output/Logs_Intermediates/" -type f \
        -name "*.vcf"  -exec gzip {} \;

    # upload rest of files
    find /home/dnanexus/out/ -type f | xargs -P ${THREADS} -n1 -I{} bash -c \
        "_upload_single_file {} analysis_folder false"
    
    duration="$SECONDS"
    echo "Uploading completed in ${sample} in $(($duration / 60))m$(($duration % 60))s"
}


_upload_gather_output() {
    : '''
    Upload the final output files

    To generate distinct output fields from the app for result files, the
    following files are uploaded to individual output fields for each sample:

        - DNA bam & bai (from StitchedRealigned step)
        - RNA bam & bai (from RnaAlignment step)
        - MSI JSON (from Msi step)
        - TMB JSON (from Tmb step)
        - MergedSmallVariants.genome.vcf (from Results)
        - annotated JSON (Nirvana annotated variants; from Results)
        - CombinedVariantOutput (from Results)
        - CopyNumberVariants.vcf (from Results)
    
    The rest of the per sample analysis and gather step files are then uploaded
    in parallel to the "analysis_folder" output field.

    '''
    SECONDS=0
    export -f _upload_single_file

    # delete the per sample Results directories since these are just copies of
    # the files in the Logs_Intermediates sub directories and are also
    # copied into the gather step final Results directory
    find /home/dnanexus/out/Analysis -maxdepth 2 -type d -name "Results" -exec rm -rf {} \;

    echo "Total files to upload: $(find /home/dnanexus/out/Results /home/dnanexus/out/Analysis -type f | wc -l)"

    # tar up all cromwell logs for faster upload
    tar -I pigz -cf /home/dnanexus/out/Logs/gather_cromwell_executions.tar.gz \
        /home/dnanexus/out/Results/cromwell-executions
    rm -rf /home/dnanexus/out/Results/cromwell-executions

    # tar up all the logs, stdout and stderr for faster upload
    logs=$(find out/Analysis out/Results \
        -name "*.log" -o \
        -name "*.out" -o \
        -name "*.stderr" -o \
        -name "*.stdout" -o \
        -name "*script*" -o \
        -name "*dsdm.json" -o \
        -name "rc" -o \
        -name "receipt" -o \
        -name "*Log*.txt" -o \
        -name "*Log*.zip"
        )

    mkdir all_logs && xargs -I{} mv {} all_logs/ <<< "$logs"
    tar -czf analysis_results_logs.tar.gz all_logs
    mv analysis_results_logs.tar.gz /home/dnanexus/out/Logs

    # mapping of paths to search for files with file patterns and output field to attribute to
    # this allows us to make distinct groups of outputs for each file type, whilst
    # retaining the full directory structure once uploaded to the project
    file_mapping=(
        "*/Analysis/*/Logs_Intermediates/StitchedRealigned/*.bam dna_bams"
        "*/Analysis/*/Logs_Intermediates/StitchedRealigned/*.bai dna_bam_index"
        "*/Analysis/*/Logs_Intermediates/RnaAlignment/*.bam rna_bams"
        "*/Analysis/*/Logs_Intermediates/RnaAlignment/*.bai rna_bam_index"
        "*/Analysis/*/Logs_Intermediates/Msi/*.msi.json msi_metrics"
        "*/Analysis/*/Logs_Intermediates/Tmb/*.tmb.json tmb_metrics"
        "*/Results/Results/*/*MergedSmallVariants.genome.vcf gvcfs"
        "*/Results/Results/*/*Annotated.json.gz annotation"
        "*/Results/Results/*/*CombinedVariantOutput.tsv cvo"
        "*/Results/Results/*/*CopyNumberVariants.vcf cnv_vcfs"
    )

    # upload each of the sets of distinct output files
    for mapping in "${file_mapping[@]}"; do
        read -r path field <<< "$mapping"

        files=$(find /home/dnanexus/out -type f -path "$path")
        xargs -P ${THREADS} -n1 -I{} bash -c "_upload_single_file {} ${field} true" <<< "${files}"
        xargs -n1 -I{} mv {} /tmp <<< $files
    done

    # upload final run level MetricsOutput.tsv as distinct output field
    metrics_file_id=$(dx upload -p /home/dnanexus/out/Results/Results/MetricsOutput.tsv --brief)
    dx-jobutil-add-output "metricsOutput" "$metrics_file_id" --class=file
    mv /home/dnanexus/out/Results/Results/MetricsOutput.tsv /tmp

    # upload rest of files
    find /home/dnanexus/out/Analysis \
         /home/dnanexus/out/Logs \
         /home/dnanexus/out/Results -type f \
        | xargs -P ${THREADS} -n1 -I{} bash -c "_upload_single_file {} analysis_folder true"

    duration=$SECONDS
    echo "Uploading took $(($duration / 60))m$(($duration % 60))s."
}


_upload_demultiplex_output() {
    : '''
    Upload output of demultiplexing into project, this will be called
    before launching per sample scatter jobs since the fastqs need to
    be provided as input to the scatter sub job
    '''
    echo "Uploading demultiplexing output"
    SECONDS=0
    export -f _upload_single_file

    # first upload fastqs to set to distinct fastqs output field,
    # then upload the rest
    fastqs=$(find "/home/dnanexus/out/DemultiplexOutput" -type f -name "*fastq.gz")
    xargs -P ${THREADS} -n1 -I{} bash -c "_upload_single_file {} fastqs true" <<< "$fastqs"
    xargs -n1 -I{} mv {} /tmp <<< $fastqs

    # tar up all cromwell logs for faster upload
    tar -I pigz -cf /home/dnanexus/out/DemultiplexOutput/demultiplex_cromwell_executions.tar.gz \
        /home/dnanexus/out/DemultiplexOutput/cromwell-executions
    rm -rf /home/dnanexus/out/DemultiplexOutput/cromwell-executions

    # upload rest of files
    find "/home/dnanexus/out/DemultiplexOutput/" -type f | xargs -P ${THREADS} -n1 -I{} bash -c \
        "_upload_single_file {} demultiplex_logs true"

    duration=$SECONDS
    echo "Uploading took $(($duration / 60))m$(($duration % 60))s."

    # get the file IDs of the fastqs we have just uploaded from demultiplexing to
    # provide to the scatter jobs
    fastq_ids=$(cat job_output.json | jq -r '.fastqs[][]')

    if [[ "$upload_demultiplex_output" == false ]]; then
        # option specified to not upload demultiplex output => remove the job_output file
        # so that these are removed from the jobs output spec
        rm job_output.json
    fi
}


_demultiplex() {
    : '''
    Run demultiplexing via bclconvert in the TSO500 local app
    '''
    echo "Beginning demultiplexing"

    if [[ -f runfolder/modified_SampleSheet.csv ]]; then
        # include/exclude specified => use samplesheet that
        # we have modified
        samplesheet_path="runfolder/modified_SampleSheet.csv"
    else
        samplesheet_path="runfolder/SampleSheet.csv"
    fi

    SECONDS=0

    /usr/bin/time -v sudo bash TSO500_ruo/TSO500_RUO_LocalApp/TruSight_Oncology_500_RUO.sh \
        --analysisFolder /home/dnanexus/out/DemultiplexOutput/ \
        --runFolder /home/dnanexus/runfolder/ \
        --resourcesFolder /home/dnanexus/TSO500_ruo/TSO500_RUO_LocalApp/resources \
        --sampleSheet /home/dnanexus/$samplesheet_path \
        --demultiplexOnly \
        --isNovaSeq
    
    duration=$SECONDS

    echo "Demultiplexing completed in $(($duration / 60))m$(($duration % 60))s"
}


_scatter() {
    : '''
    Run the initial per sample analysis from the TSO500 local app

    This function will be launched as a sub job within the scope of
    the main TSO500 job, it sets up the job environment, runs the
    local app and uploads output to the project to then continue
    on to the gather step for the final results output

    Inputs
    ---------
    sample : str
        name of sample to run single analysis for
    samplesheet : str
        file ID of samplesheet
    fastqs : str
        str of all fastq file IDs
    options : str
        str of other cmd line options to specify to local app
    '''
    # prefixes all lines of commands written to stdout with datetime
    PS4='\000[$(date)]\011'
    export TZ=Europe/London
    set -exo pipefail

    # set frequency of instance usage in logs to 30 seconds
    kill $(ps aux | grep pcp-dstat | head -n1 | awk '{print $2}')
    /usr/bin/dx-dstat 30

    # control how many operations to open in parallel for download / upload
    THREADS=$(nproc --all)

    mkdir -p /home/dnanexus/runfolder \
             /home/dnanexus/TSO500_ruo \
             /home/dnanexus/DemultiplexOutput/$sample \
             /home/dnanexus/out/Logs \
             /home/dnanexus/out/Analysis

    dx download "$samplesheet" -o "SampleSheet.csv"

    # download the sample fastqs into directory for local app
    SECONDS=0
    echo "Finding fastqs for sample ${sample} from fastqs provided to job"
    set +x  # suppress this going to the logs as its long
    details=$(xargs -n1 -P${THREADS} dx describe --json --verbose <<< $fastqs)
    names=$(jq -r '.name' <<< $details)
    sample_fqs=$(jq -r "select(.name | startswith(\"${sample}_\")) | .id" <<< $details)
    set -x

    if [[ -z $sample_fqs ]]; then
        echo "ERROR: no fastqs found for sample ${sample} from provided fastqs:"
        sed 's/ /\n/g' <<< $names  # output to separate lines for nicer log viewing
        exit 1
    fi

    echo "Total fastqs passed: $(wc -w <<< $fastqs)"
    echo "Total fastqs found for sample: $(wc -w <<< $sample_fqs)"
    echo "Sample fastqs parsed: ${sample_fqs}"
    echo $sample_fqs | xargs -n1 -P${THREADS} -I{} sh -c \
        "dx download --no-progress -o /home/dnanexus/DemultiplexOutput/$sample/ {}"
    
    duration=$SECONDS
    echo "Downloaded fastqs in $(($duration / 60))m$(($duration % 60))s"
    
    # download and unpack local app resources
    _get_tso_resources

    {
        echo "Starting analysis for ${sample}"
        SECONDS=0

        sudo bash TSO500_ruo/TSO500_RUO_LocalApp/TruSight_Oncology_500_RUO.sh \
            --analysisFolder /home/dnanexus/out/Analysis/"${sample}"_output \
            --fastqFolder /home/dnanexus/DemultiplexOutput/ \
            --resourcesFolder /home/dnanexus/TSO500_ruo/TSO500_RUO_LocalApp/resources \
            --sampleSheet /home/dnanexus/SampleSheet.csv \
            --sampleOrPairIDs "$sample" \
            ${options} 2>&1 | tee /home/dnanexus/${sample}_scatter_stdout.txt
    
        echo Analysis complete
    } || {
        # some form of error occured in running that raised non-zero exit code
        code=$?
        echo "ERROR: one or more errors occured running workflow"
        echo "Process exited with code: ${code}"
        exit 1
    }

    # final check that workflow did successfully complete as it can continue
    # after failing workflow steps without raising non-zero exit code
    success=$(grep "SingleWorkflowRunnerActor workflow finished with status 'Succeeded'" \
        /home/dnanexus/${sample}_scatter_stdout.txt)

    if [[ -z "$success" ]]; then
        echo "All workflow steps did not successfully complete, exiting now"
        exit 1
    fi

    duration="$SECONDS"
    echo "Scatter job complete for ${sample} in $(($duration / 60))m$(($duration % 60))s"

    _upload_scatter_output
}


_gather() {
    : '''
    Run the final step to gather all per sample files and generate output
    '''
    sample_output_dirs=""
    for sample in $sample_list; do
        sample_output_dirs+="/home/dnanexus/out/Analysis/${sample}_output "
    done

    if [[ -f runfolder/modified_SampleSheet.csv ]]; then
        # include/exclude previously specified for scatter jobs =>
        # use modified samplesheet to prevent gather step raising an error
        samplesheet_path="runfolder/modified_SampleSheet.csv"
    else
        samplesheet_path="runfolder/SampleSheet.csv"
    fi

    SECONDS=0
    echo "Starting gather step to aggregate per sample results"

    /usr/bin/time -v sudo bash TSO500_ruo/TSO500_RUO_LocalApp/TruSight_Oncology_500_RUO.sh \
        --analysisFolder /home/dnanexus/out/Results \
        --runFolder /home/dnanexus/runfolder \
        --resourcesFolder /home/dnanexus/TSO500_ruo/TSO500_RUO_LocalApp/resources \
        --sampleSheet /home/dnanexus/$samplesheet_path \
        --gather /home/dnanexus/out/DemultiplexOutput/ "${sample_output_dirs}" \
        --isNovaSeq 2>&1 | tee /home/dnanexus/gather_run.log

    find out/Results -type f

    # final check it was successful
    if [[ $(grep "WorkflowFailedState" gather_run.log) ]]; then
        cat gather_run.log
        dx-jobutil-report-error "Gather step did not complete successfully"
    fi

    mv gather_run.log /home/dnanexus/out/Results/gather_run.log

    duration="$SECONDS"
    echo "Gather step completed in ${sample} in $(($duration / 60))m$(($duration % 60))s"
}


main() {
    : '''
    Main entry point for the app, general outline of behaviour:
    
    - download and unpack run data and TSO500 resources zip
    - run demultiplexing and upload output
    - call _scatter function to run sub job of analysis per sample
    - download all outputs from _scatter jobs
    - call _gather function locally to generate final results
    - upload output from _gather
    '''
    # control how many operations to open in parallel for download / upload
    THREADS=$(nproc --all)

    mkdir -p /home/dnanexus/runfolder \
             /home/dnanexus/TSO500_ruo \
             /home/dnanexus/out/DemultiplexOutput \
             /home/dnanexus/out/Logs/ \
             /home/dnanexus/out/Analysis \
             /home/dnanexus/out/Results

    # download sample data (either tar files or fastqs)
    _get_input_files

    # download, unpack and load TSO500 docker image & resources
    _get_tso_resources

    # download samplesheet if provided or get from run data and parse out sample names
    _get_samplesheet && _parse_samplesheet

    if [[ $n_samples ]] || [[ $include_samples ]] || [[ $exclude_samples ]]; then
        # specified a subset of samples to run analysis for =>
        # modify the samplesheet accordingly
        _modify_samplesheet
    fi

    echo "Starting analysis"

    _demultiplex
    _upload_demultiplex_output

    if [[ "$demultiplex_only" == false ]]; then
        # start up analysis, running one instance of the local app per sample in parallel

        # Adds additional non-specified optional arguments to the command
        options=""
        if [ "$analysis_options" ]; then options+=" ${analysis_options}"; fi 
        if [ "$isNovaSeq" = true ]; then options+=" --isNovaSeq"; fi

        # format inputs for launching jobs
        samplesheet_id=$(grep -oe "file-[0-9A-Za-z]*" <<< "$samplesheet")
        tso_ruo_id=$(grep -oe "file-[0-9A-Za-z]*" <<< "$TSO500_ruo")

        # start up job for every sample in scatter mode to run analysis
        for sample in $sample_list; do
            dx-jobutil-new-job _scatter \
                -isample="$sample" \
                -isamplesheet="$samplesheet_id" \
                -iTSO500_ruo="$tso_ruo_id" \
                -ifastqs="$fastq_ids" \
                -ioptions="$options" \
                --instance-type="$scatter_instance" \
                --extra-args='{"priority": "high"}' \
                --name "_scatter [${sample}]" >> job_ids
        done

        # add to logs total size of fastqs per sample to give indication
        # how long each sub job will take relative to each sample
        _calculate_total_fq_size "$fastq_ids"

        # wait until all scatter jobs complete
        dx wait --from-file job_ids

        # download all scatter output
        _get_scatter_job_outputs

        # gather all per sample output and do final processing for output results
        _gather

        echo "All workflows complete"

        # upload output files from gather step
        _upload_gather_output
    fi

    # exit 1

    # check usage to monitor usage of instance storage
    echo "Total file system usage"
    df -h
}
