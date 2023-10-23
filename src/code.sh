#!/bin/bash

PS4='\000[$(date)]\011'
export TZ=Europe/London
set -exo pipefail


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

    set +x

    # files from sub jobs will be in the container- project context of the
    # current job ($DX_WORKSPAce-id) => search here for  all the files
    scatter_files=$(dx find data --json --verbose --path "$DX_WORKSPACE_ID:/analysis")

    # turn describe output into id:/path/to/file to download with dir structure
    files=$(jq -r '.[] | .id + ":" + .describe.folder + "/" + .describe.name'  <<< $scatter_files)

    # build aggregated directory structure and download all files
    cmds=$(for f in  $files; do \
        id=${f%:*}; path=${f##*:}; dir=$(dirname "$path"); \
        echo "'mkdir -p out/$dir && dx download --no-progress $id -o out/$path'"; done)

    echo $cmds | xargs -n1 -P${THREADS} bash -c

    set -x

    total=$(du -sh /home/dnanexus/out/analysis/ | cut -f1)
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
        samplesheet=$(find ./ -regextype posix-extended  -iregex '.*sample[-_ ]?sheet.csv$')
        echo "Using sample sheet in run directory: $samplesheet"
        mv "$samplesheet" /home/dnanexus/runfolder/SampleSheet.csv
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
    
    This will only be called after demultiplexing to not result in large
    undertermined fastqs.
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

    if [[ "$include_samples" ]]; then
        # retaining rows containing only those specified for given samples
        echo -e "-iinclude_samples specified:\n\t${include_samples}"
        include=$(sed 's/,/|/g' <<< "$include_samples")
        sample_rows=$(awk '/'"$include"'/ {print $1}' <<< "{$sample_rows}")
    fi

    if [[ "$exclude_samples" ]]; then
        # exclude rows containing only those specified for given samples
        echo -e "-iexclude_samples specified:\n\t${exclude_samples}"
        exclude=$(sed 's/,/|/g' <<< "$exclude_samples")
        sample_rows=$(awk -v exclude="$exclude" '$1 ~ /^include/' <<< "$exclude")
        sample_rows=$(awk '!/'"$exclude"'/ {print $1}' <<< "{$sample_rows}")
    fi

    # write out new samplesheet with specified rows
    mv runfolder/SampleSheet.csv runfolder/originalSampleSheet.csv
    echo "$samplesheet_header" >> runfolder/SampleSheet.csv
    echo "$sample_rows" >> runfolder/SampleSheet.csv

    sleep 1 && echo -e "New samplesheet:\n$(cat runfolder/SampleSheet.csv)" && sleep 1
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
    set +x
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


_upload_scatter_output() {
    : '''
    Upload the per sample output files in a scatter job

    To speed up the upload, logs and stdout/stderr files and  all cromwell
    logs in cromwell-executions/ are also combined into a tar files.

    To generate distinct output fields for the app for result files, the bam
    and index, gVCF and CombinedVariantOutput are first uploaded to individual
    output fields, then the rest of files are uploaded in parallel.
    '''
    SECONDS=0
    echo "Uploading sample output"
    export -f _upload_single_file

    # tar up all cromwell logs for faster upload
    tar -I pigz -cf /home/dnanexus/out/analysis/${sample}_output/${sample}_cromwell_executions.tar.gz \
        /home/dnanexus/out/analysis/${sample}_output/cromwell-executions
    rm -rf /home/dnanexus/out/analysis/${sample}_output/cromwell-executions/

    # tar up all the logs, stdout and stderr for faster upload
    logs=$(find out/analysis/${sample}_output/ \
        -name "*.log" -o \
        -name "*.out" -o \
        -name "*.stderr" -o \
        -name "*.stdout" -o\
        -name "receipt")
    mkdir all_logs && xargs -I{} mv {} all_logs/ <<< "$logs"
    tar -czf ${sample}_scatter_logs.tar.gz all_logs
    mv ${sample}_scatter_logs.tar.gz /home/dnanexus/out/analysis/${sample}_output/

    mv /home/dnanexus/${sample}_scatter_stdout.txt /home/dnanexus/out/analysis/${sample}_output/

    # move selected files to be distinct outputs, including bam, index, gVCF and CombinedVariantOutput
    bams=$(find "/home/dnanexus/out/analysis/${sample}_output/Logs_Intermediates/StitchedRealigned/" \
        -type f -name "*.bam")
    xargs -P ${THREADS} -n1 -I{} bash -c "_upload_single_file {} bams" <<< "$bams"
    xargs -n1 -I{} mv {} /tmp <<< $bams

    index=$(find "/home/dnanexus/out/analysis/${sample}_output/Logs_Intermediates/StitchedRealigned/" \
        -type f -name "*.bai")
    xargs -P ${THREADS} -n1 -I{} bash -c "_upload_single_file {} bam_index" <<< "$index"
    xargs -n1 -I{} mv {} /tmp <<< $index

    gvcf=$(find "/home/dnanexus/out/analysis/${sample}_output/Results/${sample}/" \
        -type f -name "*.genome.vcf")
    xargs -P ${THREADS} -n1 -I{} bash -c "_upload_single_file {} gvcfs" <<< "$gvcf"
    xargs -n1 -I{} mv {} /tmp <<< $gvcf

    cvo=$(find "/home/dnanexus/out/analysis/${sample}_output/Results/${sample}/" \
        -type f -name "*CombinedVariantOutput.tsv")
    xargs -P ${THREADS} -n1 -I{} bash -c "_upload_single_file {} cvo" <<< "$cvo"
    xargs -n1 -I{} mv {} /tmp <<< $cvo

    # upload rest of files
    find /home/dnanexus/out/ -type f | xargs -P ${THREADS} -n1 -I{} bash -c \
        "_upload_single_file {} analysis_folder"
    
    duration="$SECONDS"
    echo "Uploading completed in ${sample} in $(($duration / 60))m$(($duration % 60))s"
}


_upload_gather_output() {
    : '''
    Upload the output files from the final _gather step,
    these will all be in /home/dnanexus/out/Results
    '''
    echo "Total files to upload: $(find /home/dnanexus/out/Results -type f | wc -l)"
    SECONDS=0

    # link output of sub scatter jobs to parent job
    while read -r job; do
        dx-jobutil-add-output bams "${job}":bams --class=array:jobref
        dx-jobutil-add-output bam_index "${job}":bam_index --class=array:jobref
        dx-jobutil-add-output gvcfs "${job}":gvcfs --class=array:jobref
        dx-jobutil-add-output cvo "${job}":cvo --class=array:jobref
        dx-jobutil-add-output analysis_folder "${job}":analysis_folder --class=array:jobref

    done < job_ids

    export -f _upload_single_file

    # tar up all cromwell logs for faster upload
    tar -I pigz -cf /home/dnanexus/out/Results/gather_cromwell_executions.tar.gz \
        /home/dnanexus/out/Results/cromwell-executions
    rm -rf /home/dnanexus/out/Results/cromwell-executions

    # tar up all the logs, stdout and stderr for faster upload
    logs=$(find out/Results/ \
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
    tar -czf gather_logs.tar.gz all_logs
    mv gather_logs.tar.gz /home/dnanexus/out/Results/

    # upload final MetricsOutput.tsv as distinct output field
    metrics_file_id=$(dx upload -p /home/dnanexus/out/Results/Results/MetricsOutput.tsv --brief)
    dx-jobutil-add-output "metricsOutput" "$metrics_file_id" --class=file
    mv /home/dnanexus/out/Results/Results/MetricsOutput.tsv /tmp

    # upload rest of files
    find "/home/dnanexus/out/Results" -type f | xargs -P ${THREADS} -n1 -I{} bash -c \
        "_upload_single_file {} analysis_folder"

    duration=$SECONDS
    echo "Uploading took $(($duration / 60))m$(($duration % 60))s."
}


_upload_single_file(){
  : '''
  Uploads single file with dx upload and associates uploaded
  file ID to specified output field

  Arguments
  ---------
    1 : str
        path and file to upload
    2 : str
        app output field to link the uploaded file to
  '''
  local file=$1
  local field=$2
  local remote_path=$(sed s'/\/home\/dnanexus\/out\///' <<< "$file")

  file_id=$(dx upload "$file" --path "$remote_path" --parents --brief)
  dx-jobutil-add-output "$field" "$file_id" --array
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
    fastqs=$(find "/home/dnanexus/out/fastqFolder" -type f -name "*fastq.gz")
    xargs -P ${THREADS} -n1 -I{} bash -c "_upload_single_file {} fastqs" <<< "$fastqs"
    xargs -n1 -I{} mv {} /tmp <<< $fastqs

    # tar up all cromwell logs for faster upload
    tar -I pigz -cf /home/dnanexus/out/fastqFolder/demultiplex_cromwell_executions.tar.gz \
        /home/dnanexus/out/fastqFolder/cromwell-executions
    rm -rf /home/dnanexus/out/fastqFolder/cromwell-executions

    # upload rest of files
    find "/home/dnanexus/out/fastqFolder/" -type f | xargs -P ${THREADS} -n1 -I{} bash -c \
        "_upload_single_file {} analysis_folder"

    duration=$SECONDS
    echo "Uploading took $(($duration / 60))m$(($duration % 60))s."
}


_demultiplex() {
    : '''
    Run demultiplexing via bclconvert in the TSO500 local app
    '''
    echo "Beginning demultiplexing"
    SECONDS=0

    /usr/bin/time -v sudo bash TSO500_ruo/TSO500_RUO_LocalApp/TruSight_Oncology_500_RUO.sh \
        --analysisFolder /home/dnanexus/out/fastqFolder/ \
        --runFolder /home/dnanexus/runfolder/ \
        --resourcesFolder /home/dnanexus/TSO500_ruo/TSO500_RUO_LocalApp/resources \
        --sampleSheet /home/dnanexus/runfolder/SampleSheet.csv \
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
    PS4='\000[$(date)]\011'
    export TZ=Europe/London
    set -exo pipefail

    # control how many operations to open in parallel for download / upload
    THREADS=$(nproc --all)

    mkdir -p /home/dnanexus/runfolder \
             /home/dnanexus/TSO500_ruo \
             /home/dnanexus/fastqFolder/$sample \
             /home/dnanexus/out/logs/logs \
             /home/dnanexus/out/analysis

    dx download "$samplesheet" -o "SampleSheet.csv"

    # download the sample fastqs into directory for local app
    SECONDS=0
    echo "Finding fastqs for sample ${sample} from fastqs provided to job"
    set +x  # suppress this going to the logs as its v long
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
        "dx download --no-progress -o /home/dnanexus/fastqFolder/$sample/ {}"
    
    duration=$SECONDS
    echo "Downloaded fastqs in $(($duration / 60))m$(($duration % 60))s"
    
    # download and unpack local app resources
    _get_tso_resources

    {
        echo "Starting analysis for ${sample}"
        SECONDS=0

        sudo bash TSO500_ruo/TSO500_RUO_LocalApp/TruSight_Oncology_500_RUO.sh \
            --analysisFolder /home/dnanexus/out/analysis/"${sample}"_output \
            --fastqFolder /home/dnanexus/fastqFolder/ \
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
        sample_output_dirs+="/home/dnanexus/out/analysis/${sample}_output "
    done

    SECONDS=0
    echo "Starting gather step to aggregate per sample results"

    /usr/bin/time -v sudo bash TSO500_ruo/TSO500_RUO_LocalApp/TruSight_Oncology_500_RUO.sh \
        --analysisFolder /home/dnanexus/out/Results \
        --runFolder /home/dnanexus/runfolder \
        --resourcesFolder /home/dnanexus/TSO500_ruo/TSO500_RUO_LocalApp/resources \
        --sampleSheet /home/dnanexus/runfolder/SampleSheet.csv \
        --gather /home/dnanexus/out/fastqFolder/ "${sample_output_dirs}" \
        --isNovaSeq
    
    duration="$SECONDS"
    echo "Gather step completed in ${sample} in $(($duration / 60))m$(($duration % 60))s"

    # final check it was successful
    if [[ $(grep -i "transitioned to failed" out/Results/gather_run.log) ]]; then
        echo "Gather step did not complete successfully"
        cat out/Results/gather_run.log
        exit 1
    fi
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
             /home/dnanexus/out/fastqFolder \
             /home/dnanexus/out/logs/logs \
             /home/dnanexus/out/analysis

    # download sample data (either tar files or fastqs)
    _get_input_files

    # download, unpack and load TSO500 docker image & resources
    _get_tso_resources

    # download samplesheet if provided or get from run data and parse out sample names
    _get_samplesheet && _parse_samplesheet

    echo "Starting analysis"

    _demultiplex
    _upload_demultiplex_output

    if [[ $n_samples ]] || [[ $include_samples ]] || [[ $exclude_samples ]]; then
        # specified a subset of samples to run analysis for =>
        # modify the samplesheet accordingly
        _modify_samplesheet
    fi

    if [[ "$demultiplexOnly" == false ]]; then
        # start up analysis, running one instance of the local app per sample in parallel

        # get the file IDs of the fastqs we have just uploaded from demultiplexing
        fastq_ids=$(cat job_output.json | jq -r '.fastqs[][]')

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

        if [[ "$upload_demultiplex_output" == false ]]; then
            # option specified to not upload demultiplex output => delete it
            # from the parent container to not go into the output project
            echo "Removing demultiplexing output"
            dx rm -rf /fastqFolder
        fi
    fi

    # check usage to monitor usage of instance storage
    echo "Total file system usage"
    df -h
}
