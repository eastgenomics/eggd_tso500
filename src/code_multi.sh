#!/bin/bash

PS4='\000[$(date)]\011'
export TZ=Europe/London
set -exo pipefail


_get_tso_resources() {
    : '''
    Download and unpack the TSO500 resources zip file
    '''
    local zip_file=$(dx describe --json "$TSO500_ruo" | jq -r '.name')
    echo "Downloading and unpacking TSO500 resources zip file"

    dx download "$TSO500_ruo"
    unzip -q "$zip_file" -d TSO500_ruo
    chmod -R 777 TSO500_ruo/TSO500_RUO_LocalApp/
    sudo docker load --input TSO500_ruo/TSO500_RUO_LocalApp/trusight-oncology-500-ruo-dockerimage-ruo-*.tar

    duration=$SECONDS
    sleep 5
    echo "Downloading and unpacking ${zip_file} took $(($duration / 60))m$(($duration % 60))s"
}


_get_input_files() {
    : '''
    Download all input files (either fastQs or tar files from dx-streaming-upload)
    '''
    SECONDS=0
    echo "Downloading input files"

    if [ "$isFastQ" = true ]; then
        echo "This analysis is starting from FastQ inputs"

        # drop the $dnanexus_link from the file IDs
        file_ids=$(grep -Po  "file-[\d\w]+" <<< "${input_files[@]}")

        echo "$file_ids" | xargs -P ${THREADS} -n1 -I{} \
            dx download {} --no-progress -o /home/dnanexus/out/fastqFolder/
        
        total=$(du -sh /home/dnanexus/out/fastqFolder/ | cut -f1)
        
        # app requires the fastq files to be present in individual sample folders
        set +x  # disable set to limit dumping 4 extra lines per file into logs
        files=$(find /home/dnanexus/out/fastqFolder/ -name "*fastq.gz" -printf "%f\n")
        for fastq_file in $files; do 
            # FastQ files expected  in standard naming format i.e. SampleID_S6_L002_R2_001.fastq.gz
            # Substitution below transforms SampleID_S6_L002_R2_001.fastq.gz into SampleID
            samplename=$(echo $fastq_file | sed 's/_S[0-9]*_L[0-9]*_R[0-9]*_001.fastq.gz//')
            sample_dir="/home/dnanexus/out/fastqFolder/${samplename}/"
            mkdir -p "$sample_dir"
            mv "/home/dnanexus/out/fastqFolder/${fastq_file}"  "$sample_dir"
        done
        set -x
    else
        # If the analysis input is not FASTQs expect tar archive(s) of bcl files§
        # download the runfolder input, decompress and save in directory 'runfolder'
        echo " Downloading tar files"
        
        # drop the $dnanexus_link from the file IDs
        file_ids=$(grep -Po  "file-[\d\w]+" <<< "${input_files[@]}")

        echo "$file_ids" | xargs -P ${THREADS} -n1 -I{} sh -c \
            "dx cat {} | tar xzf - --no-same-owner --absolute-names -C /home/dnanexus/runfolder"

        total=$(du -sh /home/dnanexus/runfolder | cut -f1)
    fi

    duration=$SECONDS
    sleep 5
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

    output_path=$(dx describe --json "$DX_JOB_ID" | jq -r '.folder')
    scatter_files=$(dx find data --json --verbose --path "$output_path")

    # filter down files in job output folder to ensure they're from
    # one of our scatter jobs that just ran
    job_ids=$(sed 's/\n/|/g;s/|$//' job_ids)
    scatter_files=$(jq "map(select(.describe.createdBy.job | test(\"$job_ids\")))" <<< $scatter_files)

    # turn describe output into id:/path/to/file to download with dir structure
    files=$(jq -r '.[] | .id + ":" + .describe.folder + "/" + .describe.name'  <<< $scatter_files)

    # remove the beginning of the remote path (i.e. what was set for the app output)
    # to just leave the path we had in the scatter job for the output
    files=$(awk '{gsub("\/(.*?)analysis\/", ""); print}' <<< $files)

    # download all the files and build aggregated directory structure
    cmds=$(for f in  $files; do \
        id=${f%:*}; path=${f##*:}; dir=$(dirname $path); \
        echo "'mkdir -p out/analysis/$dir && dx download -r $id -o out/analysis/$path'"; done)

    echo $cmds | xargs -n1 -P${THREADS} bash -c

    # xargs -n1 -P${THREADS} -Ifile bash -c "IFS=: read -r id path <<< file; echo dx download $id -o out/analysis/$path" <<< $files

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


_parse_sample_names() {
    : '''
    Parse the sample names from the samplesheet, and optionally limit analysis to n
    number of samples if -in_samples specified
    '''
    # parse list of sample names from samplesheet
    sample_list=$(cat runfolder/SampleSheet.csv | awk '/Sample_ID/{ y=1; next }y'  | cut -d, -f1)
    echo "Samples parsed from samplesheet: ${sample_list}"

    if [[ "$n_samples" ]]; then
        echo "Limiting analysis to ${n_samples} samples"
        sample_list=$(sed '2p;d' <<< "${sample_list}")
        echo "Samples to run analysis for: ${sample_list}"
    fi
    sleep 5
}


_upload_all_output() {
    : '''
    Upload all data in /home/dnanexus/out/

    We are going to upload all outputs in the folder structure as
    output by the local app, but set the following types of output
    files to arrays of their own:

    - FASTQs
    - BAMs
    - gVCFs
    - CombinedVariantOutput tsvs
    - MetricsOutput.tsv
    
    This is to allow downstream app(s) to take in these as inputs
    since the local app outputs multiple of the same file type

    Additionally, all files matching *stdout*, *stderr* *log will be
    combined into a single tar for uploading
    '''
    echo "Total files to upload: $(find /home/dnanexus/out/ -type f | wc -l)"
    SECONDS=0

    # upload all output files in parallel and add to output spec
    export -f _upload_single_file

    set +x  # disable writing all commands to logs uploading since there are thousands

    # find each of the sets of files we want to split as distinct outputs,
    # upload them and move them out of the path to prevent uploading again
    # fastqs=$(find "/home/dnanexus/out/PATH" -type f -name "*fastq.gz")
    # xargs -P ${THREADS} -n1 -I{} bash -c "_upload_single_file {} fastqs" <<< "$fastqs"
    # mv "$fastqs" /tmp

    # bams=$(find "/home/dnanexus/out/analysis/Logs_Intermediates/StitchedRealigned/" -type f -name "*bam|*bai")
    # xargs -P ${THREADS} -n1 -I{} bash -c "_upload_single_file {} bams" <<< "$bams"
    # mv "$bams" /tmp

    # gvcfs=$(find "/home/dnanexus/out/analysis/Logs_Intermediates/VariantCaller/" -type f -name "*genome.vcf")
    # xargs -P ${THREADS} -n1 -I{} bash -c "_upload_single_file {} gvcfs" <<< "$gvcfs"
    # mv "$gvcfs" /tmp

    # cvo=$(find "/home/dnanexus/out/analysis/Logs_Intermediates/" -type f "*CombinedVariantOutput.tsv")
    # xargs -P ${THREADS} -n1 -I{} bash -c "_upload_single_file {} cvo" <<< "$cvo"
    # mv "$cvo" /tmp

    # find all the stdout / stderr log files and combine to
    # single tar to speed up uploading
    logs=$(find "/home/dnanexus/out/" -name "*stderr*|*stdout*|*.log")
    mkdir all_logs && mv "$logs" all_logs/
    tar cf all_logs.tar.gz "$logs"

    log_file_id=$(dx upload -p all_logs.tar.gz --brief)
    dx-jobutil-add-output "logs" "$log_file_id" --file

    metrics_file_id=$(dx upload -p out/analysis/Results/MetricsOutput.tsv --path "$remote_path" --brief)
    dx-jobutil-add-output "metrics" "$metrics_file_id" --file
    mv out/analysis/Results/MetricsOutput.tsv /tmp

    rm -rf /home/dnanexus/out/analysis/Logs_Intermediates/

    # upload rest of files
    find "/home/dnanexus/out/" -type f | xargs -P ${THREADS} -n1 -I{} bash -c "_upload_single_file {} analysis_folder"

    duration=$SECONDS
    echo "Uploading took $(($duration / 60))m$(($duration % 60))s."
    set -x
}


_upload_single_file(){
  : '''
  Uploads single file with dx upload and associates uploaded file ID to output spec

  Arguments
  ---------
    1 : str
        path and file to upload
    2 : str
        app output field to set the uploaded file to
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

    # upload rest of files
    find "/home/dnanexus/out/fastqFolder/" -type f | xargs -P ${THREADS} -n1 -I{} bash -c "_upload_single_file {} analysis_folder"

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
    THREADS=$(echo $(nproc --all) / 2 | bc)

    mkdir -p /home/dnanexus/runfolder \
             /home/dnanexus/TSO500_ruo \
             /home/dnanexus/fastqFolder/$sample \
             /home/dnanexus/out/logs/logs \
             /home/dnanexus/out/analysis

    dx download "$samplesheet" -o "SampleSheet.csv"

    # download the sample fastqs into directory for local app
    SECONDS=0
    set +x  # suppress this going to the logs as its v long
    details=$(xargs -n1 -P${THREADS} dx describe --json --verbose <<< $fastqs)
    sample_fqs=$(jq -r "select(.name | startswith(\"${sample}_\")) | .id" <<< $details)
    set -x

    echo "sample fastqs parsed: ${sample_fqs}"
    echo $sample_fqs | xargs -n1 -P${THREADS} -I{} sh -c "dx download --no-progress -o /home/dnanexus/fastqFolder/$sample/ {}"
    
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
            ${options} 1> /home/dnanexus/out/logs/logs/${sample}_stdout.txt
    
        echo Analysis complete for "$sample"
    } || {
        # some form of error occured in running, dump the end of each
        # files stdout/stderr to the logs for debugging and exit
        code=$?
        printf "ERROR: one or more errors occured during per sample scatter jobs"
        # printf "Process exited with code: ${code}"
        # printf "Recent logs of each sample analysis:\n \n"
        # for sample in $(echo "$sample_list"); do
        #     printf "${sample} logs from ${sample}_stdout.txt"
        #     tail -n20 "${sample}_stdout.txt"
        #     printf "\n \n"s
        #     exit 1
        # done
    }

    sleep 10  # add a wait to ensure logs all get written and aren't (quite) as messy

    duration="$SECONDS"
    echo "Scatter job complete for ${sample} in $(($duration / 60))m$(($duration % 60))s"

    SECONDS=0
    echo "Uploading sample output"

    # find all the stdout / stderr log files and combine to
    # single tar to speed up uploading
    # logs=$(find "/home/dnanexus/out/" -name "*stderr*|*stdout*|*.log")
    # mkdir all_logs && mv "$logs" all_logs/
    # tar cf ${sample}_all_logs.tar.gz "$logs"

    # log_file_id=$(dx upload -p all_logs.tar.gz --brief)
    # dx-jobutil-add-output "logs" "$log_file_id" --file

    # upload rest of files, first move all output files to have the
    # job output in the path for correct upload location
    job_output_path=$(dx describe --json "$DX_JOB_ID" | jq -r '.folder')
    mkdir /home/dnanexus/out/${job_output_path}
    mv /home/dnanexus/out/analysis/ /home/dnanexus/out/${job_output_path}/

    export -f _upload_single_file
    find /home/dnanexus/out/ -type f | xargs -P ${THREADS} -n1 -I{} bash -c "_upload_single_file {} analysis_folder"
    
    duration="$SECONDS"
    echo "Uploading completed in ${sample} in $(($duration / 60))m$(($duration % 60))s"

    find /home/dnanexus/out/ -type f
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
        --analysisFolder /home/dnanexus/out/GatheredResults \
        --runFolder /home/dnanexus/runfolder \
        --resourcesFolder /home/dnanexus/TSO500_ruo/TSO500_RUO_LocalApp/resources \
        --sampleSheet /home/dnanexus/runfolder/SampleSheet.csv \
        --gather /home/dnanexus/out/fastqFolder/ "${sample_output_dirs}" \
        --isNovaSeq
    
    duration="$SECONDS"
    echo "Gather step completed in ${sample} in $(($duration / 60))m$(($duration % 60))s"
}


main() {
    THREADS=$(nproc --all)  # control how many operations to open in parallel for download / upload

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
    _get_samplesheet && _parse_sample_names

    echo "Starting analysis"

    if [[ "$isFastQ" == false ]]; then
        # not starting from fastqs => run demultiplexing
        _demultiplex

        _upload_demultiplex_output
    fi

    if [[ "$demultiplexOnly" == false ]]; then
        # start up analysis, running one instance of the local app per sample in parallel

        if [ "$isFastQ" == true ]; then
            # get the file IDs of the fastqs provided as input
            fastq_ids=$(grep -Po  "file-[\d\w]+" <<< "${input_files[@]}")
        else
            # get the file IDs of the fastqs we have just uploaded from demultiplexing
            fastq_ids=$(cat job_output.json | jq -r '.fastqs[][]')
        fi
        
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

        # wait until all scatter jobs complete
        dx wait --from-file job_ids

        # download all scatter output
        _get_scatter_job_outputs

        # gather all per sample output and do final processing for output results
        _gather

        echo "All workflows complete"

        # upload our output files
        _upload_all_output
    fi

    # check usage to monitor usage of instance storage
    echo "Total file system usage"
    df -h

}