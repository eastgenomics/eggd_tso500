#!/bin/bash

# mem3_ssd1_v2_x96

# Output each line as it is executed (-x) and stop if any non zero exit codes are seen (-e)
PS4='$(date)\011 '
export TZ=Europe/London
set -exo pipefail


_get_tso_resources() {
    : '''
    Download and unpack the TSO500 resources zip file
    '''
    echo "Downloading and unpacking ${TSO500_ruo_name}"
    SECONDS=0

    dx download "$TSO500_ruo"
    unzip -q $TSO500_ruo_name -d TSO500_ruo
    chmod -R 777 TSO500_ruo/TSO500_RUO_LocalApp/
    sudo docker load --input TSO500_ruo/TSO500_RUO_LocalApp/trusight-oncology-500-ruo-dockerimage-ruo-*.tar

    duration=$SECONDS
    sleep 5
    echo "Downloading and unpacking ${TSO500_ruo_name} took $(($duration / 60))m$(($duration % 60))s"
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
            dx download {} --no-progress -o /home/dnanexus/fastqFolder/
        
        # The app requires the fastq files to be present in individual sample folders
        for fastq_file in $(find /home/dnanexus/fastqFolder/ -name "*fastq.gz" -printf "%f\n"); do 
            # FastQ files expected  in standard naming format i.e. SampleID_S6_L002_R2_001.fastq.gz
            # Substitution below transforms SampleID_S6_L002_R2_001.fastq.gz into SampleID
            samplename=$(echo $fastq_file  | sed 's/_S[0-9]*_L[0-9]*_R[0-9]*_001.fastq.gz//')
            sample_dir="/home/dnanexus/fastqFolder/${samplename}/"
            mkdir -p "$sample_dir"
            mv "/home/dnanexus/fastqFolder/${fastq_file}"  "$sample_dir"
        done

        sleep 5
        printf "FastQs downloaded:\n $(find fastqFolder -type f | sed 's/^/\n\t/g')\n"
        sleep 5

    else
        # If the analysis input is not FASTQs expect tar archive(s) of bcl files§
        # download the runfolder input, decompress and save in directory 'runfolder'
        echo " Downloading tar files"
        
        # drop the $dnanexus_link from the file IDs
        file_ids=$(grep -Po  "file-[\d\w]+" <<< "${input_files[@]}")

        echo "$file_ids" | xargs -P ${THREADS} -n1 -I{} sh -c \
            "dx cat {} | tar xzf - --no-same-owner --absolute-names -C /home/dnanexus/runfolder"
       
    fi

    duration=$SECONDS
    sleep 5
    echo "Downloaded $(wc -w <<< ${file_ids}) files in $(($duration / 60))m$(($duration % 60))s"
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

  We are going to upload all output in the folder structure as
  output by the local app, but set the following types of output
  files to arrays of their own:

  - FASTQs
  - BAMs
  - gVCFs
  - CombinedVariantOutput tsvs
  - MetricsOutput.tsv
  
  This is to allow downstream app(s) to take in these as inputs
  since the local app outputs multiple of the same file type

  TODO

  Need to:
  - separate out bams, gvcfs, fastqs, metricsOutput, CVO

  want to upload in same path, but set the file IDs of these to be distinct outputs

  '''
  echo "Total files to upload: $(find /home/dnanexus/out/ -type f | wc -l)"
  SECONDS=0

  # upload all output files in parallel and add to output spec
  export -f _upload_single_file 

  # find each of the sets of files we want to split as distinct outputs,
  # upload them and move them out of the path to prevent uploading again
  fastqs=$(find "/home/dnanexus/out/PATH" -type f )
  xargs -P ${THREADS} -n1 -I{} bash -c "_upload_single_file {} fastqs" <<< "$fastqs"
  mv "$fastqs" /tmp

  bams=$(find "/home/dnanexus/out/PATH" -type f )
  xargs -P ${THREADS} -n1 -I{} bash -c "_upload_single_file {} bams" <<< "$bams"
  mv "$bams" /tmp

  gvcfs=$(find "/home/dnanexus/out/PATH" -type f )
  xargs -P ${THREADS} -n1 -I{} bash -c "_upload_single_file {} gvcfs" <<< "$gvcfs"
  mv "$gvcfs" /tmp

  cvo=$(find "/home/dnanexus/out/PATH" -type f )
  xargs -P ${THREADS} -n1 -I{} bash -c "_upload_single_file {} cvo" <<< "$cvo"
  mv "$cvo" /tmp

  metrics_file_id=$(dx upload -p out/Results/MetricsOutput.tsv --path "$remote_path" --brief)
  dx-jobutil-add-output "metrics" "$metrics_file_id" --file
  mv out/Results/MetricsOutput.tsv /tmp

  # upload rest of files
  find "/home/dnanexus/out/" -type f | xargs -P ${THREADS} -n1 -I{} bash -c "_upload_single_file {}"

  duration=$SECONDS
  echo "Uploading took $(($duration / 60))m$(($duration % 60))s."
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
  local remote_path=$(sed s'/\/home\/dnanexus\/out\//' <<< "$file")

  echo "file: ${file}"

  file_id=$(dx upload -p "$file" --path "$remote_path" --brief)
  dx-jobutil-add-output "$field" "$file_id" --array
}


_demultiplex() {
    : '''
    Run demultiplexing via bclconvert in the TSO500 local app
    '''
    echo "Beginning demultiplexing"
    SECONDS=0

    /usr/bin/time -v sudo bash TSO500_ruo/TSO500_RUO_LocalApp/TruSight_Oncology_500_RUO.sh \
        --analysisFolder /home/dnanexus/fastqFolder/ \
        --runFolder /home/dnanexus/runfolder/ \
        --resourcesFolder /home/dnanexus/TSO500_ruo/TSO500_RUO_LocalApp/resources \
        --sampleSheet /home/dnanexus/runfolder/SampleSheet.csv \
        --demultiexOnly
    
    duration=$SECONDS

    echo "Demultiplexing completed in $(($duration / 60))m$(($duration % 60))s"
}


_scatter() {
    : '''
    Run the initial per sample analysis from the TSO500 local app

    This runs the TSO500 local app in gather mode, running per sample
    in parallel with a maximum number of parallel jobs defined by the
    "n_proc" input (default: 8). This will need balancing against the
    total CPU and memory available, as if the number of parallel jobs
    compute requirements peak usage exceeds the total available in
    the instance it will cause the job to fail.
    '''
    SECONDS=0
    echo "Starting scatter jobs for $(wc -w <<< $sample_list) samples with ${n_proc} processes"

    # Adds additional non-specified optional arguments to the command
    OPTIONS=""
    if [ "$analysis_options" ]; then OPTIONS+=" ${analysis_options}"; fi 
    if [ "$isNovaSeq" = true ]; then OPTIONS+=" --isNovaSeq"; fi

    {
        echo "$sample_list" | /usr/bin/time -v xargs -n1 -P"$n_proc" -I{} sh -c " \
            echo Starting analysis for {} && \
            sudo bash TSO500_ruo/TSO500_RUO_LocalApp/TruSight_Oncology_500_RUO.sh \
            --analysisFolder /home/dnanexus/out/analysis/{}_output \
            --fastqFolder /home/dnanexus/fastqFolder \
            --resourcesFolder /home/dnanexus/TSO500_ruo/TSO500_RUO_LocalApp/resources \
            --sampleSheet /home/dnanexus/runfolder/SampleSheet.csv \
            --sampleOrPairIDs {} \
            ${OPTIONS} 1> /home/dnanexus/out/logs/logs/{}_stdout.txt &&
            echo Analysis complete for {}"
    } || {
        # some form of error occured in running, dump the end of each
        # files stdout/stderr to the logs for debugging and exit
        code=$?
        printf "ERROR: one or more errors occured during per sample scatter jobs"
        printf "Process exited with code: ${code}"
        printf "Recent logs of each sample analysis:\n \n"
        for sample in $(echo "$sample_list"); do
            printf "${sample} logs from ${sample}_stdout.txt"
            tail -n20 "${sample}_stdout.txt"
            printf "\n \n"s
            exit 1
        done
    }

    sleep 10  # add a wait to ensure logs all get written and aren't as messy

    duration="$SECONDS"
    echo "Scatter steps complete for all samples in $(($duration / 60))m$(($duration % 60))s"
}


_gather() {
    : '''
    Run the final step to gather all per sample files and generate output
    '''
        sample_output_dirs=""
        for sample in $sample_list; do
            sample_output_dirs+="/home/dnanexus/out/analysis/${sample}_output "
        done

        /usr/bin/time -v sudo bash TSO500_ruo/TSO500_RUO_LocalApp/TruSight_Oncology_500_RUO.sh \
            --analysisFolder /home/dnanexus/out/GatheredResults \
            --runFolder /home/dnanexus/runFolder \
            --resourcesFolder /home/dnanexus/TSO500_ruo/TSO500_RUO_LocalApp/resources \
            --sampleSheet /home/dnanexus/runfolder/SampleSheet.csv \
            --gather /home/dnanexus/fastqFolder/ "${sample_output_dirs}" \
            --isNovaSeq
}


main() {
    THREADS=$(nproc --all)  # control how many operations to open in parallel for download / upload

    mkdir -p /home/dnanexus/runfolder \
             /home/dnanexus/fastqFolder \
             /home/dnanexus/TSO500_ruo \
             /home/dnanexus/out/logs/logs \
             /home/dnanexus/out/analysis

    # download samplesheet if provided or get from run data and parse out sample names
    _get_samplesheet

    _parse_sample_names

    # download sample data (either tar files or fastqs)
    _get_input_files

    # download, unpack and load TSO500 docker image & resources
    _get_tso_resources


    echo "Starting analysis"

    if [[ "$isFastQ" == false ]]; then
        # not starting from fastqs => run demultiplexing
        _demultiplex
    fi

    if [[ "$demultiplexOnly" == false ]]; then
        # start up analysis, running one instance of the local app per sample in parallel
        _scatter

        # gather all per sample output and do final processing for output results
        _gather

        # final check of if everything was successful
        if less /home/dnanexus/out/logs/logs/RUO_stdout.txt | grep WorkflowSucceededState ; 
        then 
            echo "LocalApp Completed Successfully"
        else 
            echo "Local App Failed"
            echo "SampleSheet Validation log:"
            less /home/dnanexus/out/analysis/analysis/Logs_Intermediates/SamplesheetValidation/SamplesheetValidation-*.log

            dx-jobutil-report-error "ERROR: Workflow failed to complete"
            exit 1
        fi
    fi

    # upload our output files
    _upload_all_output

    # check usage to monitor usage of instance storage
    echo "Total file system usage"
    df -h

}