#!/bin/bash

# Output each line as it is executed (-x) and stop if any non zero exit codes are seen (-e)
set -euxo pipefail

mark-section "download inputs"

mkdir -p /home/dnanexus/runfolder TSO500_ruo out/logs/logs out/analysis_folder

# download all inputs
dx-download-all-inputs --parallel --except input_files 

#Unload and prep docker image
unzip $TSO500_ruo_path -d TSO500_ruo 
rm $TSO500_ruo_path
# change the owner of the app
chmod 777 TSO500_ruo/TSO500_RUO_LocalApp/
# move the docker image into ~
mv TSO500_ruo/TSO500_RUO_LocalApp/trusight-oncology-500-ruo-dockerimage-ruo-*.tar /home/dnanexus/
# load docker image
sudo docker load --input /home/dnanexus/trusight-oncology-500-ruo-dockerimage-ruo-*.tar

options=()

if [ "$isFastQ" = true ];

options+=( " --fastqFolder /home/dnanexus/runfolder" )

    if [ ! "$samplesheet" ]
    then 
        dx-jobutil-report-error "Please provide a SampleSheet for analysis. :)"
    fi

then
    # Download all fastq files from the array
    for i in "${!input_files[@]}"
        do 
            dx download "${input_files[$i]}" -o runfolder/
        done
    
    cd runfolder
    # The app requires the fastq files to be present in individual sample folders
    for fastq_file in *fastq.gz; 
        do 
        if [[ -e "$fastq_file" ]]
        then
            # FastQ files expected  in standard naming format i.e. SampleID_S6_L002_R2_001.fastq.gz
            # Substitution below transforms SampleID_S6_L002_R2_001.fastq.gz into SampleID
            samplename=$(echo $fastq_file  | sed 's/_S[0-9]*_L[0-9]*_R[0-9]*_001.fastq.gz//'); 
            sample_dir="/home/dnanexus/runfolder/${samplename}"
            mkdir -p "$sample_dir"
            mv "$fastq_file"  "$sample_dir"
        else
            dx-jobutil-report-error "No FastQ files found, please check your inputs. :)"
        fi
        done

    cd /home/dnanexus

else
    options+=( "--runFolder /home/dnanexus/runfolder/" )
    # download the runfolder input, decompress and save in directory 'runfolder'
    for i in "${!input_files_name[@]}" 
    do
        # check format and decompress
        name="${input_files_name[${i}]}"

        if [[ "${name}" == *.tar.gz ]] || [[ "${name}" == *.tgz ]]; then
            dx cat "${input_files[${i}]}" | tar zxf - --no-same-owner -C /home/dnanexus/runfolder

        else
            dx-jobutil-report-error "ERROR: The input was not a .tar.gz or .tgz"

        fi
    done

    if [ "$demultiplexOnly" = true ]; then options+=("--demultiplexOnly"); fi  

fi

if [ "$isNovaSeq" = true ]; then options+=("--isNovaSeq"); fi
if [ "$samplesheet" ];then options+=("--sampleSheet $samplesheet_path");fi

# Adds additional non-specified optional arguments to the command
options+=( "${analysis_options}" )

# run the shell script, specifying the analysis folder, runfolder, samplesheet, resourcesFolder and any analysis options string given as an input
# pipe stderr into stdout and write this to a file and to screen - this allows a record of the logs to be saved and visible on screen if it goes wrong
    
    sudo bash TSO500_ruo/TSO500_RUO_LocalApp/TruSight_Oncology_500_RUO.sh \
    --analysisFolder /home/dnanexus/out/analysis_folder/analysis_folder \
    --resourcesFolder /home/dnanexus/TSO500_ruo/TSO500_RUO_LocalApp/resources \
    "${options[@]}" 2>&1 | tee /home/dnanexus/out/logs/logs/RUO_stdout.txt


# upload all outputs
dx-upload-all-outputs --parallel