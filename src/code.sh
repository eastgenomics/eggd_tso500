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
options+=( "${analysis_options}" )


if [ "$isFastQ" = true]; 
then
    # Download all fastq files from the array
    for i in "${!input_files[@]}"
        do 
            dx download "${input_files[$i]}" -o runfolder/
        done

    # run the shell script, specifying the analysis folder, fastq folder, samplesheet, resourcesFolder and any analysis options string given as an input
    # pipe stderr into stdout and write this to a file and to screen - this allows a record of the logs to be saved and visible on screen if it goes wrong
    sudo bash TSO500_ruo/TSO500_RUO_LocalApp/TruSight_Oncology_500_RUO.sh \
        --analysisFolder /home/dnanexus/out/analysis_folder/analysis_folder \
        --fastqFolder /home/dnanexus/runfolder \
        --sampleSheet "$samplesheet_path" \
        --resourcesFolder /home/dnanexus/TSO500_ruo/TSO500_RUO_LocalApp/resources \
        "${options[@]}" 2>&1 | tee /home/dnanexus/out/logs/logs/RUO_stdout.txt
else

    # download the runfolder input, decompress and save in directory 'runfolder'
    #dx cat "$input_files" | tar zxf - -C runfolder

    echo "Step 1: download input tars and unpack them"
  

    for i in ${!input_files_name[@]}; do
        # check format and decompress
        name="${input_files_name[${i}]}"

        if [[ "${name}" == *.tar.gz ]] || [[ "${name}" == *.tgz ]]; then
        dx cat "${input_files[${i}]}" | tar zxf - --no-same-owner -C /home/dnanexus/runfolder

        elif [ "${name}" == *.zip ]; then
        dx download "${input_files[$i]}" -o "${name}"
        unzip "${name}"

        elif [ "${name}" == *.rar ]; then
        dx download "${input_files[${i}]}" -o "${name}"
        unrar x "${name}"

        else
        dx-jobutil-report-error "ERROR: The input was not a .rar, .zip, .tar.gz or .tgz"

        fi
    done

    if [ "$DemultiplexOnly" = true ];
    then
        options+=("--demultiplexOnly")
    fi
    # run the shell script, specifying the analysis folder, runfolder, samplesheet, resourcesFolder and any analysis options string given as an input
    # pipe stderr into stdout and write this to a file and to screen - this allows a record of the logs to be saved and visible on screen if it goes wrong
    
    sudo bash TSO500_ruo/TSO500_RUO_LocalApp/TruSight_Oncology_500_RUO.sh \
    --analysisFolder /home/dnanexus/out/analysis_folder/analysis_folder \
    --runFolder /home/dnanexus/runfolder/ \
    --sampleSheet "$samplesheet_path" \
    --resourcesFolder /home/dnanexus/TSO500_ruo/TSO500_RUO_LocalApp/resources \
    "${options[@]}" 2>&1 | tee /home/dnanexus/out/logs/logs/RUO_stdout.txt
    

fi

# upload all outputs
dx-upload-all-outputs --parallel
