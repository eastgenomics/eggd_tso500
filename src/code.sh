#!/bin/bash

# Output each line as it is executed (-x) and stop if any non zero exit codes are seen (-e)
set -exo pipefail

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

options=""

# Analysis starting from FastQ files
if [ "$isFastQ" = true ];
then
    echo "This analysis is starting from FastQ inputs"
    options+=" --fastqFolder /home/dnanexus/runfolder"

    # Download all FastQ files from the array
    echo "Downloading FastQ inputs"
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
            exit 1
        fi
        done

    cd /home/dnanexus

else
    # If the analysis input is not FastQs expect tar archive(s) of bcl files
    options+=" --runFolder /home/dnanexus/runfolder/" 

    # download the runfolder input, decompress and save in directory 'runfolder'
    echo " Downloading tar archive(s)"
    for i in "${!input_files_name[@]}" 
    do
        # check format and decompress
        name="${input_files_name[${i}]}"

        if [[ "${name}" == *.tar.gz ]] || [[ "${name}" == *.tgz ]]; then
            dx cat "${input_files[${i}]}" | tar zxf - --no-same-owner -C /home/dnanexus/runfolder

        else
            dx-jobutil-report-error "ERROR: The input was not a .tar.gz or .tgz"
            exit 1

        fi
    done

    if [ "$demultiplexOnly" = true ]; then options+=" --demultiplexOnly"; fi  

fi


#Find the SampleSheet if provided as an input or present in the input directory
if [ "$samplesheet" ];
then 
    echo "SampleSheet was provided as input ${samplesheet_name}"
    mv -v /home/dnanexus/in/samplesheet/*.csv /home/dnanexus/runfolder/SampleSheet.csv

elif [ -f /home/dnanexus/runfolder/*heet.csv ];
then
    echo "SampleSheet is present in the run directory"
    mv -v /home/dnanexus/runfolder/*heet.csv /home/dnanexus/runfolder/SampleSheet.csv 
else
    dx-jobutil-report-error "SampleSheet not provided or found, please check your inputs."
    exit 1
fi

# Adds additional non-specified optional arguments to the command
if [ "$analysis_options" ]; then options+=" ${analysis_options}"; fi 
if [ "$isNovaSeq" = true ]; then options+=" --isNovaSeq"; fi

# run the shell script, specifying the analysis folder, runfolder, samplesheet, resourcesFolder and any analysis options string given as an input
# pipe stderr into stdout and write this to a file and to screen - this allows a record of the logs to be saved and visible on screen if it goes wrong

echo "Starting analysis"

sudo bash TSO500_ruo/TSO500_RUO_LocalApp/TruSight_Oncology_500_RUO.sh \
--analysisFolder /home/dnanexus/out/analysis_folder/analysis_folder \
--resourcesFolder /home/dnanexus/TSO500_ruo/TSO500_RUO_LocalApp/resources \
${options} 2>&1 | tee /home/dnanexus/out/logs/logs/RUO_stdout.txt

if less /home/dnanexus/out/logs/logs/RUO_stdout.txt | grep WorkflowSucceededState ; 
then 
    echo "LocalApp Completed Successfully";
else 
    echo "Local App Failed"
    echo "SampleSheet Validation log:"
    less /home/dnanexus/out/analysis_folder/analysis_folder/Logs_Intermediates/SamplesheetValidation/SamplesheetValidation-*.log

    dx-jobutil-report-error "ERROR: Workflow failed to complete"
    exit 1
fi

# upload all outputs
dx-upload-all-outputs --parallel