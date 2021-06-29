#!/bin/bash

# Output each line as it is executed (-x) and don't stop if any non zero exit codes are seen (+e)
set -x +e
mark-section "download inputs"

mkdir -p runfolder TSO500_ruo out/logs/logs out/analysis_folder

# download all inputs
dx-download-all-inputs --parallel --except run_folder

unzip $TSO500_ruo_path -d TSO500_ruo
rm $TSO500_ruo_path
# change the owner of the app
chmod 777 TSO500_ruo/TSO500_RUO_LocalApp/
# move the docker image into ~
mv TSO500_ruo/TSO500_RUO_LocalApp/trusight-oncology-500-ruo-dockerimage-ruo-*.tar /home/dnanexus/
# load docker image
sudo docker load --input /home/dnanexus/trusight-oncology-500-ruo-dockerimage-ruo-*.tar

if $isFastQ
then 
    for i in ${!run_folder[@]}
        do 
            dx download "${run_folder[$i]}" -o runfolder/
        done
    sudo bash TSO500_ruo/TSO500_RUO_LocalApp/TruSight_Oncology_500_RUO.sh --analysisFolder /home/dnanexus/out/analysis_folder/analysis_folder --fastqFolder /home/dnanexus/runfolder --sampleSheet $samplesheet_path --resourcesFolder /home/dnanexus/TSO500_ruo/TSO500_RUO_LocalApp/resources $analysis_options 2>&1 | tee /home/dnanexus/out/logs/logs/RUO_stdout.txt
else
    # download the runfolder input, decompress and save in directory 'runfolder'
    dx cat "$run_folder" | tar zxf - -C runfolder
    # run the shell script, specifying the analysis folder, runfolder, samplesheet, resourcesFolder and any analysis options string given as an input
    # pipe stderr into stdout and write this to a file and to screen - this allows a record of the logs to be saved and visible on screen if it goes wrong
    sudo bash TSO500_ruo/TSO500_RUO_LocalApp/TruSight_Oncology_500_RUO.sh --analysisFolder /home/dnanexus/out/analysis_folder/analysis_folder --runFolder /home/dnanexus/runfolder/* --sampleSheet $samplesheet_path --resourcesFolder /home/dnanexus/TSO500_ruo/TSO500_RUO_LocalApp/resources $analysis_options 2>&1 | tee /home/dnanexus/out/logs/logs/RUO_stdout.txt
fi

# upload all outputs
dx-upload-all-outputs --parallel

