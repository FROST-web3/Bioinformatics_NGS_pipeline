#!/bin/bash

change_name(){
    local input_dir=$1
    local initial_dir=$(pwd)
    
    if [ ! -d "$input_dir" ]; then
        echo "Error: Directory $input_dir does not exist\n\n\n"
        return 1
    fi

    cd "$input_dir"
    current_dir=${PWD##*/}
    
    echo "\n\n\nProcessing directory: $current_dir "

    if [ -f *.narrowPeak ]; then
        mv *.narrowPeak "$current_dir.narrowPeak"
        echo "Renamed: .narrowPeak -> $current_dir.narrowPeak"
    else
        echo "Warning: No *.narrowPeak file found"
    fi

    if [ -f *.broadPeak ]; then
        mv *.broadPeak "$current_dir.broadPeak"
        echo "Renamed: .broadPeak -> $current_dir.broadPeak"
    else
        echo "Warning: No *.broadPeak file found"
    fi

    if [ -f *.out_treat_pileup.bdg.bw ]; then
        mv *.out_treat_pileup.bdg.bw "$current_dir.bw"
        echo "Renamed: .out_treat_pileup.bdg.bw -> $current_dir.bw"
    else
        echo "Warning: No *.out_treat_pileup.bdg.bw file found"
    fi

   if [ -f *.gappedPeak ]; then
        mv *.gappedPeak "$current_dir.gappedPeak"
        echo "Renamed: .gappedPeak -> $current_dir.gappedPeak"
    else
        echo "Warning: No *.gappedPeak file found"
    fi

   if [ -f *.xls ]; then
        mv *.xls "$current_dir.xls"
        echo "Renamed: .xls -> $current_dir.xls"
    else
        echo "Warning: No *.xls file found"
    fi


    cd "$initial_dir"
}

main(){
    for folder in "$@"; do
        change_name "$folder"
    done
}

main "$@"