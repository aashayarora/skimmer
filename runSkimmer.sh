#!/usr/bin/env bash
set -euo pipefail

# Configuration
CPUS=1
MEMORY=2G

X509_USER_PROXY="${X509_USER_PROXY:-/tmp/x509up_u$(id -u)}"

# Display usage information
usage() {
    cat <<EOF
Usage: $0 -i <input_list> [-s]
    -i  file with one dataset path per line (required)
    -s  submit signal
    -h  show this help message
    -v  version tag for the skimmer (required)
EOF
    exit 1
}

# Parse command line options
parse_options() {
    INPUT_LIST=""
    IS_SIG=0

    while getopts ":i:shv:" opt; do
        case "$opt" in
            i) INPUT_LIST=$OPTARG ;;
            s) IS_SIG=1 ;;
            v) VERSION=$OPTARG ;;
            h) usage ;;
            *) usage ;;
        esac
    done

    # Validate args
    if [[ -z "${INPUT_LIST}" ]]; then
        echo "Error: -i is required."
        usage
    fi

    if [[ -z "${VERSION}" ]]; then
        echo "Error: -v is required."
        usage
    fi

    # Set OUTPUT_TAG and LOGDIR after VERSION is validated
    OUTPUT_TAG="${VERSION}"
    LOGDIR="logs/${OUTPUT_TAG}"
}

# Get list of files for a dataset
get_file_list() {
    local ds=$1 
    local is_sig=$2 
    local out=$3

    echo "Getting file list for dataset: $ds"
    if (( is_sig )); then
        # Handle multiple directories separated by commas
        IFS=',' read -ra DIRS <<< "$ds"
        > "$out"  # Clear output file
        for dir in "${DIRS[@]}"; do
            dir=$(echo "$dir" | xargs)  # Trim whitespace
            if [[ -d "$dir" ]]; then
                ls -1 "${dir}"/*.root >> "$out" 2>/dev/null || 
                    { echo "Warning: No ROOT files found in $dir"; }
            else
                echo "Warning: Directory $dir does not exist"
            fi
        done
        
        if [[ ! -s "$out" ]]; then
            echo "Error: No ROOT files found in any of the specified directories"
            return 1
        fi
    else
        dasgoclient --query="file dataset=${ds}" >"$out" || 
            dasgoclient --query="file dataset=${ds}" >"$out" || 
            { echo "Error: Failed to query DAS for dataset $ds"; return 1; }
    fi
    
    local count=$(wc -l < "$out")
    echo "Found $count files"
}

# Create and submit a condor job
submit_job() {
    local file_list=$1 
    local sigflag=$2
    local proxy_path=$3
    
    local job_name=$(basename "$file_list")
    local jdl=(condor_${job_name}.jdl) || { echo "Error: Failed to create temp file"; return 1; }    

    echo "Creating job description file: $jdl"
    cat >"$jdl" <<EOF
universe                = Vanilla
request_cpus            = ${CPUS}
request_memory          = ${MEMORY}
executable              = executable.py
transfer_executable     = True
transfer_input_files    = keep_and_drop_skim.txt, truthSelections.h, jetId.h
arguments               = ${proxy_path} \$(FILE) ${sigflag} ${OUTPUT_TAG}
log                     = ${LOGDIR}/\$(Cluster).log 
output                  = ${LOGDIR}/\$(Cluster).out
error                   = ${LOGDIR}/\$(Cluster).err
+JobFlavour             = "tomorrow"

queue FILE from ${file_list}
EOF
    
    echo "Submitting job to condor"
    condor_submit "$jdl"
    rm -f "$jdl"  # Clean up the JDL file after submission
}

# Process each dataset and submit jobs
process_datasets() {
    local input_list=$1 
    local is_sig=$2
    local proxy_path=$3

    mkdir -p "$LOGDIR"
    echo "Processing datasets from: $input_list"
    
    while read -r ds; do
        [[ -z "$ds" || "$ds" =~ ^# ]] && continue  # Skip empty lines and comments
        
        echo "Processing dataset: $ds"
        local tmpfl="tmp_filelist_$(uuidgen)"
        
        if get_file_list "$ds" "$is_sig" "$tmpfl"; then
            submit_job "$tmpfl" "$is_sig" "$proxy_path"
        else
            echo "Warning: Skipping dataset $ds due to errors"
        fi

        rm -f "$tmpfl"
    done <"$input_list"

    echo "All jobs submitted."
}

# Main function
main() {
    # check if x509 proxy exists, and then copy it to afs
    proxy_path="/afs/cern.ch/user/$(whoami | head -c 1)/$(whoami)/x509up_u$(id -u)"
    if [[ -f "$X509_USER_PROXY" ]]; then
        cp "$X509_USER_PROXY" $proxy_path
    else
        echo "Error: X509 user proxy not found at $X509_USER_PROXY"
        exit 1
    fi

    local options=$(parse_options "$@")
    read -r input_list is_sig <<< "$options"
    parse_options "$@"
    
    process_datasets "$INPUT_LIST" "$IS_SIG" "$proxy_path"
}
# Execute main with all arguments
main "$@"
