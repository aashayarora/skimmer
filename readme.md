# NanoAOD Skimmer# NanoAOD Skimmer

This project contains scripts to run a NanoAOD skimmer using HTCondor. The skimmer selects events based on specific criteria and produces output files with the selected events.

## Requirements

- CMSSW environment (tested with CMSSW_14_1_4)
- Python 3
- ROOT
- HTCondor
- XRootD client

## Setup

1. Ensure you have a CMSSW environment set up. If not, set it up using the following commands:

    ```bash
    export SCRAM_ARCH=slc7_amd64_gcc820
    cmsrel CMSSW_14_1_4
    cd CMSSW_14_1_4/src
    cmsenv
    ```

## Running the Skimmer

To run the skimmer, use the `runSkimmer.sh` script. This script will create a tarball of the skimming code, generate a list of files from the specified datasets, and submit jobs to HTCondor.

### Usage

```bash
./runSkimmer <name-of-dataset>
```

`<input_datasets>: Comma-separated list of input datasets.`