# NanoAOD Skimmer for VBS VVH Analysis

**Note: This is designed to run only on lxplus/CERN condor pool**

## Files

- `runSkimmer.sh` - What you need to run (creates condor submission scripts and submits them)
- `executable.py` - ROOT skimmer with VBS VVH selection + jet ID
- `keep_and_drop_skim.txt` - Branch keep/drop configuration
- `datasets-*.txt` - list of samples (for bkg/data: as they appear on DAS, for sig: path to directory)

## Usage

```bash
# Background/data samples
./runSkimmer.sh -i datasets-2024.txt

# Signal samples
./runSkimmer.sh -i datasets-signal.txt -s
```

Make sure to change `OUTPUT_XRD` in executable to your own directory before running this. 
Make sure the `OUTPUT_TAG` is set in `runSkimmer.sh`, including the skims version, e.g. "v2".  