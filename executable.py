#!/usr/bin/python3
# -*- coding: utf-8 -*-

import re
import os
import sys
import glob
import time
import shutil
import subprocess
from argparse import ArgumentParser

import ROOT as r

# Constants
CONDOR_OUTPUT_DIR = "output"
XROOTD_REDIRECTOR = "root://xrootd-cms.infn.it/"
OUTPUT_XRD = "davs://redirector.t2.ucsd.edu:1095//store/user/aaarora/skims"
CMSSW_VERSION = 'CMSSW_14_1_4'
MAX_RETRIES = 10
SLEEP_DURATION = 60  # 1 minute in seconds

class Skimmer():
    def __init__(self, inFiles, outDir, keepDropFile, n_workers=1):
        self.inFiles = inFiles
        self.outDir = outDir
        self.keepDropFile = keepDropFile
        r.EnableImplicitMT(n_workers)

        self.df = r.RDataFrame("Events", self.inFiles)
    
    def analyze(self):
        self.df = self.df.Define("tight_mu_mask", "Muon_pt > 35. && abs(Muon_eta) < 2.4 && Muon_tightId") \
            .Define("tight_ele_mask", "Electron_pt > 35. && abs(Electron_eta) < 2.5 && Electron_cutBased >= 4") \
            .Filter("Sum(tight_mu_mask) + Sum(tight_ele_mask) < 2") \
            .Define("fatjet_mask", "FatJet_pt > 200 && FatJet_msoftdrop > 10 && FatJet_mass > 10") \
            .Filter("(Sum(fatjet_mask) >= 1)")

        return self.df.Count().GetValue()

    def Snapshot(self, tag):
        all_cols = [str(col) for col in self.df.GetColumnNames()]
        keep_cols = {col: 0 for col in all_cols}
        comment = re.compile(r"#.*")
        ops = []
        with open(self.keepDropFile, 'r') as f:
            for line in f:
                # convert to python regex
                if len(line) == 0 or line[0] == '#': 
                    continue
                line = re.sub(comment, "", line)
                (op, sel) = line.split()
                if op == "keep":
                    ops.append((sel, 1))
                elif op == "drop":
                    ops.append((sel, 0))
        
        for bre, stat in ops:
            try:
                re.compile(bre)
                for n in all_cols:
                    if re.match(bre, n):
                        keep_cols[n] = stat
            except re.error:
                keep_cols[bre] = stat

        keep_cols = [k for k, v in keep_cols.items() if v == 1]

        self.df.Snapshot("Events", self.outDir + "/" + tag + ".root", keep_cols)
        
        snap_opts = r.RDF.RSnapshotOptions()
        snap_opts.fMode = "UPDATE"

        runs_df = r.RDataFrame("Runs", self.inFiles)
        cols = [str(col) for col in runs_df.GetColumnNames()]
        runs_df.Snapshot("Runs", self.outDir + "/" + tag + ".root", cols, snap_opts)


def run_skimmer(input_file, output_dir, n_workers):
    print(f"Running skimmer on {input_file}")
    os.makedirs(output_dir, exist_ok=True)
    
    inFiles = [XROOTD_REDIRECTOR + input_file if input_file.startswith('/store') else 'file://' + input_file]
    keepDropFile = "keep_and_drop_skim.txt"
    
    skimmer = Skimmer(inFiles, output_dir, keepDropFile, n_workers=n_workers)
    passed = skimmer.analyze()
    if passed:
        skimmer.Snapshot("skim")
        return True
    else:
        print("No entries in output")
        return False


def merge_skims(output_dir):
    skim_files = glob.glob(f"{output_dir}/*")
    
    if len(skim_files) == 0:
        print("No output files to merge; exiting...")
        return True
    elif len(skim_files) == 1:
        shutil.move(skim_files[0], f"{output_dir}/output.root")
        return True
    else:
        merge_cmd = ["hadd", f"{output_dir}/output.root"] + skim_files
        print(" ".join(merge_cmd))
        result = subprocess.run(merge_cmd)
        return result.returncode == 0


def determine_output_paths(input_file, is_signal):
    if not is_signal:
        era = input_file.split('/')[3]
        sample_name = input_file.split('/')[4]
        campaign = input_file.split('/')[6]
    else:
        era = input_file.split('/')[6]
        sample_name = input_file.split('/')[7]
        campaign = "private"
        
    output_dir = f"{OUTPUT_XRD}/{era}/{campaign}/{sample_name}"
    return output_dir


def copy_output_file(source, destination):
    print(f"Copying {source} to {destination}")
    
    # Create destination directory
    subprocess.run(["gfal-mkdir", "-p", os.path.dirname(destination)])
    
    # Copy with retries
    for i in range(1, MAX_RETRIES + 1):
        print(f"Attempt {i}")
        result = subprocess.run(["gfal-copy", "-f", source, destination])
        if result.returncode == 0:
            return True
        
        print(f"Failed to copy {source} to {destination}; sleeping for 60s")
        time.sleep(SLEEP_DURATION)
    
    return False


if __name__ == "__main__":
    parser = ArgumentParser(description='Run the NanoAOD skimmer with file transfer.')
    parser.add_argument('proxy', help="Path to the X509 proxy")
    parser.add_argument('input_file', help="Input file path")
    parser.add_argument('job_id', help="Job ID")
    parser.add_argument('is_signal', help='Flag indicating if this is a signal sample', type=int)
    parser.add_argument('nworkers', help='Number of workers to use', type=int)
    args = parser.parse_args()
    
    # Set up X509 proxy
    os.environ['X509_USER_PROXY'] = args.proxy

    # Run the skimmer
    success = run_skimmer(args.input_file, CONDOR_OUTPUT_DIR, args.nworkers)
    
    # Retry once if failed
    if not success:
        print("Skimmer failed; retrying one more time...")
        success = run_skimmer(args.input_file, CONDOR_OUTPUT_DIR, args.nworkers)
    
    # Merge results
    merge_skims(CONDOR_OUTPUT_DIR)
    
    # Determine output paths
    output_dir = determine_output_paths(args.input_file, args.is_signal)
    
    # Copy the output file
    copy_src = os.path.join(os.getcwd(), f"{CONDOR_OUTPUT_DIR}/output.root")
    copy_dest = f"{output_dir}/output_{args.job_id}.root"
    
    success = copy_output_file(copy_src, copy_dest)
    if not success:
        print(f"Failed to copy output file after {MAX_RETRIES} attempts")
        sys.exit(1)
    
    sys.exit(0)