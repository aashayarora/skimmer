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

subprocess.run("python3 -m pip install --user --no-binary=correctionlib correctionlib", shell=True, check=True)
import importlib
correctionlib = importlib.import_module("correctionlib")
correctionlib.register_pyroot_binding()

r.gInterpreter.Declare('#include "truthSelections.h"')
r.gInterpreter.Declare('#include "jetId.h"')

# Constants
CONDOR_OUTPUT_DIR = "output"
XROOTD_REDIRECTOR = "root://xrootd-cms.infn.it/"
OUTPUT_XRD = "davs://redirector.t2.ucsd.edu:1095//store/user/USER_UAF_DIR/skim/" # Change to user's skim directory on UAF
MAX_RETRIES = 10
SLEEP_DURATION = 60  # 1 minute in seconds


class Skimmer():
    def __init__(self, inFiles, outDir, keepDropFile, is_signal):
        self.inFiles = inFiles
        self.outDir = outDir
        self.keepDropFile = keepDropFile
        self.is_signal = is_signal

        self.df = r.RDataFrame("Events", self.inFiles)
        columns = self.df.GetColumnNames()
        for col in columns:
            if col.startswith("Muon_") or col.startswith("Electron_") or col.startswith("Jet_") or col.startswith("FatJet_"):
                colType = self.df.GetColumnType(col)
                if colType == "ROOT::VecOps::RVec<Float_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('Float_t')())
                elif colType == "ROOT::VecOps::RVec<Int_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('Int_t')())
                elif colType == "ROOT::VecOps::RVec<UShort_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('UShort_t')())
                elif colType == "ROOT::VecOps::RVec<Bool_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('Bool_t')())
                elif colType == "ROOT::VecOps::RVec<UChar_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('UChar_t')())
                elif colType == "ROOT::VecOps::RVec<Double_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('Double_t')())
                elif colType == "ROOT::VecOps::RVec<Long64_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('Long64_t')())
                elif colType == "ROOT::VecOps::RVec<Short_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('Short_t')())
                else:
                    print(f"Unknown column type: {colType} for column {col}")

        self.df.Display(["Electron_pt"]).Print()

    @staticmethod
    def genSelection(df):
        # define truth event from gen particles
        df = df.Define("truth_indices", "getTruthEventInfo(GenPart_pdgId, GenPart_status, GenPart_genPartIdxMother, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, false)") \
            .Define("gen_h_idx", "truth_indices[0]") \
            .Define("gen_b1_idx", "truth_indices[1]") \
            .Define("gen_b2_idx", "truth_indices[2]") \
            .Define("gen_v1_idx", "truth_indices[3]") \
            .Define("gen_v1q1_idx", "truth_indices[4]") \
            .Define("gen_v1q2_idx", "truth_indices[5]") \
            .Define("gen_v2_idx", "truth_indices[6]") \
            .Define("gen_v2q1_idx", "truth_indices[7]") \
            .Define("gen_v2q2_idx", "truth_indices[8]") \
            .Define("gen_vbs1_idx", "truth_indices[9]") \
            .Define("gen_vbs2_idx", "truth_indices[10]") 
        
        return df
        
    def analyze(self):
        self.df = self.df.Define("__tight_mu_mask", "Muon_pt > 35. && abs(Muon_eta) < 2.4 && Muon_tightId") \
            .Define("__tight_ele_mask", "Electron_pt > 35. && abs(Electron_eta) < 2.5 && Electron_cutBased >= 4") \
            .Define("__n_tight_leptons", "Sum(__tight_mu_mask) + Sum(__tight_ele_mask)") \
            .Define("__fatjet_mask", "FatJet_pt > 200 && FatJet_msoftdrop > 10") \
            .Define("__n_fatjets", "Sum(__fatjet_mask)") \
            .Filter("(__n_fatjets + __n_tight_leptons) >= 1")

        self.df = self.df.Define("Jet_multiplicity", "Jet_chMultiplicity + Jet_neMultiplicity") \
            .Define("FatJet_multiplicity", "FatJet_chMultiplicity + FatJet_neMultiplicity")

        # Run 2 years use cut-based jet ID (no multiplicity arg, no FatJet equivalent)
        # Run 3 years use JSON-based jet ID (with multiplicity arg)
        year = self.sample_year
        if year in ("2016", "2017", "2018"):
            jet_id_func = f"evalJetID{year if year == '2016' else '20172018'}"
            self.df = self.df \
                .Define("Jet_jetId", f"{jet_id_func}(Jet_eta, Jet_chHEF, Jet_neHEF, Jet_chEmEF, Jet_neEmEF, Jet_muEF, Jet_chMultiplicity, Jet_neMultiplicity)") \
                .Define("FatJet_jetId", f"{jet_id_func}(FatJet_eta, FatJet_chHEF, FatJet_neHEF, FatJet_chEmEF, FatJet_neEmEF, FatJet_muEF, FatJet_chMultiplicity, FatJet_neMultiplicity)")
        else:
            self.df = self.df \
                .Define("Jet_jetId", f"evalJetID{year}(Jet_eta, Jet_chHEF, Jet_neHEF, Jet_chEmEF, Jet_neEmEF, Jet_muEF, Jet_chMultiplicity, Jet_neMultiplicity, Jet_multiplicity)") \
                .Define("FatJet_jetId", f"evalFatJetID{year}(FatJet_eta, FatJet_chHEF, FatJet_neHEF, FatJet_chEmEF, FatJet_neEmEF, FatJet_muEF, FatJet_chMultiplicity, FatJet_neMultiplicity, FatJet_multiplicity)")

        if self.is_signal:
            self.df = self.genSelection(self.df)

        # Run3 event filters
        self.df = self.df.Filter("Flag_goodVertices && "
            "Flag_globalSuperTightHalo2016Filter && "
            "Flag_EcalDeadCellTriggerPrimitiveFilter && "
            "Flag_BadPFMuonFilter && "
            "Flag_BadPFMuonDzFilter && "
            "Flag_hfNoisyHitsFilter &&"
            "Flag_eeBadScFilter && "
            "Flag_ecalBadCalibFilter")

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

    @property
    def sample_year(self):
        fname = self.inFiles[0]
        if re.search(r'Run2016|Summer16|UL2016', fname):
            return "2016"
        elif re.search(r'Run2017|Fall17|UL2017', fname):
            return "2017"
        elif re.search(r'Run2018|Autumn18|UL2018', fname):
            return "2018"
        elif re.search(r'v15', fname):
            return "2024"
        else:
            raise ValueError("Could not determine sample year from filename")


def run_skimmer(input_file, output_dir, is_signal):
    print(f"Running skimmer on {input_file}")
    os.makedirs(output_dir, exist_ok=True)
    
    inFiles = [XROOTD_REDIRECTOR + input_file if input_file.startswith('/store') else 'file://' + input_file]
    keepDropFile = "keep_and_drop_skim.txt"
    
    skimmer = Skimmer(inFiles, output_dir, keepDropFile, is_signal)
    passed = skimmer.analyze()
    if passed:
        skimmer.Snapshot("output")
        return True
    else:
        print("No entries in output")
        return False

def determine_output_paths(input_file, is_signal, output_tag):
    if not is_signal:
        sub_output_dir = "/".join(input_file.split('/')[3:5] + input_file.split('/')[8:])
    else:
        sub_output_dir = "/".join(input_file.split('/')[7:])

    output_dir = f"{OUTPUT_XRD}/skims_{output_tag}/{sub_output_dir}"
    return output_dir

def check_output_liveness(file):
    with r.TFile.Open(file) as f:
        t = f.Get("Events")
        nevts = t.GetEntries()
        for i in range(0,t.GetEntries(),1):
            if t.GetEntry(i) < 0:
                return False
        return True

def copy_output_file(source, destination):
    print(f"Copying {source} to {destination}")
    
    subprocess.run(["gfal-mkdir", "-p", os.path.dirname(destination)])

    # check if output is good
    if not check_output_liveness(source):
        print(f"Output file {source} is corrupted; not copying")
        return False
    
    result = subprocess.run(["gfal-copy", "-f", source, destination])

    if result.returncode != 0:    
        print(f"Failed to copy {source} to {destination}; sleeping for 60s")
        return False

    # check copied file liveness
    if not check_output_liveness(destination):
        print(f"Copied file {destination} is corrupted")
        subprocess.run(["gfal-rm", destination])
        print(f"Removed corrupted file {destination}")
        return False

    return True

if __name__ == "__main__":
    parser = ArgumentParser(description='Run the NanoAOD skimmer with file transfer.')
    parser.add_argument('proxy', help="Path to the X509 proxy")
    parser.add_argument('input_file', help="Input file path")
    parser.add_argument('is_signal', help='Flag indicating if this is a signal sample', type=int)
    parser.add_argument('output_tag', help='Output tag, including version of skims eg. v2', type=str)
    args = parser.parse_args()
    
    os.environ['X509_USER_PROXY'] = args.proxy
    
    success = run_skimmer(args.input_file, CONDOR_OUTPUT_DIR, args.is_signal)
    
    if not success:
        print("Skimmer failed; retrying one more time...")
        success = run_skimmer(args.input_file, CONDOR_OUTPUT_DIR, args.is_signal)

    if not success:
        raise ValueError("Skimmer failed twice; exiting...")

    
    copy_src = os.path.join(os.getcwd(), f"{CONDOR_OUTPUT_DIR}/output.root")
    copy_dest = determine_output_paths(args.input_file, args.is_signal, args.output_tag)
    
    for attempt in range(MAX_RETRIES + 1):
        success = copy_output_file(copy_src, copy_dest)
        if success:
            break
        
        if attempt < MAX_RETRIES:
            print(f"Retrying copy attempt {attempt + 1} of {MAX_RETRIES}...")
            time.sleep(SLEEP_DURATION)

    if not success:
        raise ValueError(f"Failed to copy output file after {MAX_RETRIES} attempts; exiting...")
        
    sys.exit(0)