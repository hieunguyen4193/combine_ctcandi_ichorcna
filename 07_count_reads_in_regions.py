import pandas as pd
import numpy as np
import pathlib 
import os
import matplotlib.pyplot as plt
import seaborn as sns
import random
from tqdm import tqdm
import warnings
import pandas as pd
import argparse

from dataset_paths import *

warnings.filterwarnings('ignore')

def main(args):
    output_version = args.output_version
    dataset_name = args.dataset_name
    input_cancer_class = args.input_cancer_class
    outdir = args.outdir
    mode = args.mode

    print(f"output_version: {output_version}")
    print(f"dataset_name: {dataset_name}")
    print(f"input_cancer_class: {input_cancer_class}")
    print(f"outdir: {outdir}")
    print(f"mode: {mode}")

    PROJECT = "combine_ctcandi_ichorcna"
    thres_hypo = 0.3
    thres_hyper = 0.6

    path_to_input = dataset_paths[dataset_name]

    path_to_main_output = os.path.join(outdir, PROJECT, output_version, dataset_name)
    if mode == "all":
       path_to_03_output = os.path.join(outdir, PROJECT, output_version, "03_output", input_cancer_class)
    elif mode == "hypo":
       path_to_03_output = os.path.join(outdir, PROJECT, output_version, "03_output_all_hypo", input_cancer_class)
    elif mode == "hyper":
        path_to_03_output = os.path.join(outdir, PROJECT, output_version, "03_output_all_hyper", input_cancer_class)
    else: 
        raise ValueError("mode must be either all, hypo or hyper")

    path_to_07_output = os.path.join(path_to_main_output, f"07_output_{mode}", input_cancer_class)
    os.system(f"mkdir -p {path_to_07_output}")

    regiondf = pd.read_excel(os.path.join(path_to_03_output, "countDMPs.xlsx"))
    regiondf["hypo_or_hyper"] = regiondf[["hyper", "hypo"]].apply(lambda x: "hyper" if x[0] > x[1] else "hypo", axis = 1)

    #####---------------------------------------------------------------------#####
    ##### Helper functions
    #####---------------------------------------------------------------------#####
    def assign_read_type(x, thres_hypo, thres_hyper):
        if x < thres_hypo:
            return "hypo"
        elif x > thres_hyper:
            return "hyper"
        else:
            return "none"
    def check_read_inside_region(start, seq, region):
            read_end = start + len(seq)
            region_start = int(region.split(":")[1].split("-")[0])
            region_end = int(region.split(":")[1].split("-")[1])
            if start >= region_start and read_end <= region_end:
                return "in"
            else: 
                return "overlap"

    all_read_files = [item for item in pathlib.Path(path_to_input).glob("*.csv")]
    print(f"Number of sample in this dataset: {len(all_read_files)}")

    for file in tqdm(all_read_files):
        tmpdf = pd.read_csv(file, index_col = [0])
        tmpdf["read_overlap_rate"] = tmpdf[["start", "seq", "region"]].apply(lambda x: check_read_inside_region(x[0], x[1], x[2]), axis = 1)

        ##### keep only reads that are completely inside the region
        tmpdf = tmpdf[tmpdf["read_overlap_rate"] == "in"]

        ##### assign read type: hyper or hypo reads based on the given thresholds
        tmpdf["read_classification"] = tmpdf["alpha"].apply(lambda x: assign_read_type(x, thres_hypo, thres_hyper))

        ##### considers only regions that are tested with the TCGA data
        tmpdf["region"] = tmpdf["region"].apply(lambda x: x.replace(":", "_").replace("-", "_"))
        tmpdf = tmpdf[tmpdf["region"].isin(regiondf.Var1.unique())]

        ##### count hypo and hyper reads in each region
        resdf = tmpdf.groupby(["region", "read_classification"]).seq.count().reset_index().pivot_table(index = "region", columns = "read_classification", values = "seq").reset_index().fillna(0)

        ##### get the region type from TCGA test results
        resdf["region_type"] = resdf["region"].apply(lambda x: regiondf[regiondf.Var1 == x].hypo_or_hyper.values[0])

        ##### assign candi reads equal to number of hypo or hyper reads, depending on the region type
        resdf["candi_reads"] = resdf[["region_type", "hyper", "hypo"]].apply(lambda x: x[1] if x[0] == "hyper" else x[2], axis = 1)

        ##### save the results
        resdf.to_csv(os.path.join(path_to_07_output, "{}.candi_reads.csv".format(file.name.split(".")[0])), index = False)
        tmpdf.to_csv(os.path.join(path_to_07_output, file.name.replace(".sorted.csv", ".read_classification.csv")))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process read data and classify reads.")
    parser.add_argument("--output_version", type=str, required=True, help="output version of the analysis")
    parser.add_argument("--dataset_name", type=str, required=True, help="Name of the dataset (LOD, spike in, report 4, validation, ...)")
    parser.add_argument("--outdir", type=str, required=True, help="path to main output directory")
    parser.add_argument("--input_cancer_class", type=str, required=True, help="input cancer class")
    parser.add_argument("--mode", type=str, required=True, help="choose all or hypo or hyper only")
    
    args = parser.parse_args()
    
    main(args)
    
##### Example command:
# output_version="20241229";
# outdir="/media/hieunguyen/HNSD_mini/outdir";
# for mode in all hypo hyper;do \
#     for input_cancer_class in Liver Lung Breast CRC pan_cancer;do \
#         for dataset_name LOD SPIKE_IN REPORT4 VALIDATION;do \
#             echo -e "Working on mode: " $mode ", single/multi cancer tpe: " $input_cancer_class ", dataset " $dataset_name;
#             python 07_count_reads_in_regions.py --output_version $output_version --dataset_name $dataset_name --input_cancer_class $input_cancer_class --outdir $outdir --mode $mode;\
#                 done;done;done

