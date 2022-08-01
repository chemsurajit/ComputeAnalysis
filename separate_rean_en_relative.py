import argparse
from chunk import Chunk
from genericpath import isfile
import os,sys
import pandas as pd
import numpy as np
import re
import separate_rean_en



def get_arguments():
    parser = argparse.ArgumentParser(
        description="Extract the relative reaction energies."
    )
    parser.add_argument(
        "-csv_dir", "--csv_dir",
        type=str,
        required=True,
        help="Directory where all the CSV files are saved."
    )
    parser.add_argument(
        "-g4mp2_csv", "--g4mp2_csv",
        required=False,
        default="./",
        help="CSV where G4MP2 reaction energies are. Required only the -relative present."
    )
    parser.add_argument(
        "-out_dir", "--out_dir",
        type=str,
        required=False,
        default="./",
        help="Directory where the extracted energies will be saved."
    )
    return parser.parse_args()


def save_relative_csv_for_functionals(csv_files, g4mp2_csv, dft_functional_names):
    chunksize = 10000
    for func in dft_functional_names:
        func_dfs = []
        for csv_file in csv_files:
            for chunk in pd.read_csv(csv_file, chunksize=chunksize, usecols=[func]):
                func_dfs.append(chunk)
            break
        break
    pd.concat(func_dfs).to_csv("test.csv")
    return



def main():
    nbins = 50
    args = get_arguments()
    csv_files = separate_rean_en.get_csv_files(args.csv_dir)
    dft_functional_names = ["PBE_TZP", "PBE_DZP", "PBE_SZ", "B3LYP_TZP", "M06-2X_TZP", "GFNXTB"]
    if os.path.isfile(args.g4mp2_csv):
        save_relative_csv_for_functionals(csv_files, 
                                          args.g4mp2_csv, 
                                          dft_functional_names, 
                )
    print("All finished")
    return


if __name__ == "__main__":
    main()