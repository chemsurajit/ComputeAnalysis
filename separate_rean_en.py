import argparse
from chunk import Chunk
import os
import pandas as pd
import numpy as np


def get_arguments():
    parser = argparse.ArgumentParser(
        description="Extract the reaction energies (all 3 basis) from the csv files specific to a DFT functiona."
    )
    #parser.add_argument(
    #    "-dft_functional", "--dft_functional",
    #    type=str,
    #    required=True,
    #    help="Name of the DFT functional which energy need to be extracted."
    #)
    parser.add_argument(
        "-csv_directory", "--csv_directory",
        type=str,
        required=True,
        help="Directory where all the CSV files are saved."
    )
    parser.add_argument(
        "-out_dir", "--out_dir",
        type=str,
        required=False,
        default="./",
        help="Directory where the extracted energies will be saved."
    )
    return parser.parse_args()


def get_csv_files(directory):
    """
    Return a list of csv files with absolute path from the directory.
    """
    print("csv files located in: ", directory)
    csv_files = []
    for files in os.listdir(directory):
        abs_filepath = os.path.abspath(os.path.join(directory, files))
        if os.path.isfile(abs_filepath):
            if abs_filepath.endswith(".csv"):
                csv_files.append(abs_filepath)
    print("csv files:")
    print(*csv_files, sep='\n')
    return csv_files


def save_to_txt_file(out_dir, csv_files, functional):
    """
    This function save the numpy analysis to txt files based on the functional
    with the three bases (except gfnxtb)
    """
    nbins = 50
    chunksize = 1000000
    print("Processing data for: %s" % functional)
    if functional.strip() == "GFNXTB":
        column_names = ["GFNXTB"]
    else:
        column_names = [
            functional + "_SZ",
            functional + "_DZP",
            functional + "_TZP"
        ]
    min_tzp = np.inf
    max_tzp = -np.inf
    min_dzp = np.inf
    max_dzp = -np.inf
    min_sz = np.inf
    max_sz = -np.inf
    min_xtb = np.inf
    max_xtb = -np.inf
    # first determine the max and min for all the bases or xtb
    for csv_file in csv_files:
        print("File started reading for min, max: ", csv_file)
        for chunk in pd.read_csv(csv_file, 
                                chunksize=chunksize
                                ):
            if functional.strip() == "GFNXTB":
                min_xtb = np.minimum(chunk[column_names].min(), min_xtb)
                max_xtb = np.maximum(chunk[column_names].max(), max_xtb)
            else:
                min_tzp = np.minimum(chunk[functional + "_TZP"].min(), min_tzp)
                max_tzp = np.maximum(chunk[functional + "_TZP"].max(), max_tzp)
                min_dzp = np.minimum(chunk[functional + "_DZP"].min(), min_dzp)
                max_dzp = np.maximum(chunk[functional + "_DZP"].max(), max_dzp)
                min_sz = np.minimum(chunk[functional + "_SZ"].min(), min_sz)
                max_sz = np.maximum(chunk[functional + "_SZ"].max(), max_sz)
        print("File reading complete.")
    print("The min, max for TZP: " , min_tzp, max_tzp)
    print("The min, max for DZP: ", min_dzp, max_dzp)
    print("The min, max for SZ: ", min_sz, max_sz)
    print("The min, max for XTB: ", min_xtb, max_xtb)
    # now calculate the bins and counts.
    if functional.strip() == "GFNXTB":
        bin_edges_xtb = np.linspace(min_xtb, max_xtb, nbins+1)
        count_xtb = np.zeros(nbins, np.int32)
    else:
        bin_edges_tzp = np.linspace(min_tzp, max_tzp, nbins+1)
        count_tzp = np.zeros(nbins, np.int32)
        bin_edges_dzp = np.linspace(min_dzp, max_dzp, nbins+1)
        count_dzp = np.zeros(nbins, np.int32)
        bin_edges_sz = np.linspace(min_sz, max_sz, nbins+1)
        count_sz = np.zeros(nbins, np.int32)
    # calculate the bins
    for csv_file in csv_files:
        print("File started reading for hist count: ", csv_file)
        for chunk in pd.read_csv(csv_file, chunksize=chunksize):
            if functional.strip() == "GFNXTB":
                subtotal_xtb, e_xtb = np.histogram(Chunk[functional], bins=bin_edges_xtb)
                count_xtb += subtotal_xtb
            else:
                subtotal_tzp, e_tzp = np.histogram(chunk[functional+"_TZP"], bins=bin_edges_tzp)
                subtotal_dzp, e_dzp = np.histogram(chunk[functional+"_DZP"], bins=bin_edges_dzp)
                subtotal_sz, e_sz = np.histogram(chunk[functional+"_SZ"], bins=bin_edges_sz)
                count_tzp += subtotal_tzp
                count_dzp += subtotal_dzp
                count_sz += subtotal_sz
        print("File reading complete.")
    print("Histogram for all data complete.")
    if functional.strip() == "GFNXTB":
        np.savetxt(functional+"_counts.txt", count_xtb, fmt="%d")
        np.savetxt(functional+"_binedges.txt", bin_edges_xtb, fmt="%f")
    else:
        np.savetxt(functional+"_counts_TZP.txt", count_tzp, fmt="%d")
        np.savetxt(functional+"_binedges_TZP.txt", bin_edges_tzp, fmt="%f")
        np.savetxt(functional+"_counts_DZP.txt", count_dzp, fmt="%d")
        np.savetxt(functional+"_binedges_DZP.txt", bin_edges_dzp, fmt="%f")
        np.savetxt(functional+"_counts_SZ.txt", count_sz, fmt="%d")
        np.savetxt(functional+"_binedges_SZ.txt", bin_edges_sz, fmt="%f")
    print("Histogram analysis complete for functional: ", functional)
    return


def save_reaction_histogram(out_dir, csv_files, dft_functional_names):
    """
    This function will save the histogram analysis of the dft functionals names.
    """
    for functional in dft_functional_names:
        print("Doing histogram for functional: ", functional)
        save_to_txt_file(out_dir, csv_files, functional)
        print("Saved histogram analysis for functional: ", functional)
    return


def main():
    args = get_arguments()
    dft_functional_names = ["PBE", "BLYP", "BP", "TPSSH", "B3LYP(VWN5)", "M06-2X", "MPW1PW", "GFNXTB"]
    csv_files = get_csv_files(args.csv_directory)
    save_reaction_histogram(args.out_dir, csv_files, dft_functional_names)
    print("All finished")
    return


if __name__ == "__main__":
    main()
