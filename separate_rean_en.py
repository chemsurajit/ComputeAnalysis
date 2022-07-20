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


def save_to_txt_file(out_dir, csv_files, columns_to_read):
    """
    This function save the numpy analysis to txt files based on the functional
    with the three bases (except gfnxtb)
    """
    nbins = 50
    chunksize = 1000000
    gmin = {}
    gmax = {}
    for column in columns_to_read:
        gmin[column] = np.inf
        gmax[column] = -np.inf
    gmin = pd.Series(gmin)
    gmax = pd.Series(gmax)
    # first determine the max and min for all the bases or xtb
    for csv_file in csv_files:
        print("File started reading for min, max: ", csv_file)
        for chunk in pd.read_csv(csv_file, 
                                chunksize=chunksize
                                ):
            min_chunk = chunk[columns_to_read].min()
            max_chunk = chunk[columns_to_read].max()
            gmin = np.minimum(min_chunk, gmin)
            gmax = np.maximum(max_chunk, gmax)
        print("File reading complete.")
    gmin = gmin.to_dict()
    gmax = gmax.to_dict()
    print("MINIMA:")
    [print(k, ":", v) for k, v in gmin.items()]
    print("MAXIMA:")
    [print(k, ":", v) for k, v in gmax.items()]
    # now calculate the bins and counts.
    bin_edges_dict = {}
    counts_dict = {}
    for column in columns_to_read:
        bin_edges_dict[column] = np.linspace(gmin[column], gmax[column], nbins+1)
        counts_dict[column] = np.zeros(nbins, np.int64)
    print("Created the empty numpy arrays for the functionals.")
    print("Now the counting part will be performed.")
    # calculate the bins
    for csv_file in csv_files:
        print("File started reading for hist count: ", csv_file)
        nchunk = 0
        for chunk in pd.read_csv(csv_file, chunksize=chunksize):
            selected_chunk = chunk[columns_to_read]
            for column in columns_to_read:
                subtotal_col, e_col = np.histogram(selected_chunk[column], bins=bin_edges_dict[column])
                counts_dict[column] += subtotal_col
            nchunk += 1
            if nchunk % 10 == 0:
                print("Nchunk read: ", nchunk)
        print("Analysis complete for file: ", csv_file)
    # Now saving the results to txt files
    print("Now saving the results to txt files.")
    for column in columns_to_read:
        count_file = column + "_counts.txt"
        binedges_file = column + "_binedges.txt"
        np.savetxt(count_file, counts_dict[column], fmt="%d")
        np.savetxt(binedges_file, bin_edges_dict[column], fmt="%f")
    print("Saving file is completed.")
    #         if functional.strip() == "GFNXTB":
    #            subtotal_xtb, e_xtb = np.histogram(Chunk[functional], bins=bin_edges_xtb)
    #            count_xtb += subtotal_xtb
    #         else:
    #            subtotal_tzp, e_tzp = np.histogram(chunk[functional+"_TZP"], bins=bin_edges_tzp)
    #            subtotal_dzp, e_dzp = np.histogram(chunk[functional+"_DZP"], bins=bin_edges_dzp)
    #            subtotal_sz, e_sz = np.histogram(chunk[functional+"_SZ"], bins=bin_edges_sz)
    #            count_tzp += subtotal_tzp
    #            count_dzp += subtotal_dzp
    #            count_sz += subtotal_sz
    #     print("File reading complete.")
    # print("Histogram for all data complete.")
    # if functional.strip() == "GFNXTB":
    #     np.savetxt(functional+"_counts.txt", count_xtb, fmt="%d")
    #     np.savetxt(functional+"_binedges.txt", bin_edges_xtb, fmt="%f")
    # else:
    #     np.savetxt(functional+"_counts_TZP.txt", count_tzp, fmt="%d")
    #     np.savetxt(functional+"_binedges_TZP.txt", bin_edges_tzp, fmt="%f")
    #     np.savetxt(functional+"_counts_DZP.txt", count_dzp, fmt="%d")
    #     np.savetxt(functional+"_binedges_DZP.txt", bin_edges_dzp, fmt="%f")
    #     np.savetxt(functional+"_counts_SZ.txt", count_sz, fmt="%d")
    #     np.savetxt(functional+"_binedges_SZ.txt", bin_edges_sz, fmt="%f")
    # print("Histogram analysis complete for functional: ", functional)
    return


def save_reaction_histogram(out_dir, csv_files, dft_functional_names):
    """
    This is the main function which will save the histogram analysis of the dft functionals names.
    """
    columns_to_read = []
    for functional in dft_functional_names:
        if functional.strip() == "GFNXTB":
            columns_to_read.append(functional)
        else:
            columns_to_read += [functional+"_TZP", functional+"_DZP", functional+"_SZ"]
    print("The following columns will be analysed: ", *functional, sep="\n")
    save_to_txt_file(out_dir, csv_files, columns_to_read)
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
