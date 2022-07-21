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
        "-g4mp2", "--g4mp2",
        action="store_true",
        help="If the histogram for only the G4MP2 is need. Don't run for other functional."
    )
    parser.add_argument(
        "-out_dir", "--out_dir",
        type=str,
        required=False,
        default="./",
        help="Directory where the extracted energies will be saved."
    )
    return parser.parse_args()


def save_reaction_histogram_g4mp2(out_dir, csv_files, nbins=None):
    """
    function which saves the histograms for the G4MP2 reaction energies from csv files.
    The csv file needs to have a column name G4MP2
    """
    if nbins is None:
        nbins = 50
    chunksize = 1000000
    gmin = np.inf
    gmax = -np.inf
    nreactions = 0
    for csv_file in csv_files:
        print("File started reading for min, max: ", csv_file)
        for chunk in pd.read_csv(csv_file, chunksize=chunksize):
            nreactions += len(chunk.index)
            min_chunk = chunk["G4MP2"].min()
            max_chunk = chunk["G4MP2"].max()
            gmin = np.minimum(min_chunk, gmin)
            gmax = np.maximum(max_chunk, gmax)
        print("File reading complete")
    print("min and max value of G4MP2:", gmin, gmax)
    print("Total number of reactions: ", nreactions)
    print("Now the counting part will be performed.")
    bin_edges = np.linspace(gmin, gmax, nbins+1)
    counts = np.zeros(nbins, np.int64)
    for csv_file in csv_files:
        nchunk = 0
        for chunk in pd.read_csv(csv_file, chunksize=chunksize):
            subtotal_count, e_col = np.histogram(chunk["G4MP2"], bins=bin_edges)
            counts += subtotal_count
            nchunk += 1
            if nchunk%10 == 0:
                print("Nchunk: ", nchunk)
        print("File analysis complete")
    print("Read all the csv files. Now save the data to txt file")
    np.savetxt("G4MP2_counts.txt", counts, fmt="%d")
    np.savetxt("G4MP2_binedges.txt", bin_edges, fmt="%f")
    print("Saving file complete.")
    print("Analysis complete.")
    return


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


def save_to_txt_file(out_dir, csv_files, columns_to_read, nbins=None):
    """
    This function save the numpy analysis to txt files based on the functional
    with the three bases (except gfnxtb)
    """
    if nbins is None:
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
    return


def save_reaction_histogram(out_dir, csv_files, dft_functional_names):
    """
    This is the main function which will save the histogram analysis of the dft functionals names.
    """
    columns_to_read = []
    for functional in dft_functional_names:
        if functional.strip() == "GFNXTB":
            columns_to_read += [functional]
        else:
            columns_to_read += [functional+"_TZP", functional+"_DZP", functional+"_SZ"]
    print("The following columns will be analysed: ")
    print(*columns_to_read, sep="\n")
    save_to_txt_file(out_dir, csv_files, columns_to_read)
    return


def main():
    nbins = 50
    args = get_arguments()
    dft_functional_names = ["PBE", "BLYP", "BP", "TPSSH", "B3LYP(VWN5)", "M06-2X", "MPW1PW", "GFNXTB"]
    csv_files = get_csv_files(args.csv_directory)
    if args.g4mp2:
        save_reaction_histogram_g4mp2(args.out_dir, csv_files, nbins=nbins)
    else:
        save_reaction_histogram(args.out_dir, csv_files, dft_functional_names, nbins=nbins)
    print("All finished")
    return


if __name__ == "__main__":
    main()
