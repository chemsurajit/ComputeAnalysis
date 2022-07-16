import argparse
import os
import pandas as pd


def get_arguments():
    parser = argparse.ArgumentParser(
        description="Extract the reaction energies (all 3 basis) from the csv files specific to a DFT functiona."
    )
    parser.add_argument(
        "-dft_functional", "--dft_functional",
        type=str,
        required=True,
        help="Name of the DFT functional which energy need to be extracted."
    )
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
        abs_filepath = os.path.abspath(files)
        if os.path.isfile(abs_filepath):
            if abs_filepath.endswith(".csv"):
                csv_files.append(abs_filepath)
    print("csv files:")
    print(csv_files)
    return csv_files


def save_reaction_energies(out_dir, csv_files, dft_functional_name):
    if "G4MP2" in dft_functional_name:
        column_names = ["G4MP2"]
    else:
        column_names = [dft_functional_name + "_TZP", 
                    dft_functional_name+"_DZP", 
                    dft_functional_name+"_SZ"]
    print("DFT functional name: ", dft_functional_name)
    print("Now will save the energy values.")
    chunksize = 100000
    out_csv = os.path.join(out_dir, dft_functional_name+".csv")
    n = 0
    for csv_file in csv_files:
        print("CSV writing started with file: ", csv_file)
        for chunk in pd.read_csv(csv_file, usecols=column_names, chunksize=chunksize):
            if chunk:
                if n == 0:
                    chunk.to_csv(out_csv, mode="w", index=False)
                else:
                    chunk.to_csv(out_csv, mode="a", index=False, header=False)
                n += 1
            if n%10 == 0:
                print("Nchunk: ", n)
        print("Writing complete for: ", csv_file)
    return


def main():
    args = get_arguments()
    dft_functional_name = args.dft_functional
    csv_files = get_csv_files(args.csv_directory)
    save_reaction_energies(args.out_dir, csv_files, dft_functional_name)
    return


if __name__ == "__main__":
    main()