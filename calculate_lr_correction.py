import argparse
import csv
import os
import glob
import pandas as pd


def get_csv_files(csv_dir):
    """This function returns a list of csv files with matching pattern Reaction_*.csv."""
    files = []
    cwd = os.getcwd()
    os.chdir(csv_dir)
    for ifile in glob.glob('Reactions_*.csv'):
        files.append(os.path.abspath(ifile))
    os.chdir(cwd)
    return files


def get_energy_data(csv_files, dft_functional, bonds_coeffs_dict = None):
    """Returns:
               dE, dE_corrected
    """
    dE = []
    dE_corrected = []
    reactindex = []
    pdtindex = []
    chunksize = 100000
    bonds_list = list(bonds_coeffs_dict.keys())
    print("bonds list: ", bonds_list)
    print("Starting to read csv files")
    for csvs in csv_files:
        nchunk =0
        print("Reading csv file: ", csvs)
        for chunk in pd.read_csv(csvs, chunksize=chunksize):
            nchunk += 1
            df = chunk.dropna(axis=0, how='any')
            # remove all the reactions where there is no bond changes:
            df = df.loc[(df[bonds_list].abs().sum(axis=1) != 0)]
            dE_np = df.loc[:, ("G4MP2")] - df.loc[:, (dft_functional.upper())].to_numpy()
            dE += dE_np.to_list()
            reactindex += df["reactindex"].to_list()
            pdtindex += df["pdtindex"].to_list()
            for bond, coef in bonds_coeffs_dict.items():
                bond_value = df[bond].to_numpy()
                dE_np -= coef * bond_value
            dE_corrected += dE_np.to_list()
            if nchunk%10 == 0:
                print("Nchunk: ", nchunk)
    return reactindex, pdtindex, dE, dE_corrected


def save_to_csv(dE, dE_corrected, reactindex, pdtindex, dft_functional, output=None):
    print("Writing energy correction data to csv file: ", output)
    with open(output, 'w') as fp:
        writer = csv.writer(fp)
        writer.writerow(["reactindex", "pdtindex", "dE_" + str(dft_functional), "dE_corrected_" + str(dft_functional)])
        writer.writerows(zip(reactindex, pdtindex, dE, dE_corrected))
    return


def main():
    parser = argparse.ArgumentParser("Program to plot the correction from the result of linear regression.")
    parser.add_argument(
        '-csv_directory', '--csv_directory',
        help="Directory where the Reactions_n.csv files are present",
        required=True
    )
    parser.add_argument(
        '-lr_coeff', '--lr_coeff',
        help="CSV file that contains the coefficients from LR. Format: bonds,coeffs",
        required=True
    )
    parser.add_argument(
        '-dft_functional', '--dft_functional',
        help="Name of the csv function for which the LR was performed.",
        choices=["PBE", "B3LYP-D", "M06-2X", "GFNXTB"],
        required=True
    )
    parser.add_argument(
        '-name', '--name',
        help="Name of this run. Output files will be generated based on this name.",
        required=True
    )
    args = parser.parse_args()
    output = args.name
    dft_functional = args.dft_functional
    csv_files = get_csv_files(args.csv_directory)
    bonds_lr_pd = pd.read_csv(args.lr_coeff, index_col=False)
    bonds_coeffs_dict = dict(zip(bonds_lr_pd.bonds, bonds_lr_pd.coefficients))
    reactindex, pdtindex, dE, dE_corrected = get_energy_data(csv_files, dft_functional, bonds_coeffs_dict=bonds_coeffs_dict)
    save_to_csv(dE, dE_corrected, reactindex, pdtindex, dft_functional, output=output+".csv")
    pass


if __name__ == "__main__":
    main()
