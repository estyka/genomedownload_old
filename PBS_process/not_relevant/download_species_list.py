from ftplib import FTP
import pickle
import sys
import os

NUMBER_OF_TRIES = 500

SERVER_PATH = r"ftp.ncbi.nlm.nih.gov"
SERVER_BACTERIA_LOCATION = r"/genomes/refseq/bacteria/"


def download_species_list(output_path):
    print('starting create_process_update_bacteria_list')
    # for i in range(NUMBER_OF_TRIES):
    #     try:
    #         print(f'starting num_try {i}')

    ftp = FTP(SERVER_PATH)
    ftp.login()
    print(f'logged in {ftp}')

    ftp.cwd(SERVER_BACTERIA_LOCATION)
    print(f'changed dir to {SERVER_BACTERIA_LOCATION}')

    all_ncbi_bacterias = ftp.nlst(SERVER_BACTERIA_LOCATION)
    ftp.quit()  # This is the “polite” way to close a connection
    print("Successfully quited ftp connection")

    print(f'updating bacteria list in {output_path}')
    with open(output_path, 'wb') as f:  # TODO: make sure this overwrites existing data
        pickle.dump(all_ncbi_bacterias, f)
    print(f'finished updating bacteria list')

        # except Exception as e:
        #     print(e)
        #     print(f'finished try number: {i} out of {NUMBER_OF_TRIES}')
        #
        # else:
        #     print("Exiting")
        #     break

            
def main():
    args = sys.argv[1:]
    download_species_list(args[0])


if __name__ == "__main__":
    main()