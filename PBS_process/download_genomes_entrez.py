from Bio import Entrez
import os, shutil, sys
from datetime import datetime
from ftplib import FTP
from helpers import timeit
import traceback

"""
TODO:
- test more to see that it downloads like assembly ncbi
- downloading only latest_refseq
- test Rickettsiales - works on flask, but not on apache
"""

SERVER_PATH = r"ftp.ncbi.nlm.nih.gov"
NUMBER_OF_TRIES = 5  # change to bigger number once stable
count_failed_download = [0]

def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record


def init_dirs(save_folder):
    results_folder = os.path.join(os.path.normpath(save_folder), "results")
    stats_folder = os.path.join(results_folder, os.path.normpath("stats"))
    fasta_folder = os.path.join(results_folder, os.path.normpath("fasta"))

    for dir_path in [results_folder, stats_folder, fasta_folder]:
        os.makedirs(dir_path, exist_ok=True)
        print("created dir %s" % dir_path)

    return results_folder, fasta_folder, stats_folder


# TODO: add try except here and throw exception if failed to download more than X files
def download_files(url_genbank, save_folder, is_fasta=False):
    label = os.path.basename(url_genbank)
    if is_fasta:
        file_ending = 'genomic.fna.gz'
    else:  # stats
        file_ending = 'assembly_stats.txt'
    file2download = f'{url_genbank}/{label}_{file_ending}'.split(SERVER_PATH)[-1]  # get the download link
    local_filename = os.path.join(save_folder, os.path.basename(os.path.normpath(file2download)))

    # don't download if file already been downloaded (could happen bc of try error)
    if os.path.isfile(local_filename):
        return local_filename

    file = open(local_filename, 'wb')
    try:  # creating a new to ftp connection per file (to prevent a timeout)
        ftp = FTP(SERVER_PATH)
        ftp.login()
        print('logged in')
        ftp.retrbinary('RETR ' + file2download, file.write)
        file.close()
        print(f'successfully downloaded file2download: {file2download}')
        try:  # sometimes quits on its own
            ftp.quit()  # This is the “polite” way to close a connection
            print("Successfully quited ftp connection")
        except Exception as e:
            print(f'failed to close ftp connection. Error: {e}')
            print(traceback.format_exc())
    except Exception as e:  # catch failed to connect/ failed to
        count_failed_download[0] += 1
        print(f'error downloading file2download: {file2download}, url_genbank: {url_genbank}, Error: {e}')
        print(e)
        print(traceback.format_exc())
        file.close()
        os.remove(local_filename)
    return local_filename


def get_assembly_ids(term):
    Entrez.email = "estykatzeff@mail.tau.ac.il"
    handle = Entrez.esearch(db="assembly", term=term, retmax='10000')
    record = Entrez.read(handle)
    ids = record['IdList']
    return ids

# Download genbank assemblies for a given search term.
@timeit
def get_assemblies(save_folder, term):
    start_all = datetime.now()
    ids = get_assembly_ids(term)
    print(f'found {len(ids)} ids')
    refseq_assemblies_count = 0

    results_folder, fasta_folder, stats_folder = init_dirs(save_folder)

    for id in ids:
        if count_failed_download[0] > NUMBER_OF_TRIES:
            print(f'passed number of maximum download attemps: {count_failed_download[0]}')
            return
        is_refseq = False
        summary = get_assembly_summary(id)
        url_genbank = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']  # get ftp url
        if url_genbank == '':
            print(f'there is no exisiting assembly for: {url_genbank}')
            continue

        # download stats
        stats_file_path = download_files(url_genbank, stats_folder)

        # checking if assembly is refseq
        if os.path.isfile(stats_file_path):
            with open(stats_file_path) as stats_file:
                for line in stats_file:
                    if line.startswith("# RefSeq assembly accession"):
                        is_refseq = True
                        refseq_assemblies_count += 1

        if is_refseq: # only download assemblies which are refseq (Future: add this as option?)
            fasta_file_path = download_files(url_genbank, fasta_folder, is_fasta=True)  # download assemblies
        else:
            print(f'{url_genbank} is not a refseq. skipping.')

    end = datetime.now()
    print(f'average download time per file: {(end - start_all) / len(ids)}')
    print(f'time: {end - start_all}')
    print(f'refseq_assemblies_count: {refseq_assemblies_count}')

    # zipping the fasta folder
    shutil.make_archive(os.path.join(results_folder, "fasta"), 'zip', fasta_folder)
    print("Successfully archived fasta folder")


def run(save_folder, term, filter_by_level=False, filter_by_date=False):
    get_assemblies(save_folder, term)


def main():
    args = sys.argv[1:]
    run(args[0], args[1], filter_by_level=False, filter_by_date=False)


if __name__ == "__main__":
    main()

# path_to_save_dir = r"C:\Users\97252\Documents\year_4\bio_project_data\download_results"
# bacteria2search = "Xanthomonas campestris"
# # bacteria2search = "Xanthomonas citri"
# run(path_to_save_dir, bacteria2search, filter_by_level=False, filter_by_date=False)