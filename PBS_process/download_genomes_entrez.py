"""
TODO:
- test more to see that it downloads like assembly ncbi
- downloading only latest_refseq
- test Rickettsiales
"""

from Bio import Entrez
import os, shutil
from datetime import datetime
from ftplib import FTP

SERVER_PATH = r"ftp.ncbi.nlm.nih.gov"
NUMBER_OF_TRIES = 500


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


# def download_fastas(ftp, url_genbank, fasta_folder):
#     label = os.path.basename(url_genbank)
#     assembly_file2download = (url_genbank + "/" + label + '_genomic.fna.gz').split(SERVER_PATH)[
#         -1]  # get the fasta link - change this to get other formats
#     # assembly_file2download = assembly_file2download.split(SERVER_PATH)[-1]
#     local_filename = os.path.join(fasta_folder, os.path.basename(os.path.normpath(assembly_file2download)))
#     file = open(local_filename, 'wb')
#     ftp.retrbinary('RETR ' + assembly_file2download, file.write)
#
#
# def download_stats():
#     label = os.path.basename(url_genbank)
#     assembly_file2download = url_genbank + "/" + label + '_assembly_stats.txt'  # get the fasta link - change this to get other formats
#     assembly_file2download = assembly_file2download.split(SERVER_PATH)[-1]
#     local_filename = os.path.join(stats_folder, os.path.basename(os.path.normpath(assembly_file2download)))
#     file = open(local_filename, 'wb')
#     ftp.retrbinary('RETR ' + assembly_file2download, file.write)


def download_files(ftp, url_genbank, save_folder, is_fasta=False):
    label = os.path.basename(url_genbank)
    if is_fasta:
        file_ending = 'genomic.fna.gz'
    else: # stats ending
        file_ending = 'assembly_stats.txt'
    file2download = f'{url_genbank}/{label}_{file_ending}'.split(SERVER_PATH)[-1]  # get the download link
    local_filename = os.path.join(save_folder, os.path.basename(os.path.normpath(file2download)))
    file = open(local_filename, 'wb')
    print(file2download)
    try:
        ftp.retrbinary('RETR ' + file2download, file.write)
        file.close()
    except Exception as e:
        print(f'error downloading file2download: {file2download}, url_genbank: {url_genbank}, Error: {e}')
        file.close()
        os.remove(local_filename)

# Download genbank assemblies for a given search term.
def get_assemblies(save_folder, term):
    start_all = datetime.now()
    Entrez.email = "estykatzeff@mail.tau.ac.il"
    handle = Entrez.esearch(db="assembly", term=term, retmax='200')
    record = Entrez.read(handle)
    ids = record['IdList']
    print(f'found {len(ids)} ids')
    refseq_assemblies_count = 0

    ftp = FTP(SERVER_PATH)
    ftp.login()
    print('logged in')

    results_folder, fasta_folder, stats_folder = init_dirs(save_folder)

    for id in ids:
        summary = get_assembly_summary(id)

        is_refseq = 'latest_refseq' in summary['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']
        refseq_assemblies_count += int(is_refseq)

        if is_refseq: # only download assemblies which are refseq (Future: add this as option?)
            url_genbank = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']  # get ftp url
            # url_stats_rpt = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_Stats_rpt']  #get ftp url
            # Optional keys are: 'FtpPath_GenBank', 'FtpPath_RefSeq', 'FtpPath_Assembly_rpt', 'FtpPath_Stats_rpt'
            if url_genbank == '':
                print(f'there is no exisiting assembly for id: {id}')
                continue

            print(url_genbank)
            # download assemblies
            download_files(ftp, url_genbank, fasta_folder, is_fasta=True)

            # download stats
            download_files(ftp, url_genbank, stats_folder)

        else:
            print(f'assembly id {id} is not a refseq. skipping.')

    ftp.quit()  # This is the “polite” way to close a connection
    print("Successfully quited ftp connection")

    end = datetime.now()
    print(f'average download time per file: {(end - start_all) / len(ids)}')
    print(f'time: {end - start_all} seconds')
    print(f'refseq_assemblies_count: {refseq_assemblies_count}')

    shutil.make_archive(os.path.join(results_folder, "fasta"), 'zip', fasta_folder)
    print("Successfully archived fasta folder")


def run(save_folder, term, filter_by_level=False, filter_by_date=False):
    get_assemblies(save_folder, term)
    # for i in range(NUMBER_OF_TRIES):
    #     try:
    #         get_assemblies(save_folder, term)
    #     except Exception as e:
    #         print(e)
    #         print(f'finished try number: {i} out of {NUMBER_OF_TRIES}')
    #     else:
    #         print("Exiting")
    #         break


# def main():
#     args = sys.argv[1:]
#     run(args[0], args[1], filter_by_level=False, filter_by_date=False)
#
#
# if __name__ == "__main__":
#     main()

path_to_save_dir = r"C:\Users\97252\Documents\year_4\bio_project_data\download_results"
bacteria2search = "Xanthomonas campestris"
run(path_to_save_dir, bacteria2search, filter_by_level=False, filter_by_date=False)
