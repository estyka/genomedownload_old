# https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Xanthomonas_oryzae/representative/GCF_001021915.1_ASM102191v1/
"""
next steps:
- check if the downloading time from ncbi is similar time - if not - see if there is a way of downloading the whole folder
- make sure the codes comes to an end (quits the ftp connection) - works on bacteria with small database, but for some reason doesn't end for x.oryzae. Even though downloads everything (I think) - it gets errors at the end of the run...
-> idea of how to check this out - count all assemblies and check with debug while seeing what happens when gets to the designated number
- connect the stats code - to download only relevant ones...
#example: python .\download_files_from_ncbi.py C:\Users\97252\Documents\year_4\bio_project_data\download_results Xanthomonas_albilineans
"""

from ftplib import FTP
import os.path
import os, sys, time, re
import pandas as pd
from datetime import datetime

SERVER_PATH = r"ftp.ncbi.nlm.nih.gov"
SERVER_BACTERIA_LOCATION = r"/genomes/refseq/bacteria/"


ALL_ASSEMBLY_PREFIX = r"latest_assembly_versions/"
NUMBER_OF_TRIES = 500
NUM_OF_ERRORS = 50
ASSEMBLY_LEVEL_DIC = {"complete genome": 1, "chromosome": 2, "scaffold": 3, "contig": 4}  # the smaller the number the better


def get_assemblies_per_org_dict(stats_folder):
    org_dict = {}
    for stats_filename in os.listdir(stats_folder):
        if stats_filename.endswith("stats.txt"):
            cur_dict = {}
            with open(os.path.join(stats_folder, stats_filename)) as stats_file:
                for line in stats_file:
                    if line.startswith("# Organism name:"):
                        cur_dict["organism_name"] = line.split(":")[-1].lstrip().strip("\n")
                    elif line.startswith("# Infraspecific name:  strain="):
                        cur_dict["strain"] = line.split("strain=")[-1].strip("\n")
                    elif line.startswith("# Date:"):
                        cur_dict["date"] = line.split()[-1].strip("\n")
                    elif line.startswith("# Assembly level:"):
                        cur_dict["assembly_level"] = line.split(":")[-1].lstrip().strip("\n")
                        cur_dict["assembly_level_int"] = ASSEMBLY_LEVEL_DIC[cur_dict["assembly_level"].lower()]
                    elif line.startswith("# RefSeq assembly accession:"):
                        cur_dict["refseq_accesion_id"] = line.split()[-1].strip("\n")
                        break
            if not cur_dict.get("organism_name") or not cur_dict.get("date") or not cur_dict.get(
                    "refseq_accesion_id") or not cur_dict.get("assembly_level"):
                print("ERROR: file= %s has missing info" % stats_filename)
                break

            org_full_name = get_full_name(cur_dict)

            if not org_dict.get(org_full_name):
                org_dict[org_full_name] = []
            org_dict[org_full_name].append(
                {"date": cur_dict["date"], "refseq_accesion_id": cur_dict["refseq_accesion_id"],
                 "assembly_level": cur_dict["assembly_level"], "assembly_level_int": cur_dict["assembly_level_int"]})
    return org_dict


def get_full_name(cur_dict):
    org_name = re.sub(r"\([^()]*\)", "", cur_dict["organism_name"]).strip()  # remove whatever is in parenthesis
    if (not cur_dict.get("strain")) or (cur_dict["strain"] in org_name):
        return org_name
    else:
        return org_name + " " + cur_dict["strain"]


def get_best_assemblies_per_org_df(stats_folder):
    organism_dict = get_assemblies_per_org_dict(stats_folder)
    df_dict = {"organism name": [], "number of assemblies": [], "last assembly level": [], "if best assembly level": [],
               "last RefSeq accession ID": [], "last assembly date": []}
    for org, items in organism_dict.items():
        first_item = items[0]  # default value
        last_date = datetime.strptime(first_item["date"], "%Y-%m-%d")
        last_accession_id = first_item["refseq_accesion_id"]
        last_assembly_level = first_item["assembly_level"]
        last_assembly_level_int = first_item["assembly_level_int"]
        min_assembly_level_int = first_item["assembly_level_int"]

        if len(items) > 1:  # more than one assembly of the genome
            for item in items:
                date = datetime.strptime(item["date"], "%Y-%m-%d")
                if date > last_date:
                    last_date = date
                    last_accession_id = item["refseq_accesion_id"]
                    if item["assembly_level_int"] <= min_assembly_level_int:
                        min_assembly_level_int = item["assembly_level_int"]
                    else:
                        print(item["assembly_level"], last_assembly_level)
                    last_assembly_level_int = item["assembly_level_int"]
                    last_assembly_level = item["assembly_level"]

        df_dict["organism name"].append(org)
        df_dict["last assembly date"].append(last_date)
        df_dict["number of assemblies"].append(len(items))
        df_dict["last RefSeq accession ID"].append(last_accession_id)
        df_dict["last assembly level"].append(last_assembly_level)
        df_dict["if best assembly level"].append(last_assembly_level_int == min_assembly_level_int)
    return pd.DataFrame(df_dict)


# TODO: turn this into general function for downloading stats or assemblies
def download_fastas(ftp, bacteria_folders_clean, save_folder, accession_ids_list, bacteria_input):
    start = time.time()
    print(datetime.now())

    num_already_downloaded = 0
    num_errors_in_try = 0

    print(f"Starting to download assemblies of bacteria {bacteria_input}")
    for i, bacteria in enumerate(bacteria_folders_clean):
        file = None
        all_assembly_folder = bacteria + "/" + ALL_ASSEMBLY_PREFIX  # find out if there is a way of downloading this whole folder in one go (even better - of downloading only the files in this folder that end with _genomic.fna.gz)
        assemblies = ftp.nlst(all_assembly_folder)
        for j in range(len(assemblies)):
            file_basename = os.path.basename(assemblies[j])
            inlist = False
            for acc in accession_ids_list:
                if file_basename.startswith(acc):
                    inlist = True
                    try:
                        file2download = f'{assemblies[j]}/{os.path.split(assemblies[j])[-1]}_genomic.fna.gz'
                        local_filename = os.path.join(save_folder, os.path.basename(os.path.normpath(file2download)))
                        file = open(local_filename, 'wb')
                        ftp.retrbinary('RETR ' + file2download, file.write)
                        # print(f'finished downloading file {local_filename}')
                    except:
                        print(f'error downloading fasta file for {assemblies[j]} of bacteria {bacteria}')
                        num_errors_in_try += 1
                    file.close()
                    num_already_downloaded += 1
                    if num_already_downloaded % 100 == 0:
                        end = time.time()
                        print(f"Runtime of {num_already_downloaded} downloads is {(end - start) / 60} minutes")
            if not inlist:
                print("filename %s was not in list" % file_basename)
        if num_errors_in_try > NUM_OF_ERRORS:  # if there are more than NUM_OF_ERRORS unsuccessful downloads then try running all over again
            break
    print(f"Runtime of all downloads is {(time.time() - start) / 60} minutes")
    return num_already_downloaded


def download_stats(ftp, bacteria_folders_clean, save_folder, bacteria_input):
    start = time.time()
    print(datetime.now())

    num_already_downloaded = 0
    num_errors_in_try = 0

    print(f"Starting to download assemblies stats of bacteria {bacteria_input}")

    for i, bacteria in enumerate(bacteria_folders_clean):
        file = None
        all_assembly_folder = bacteria + "/" + ALL_ASSEMBLY_PREFIX  # find out if there is a way of downloading this whole folder in one go (even better - of downloading only the files in this folder that end with _genomic.fna.gz)
        assemblies = ftp.nlst(all_assembly_folder)
        for j in range(len(assemblies)):
            try:
                file2download = f'{assemblies[j]}/{os.path.split(assemblies[j])[-1]}_assembly_stats.txt'
                local_filename = os.path.join(save_folder, os.path.basename(os.path.normpath(file2download)))
                file = open(local_filename, 'wb')
                ftp.retrbinary('RETR ' + file2download, file.write)
                # print(f'finished downloading file {local_filename}')
            except:
                print(f'error downloading stats file for {assemblies[j]} of bacteria {bacteria}')
                num_errors_in_try += 1
            file.close()
            num_already_downloaded += 1
            if num_already_downloaded % 100 == 0:
                end = time.time()
                print(f"Runtime of {num_already_downloaded} downloads is {(end - start) / 60} minutes")

        if num_errors_in_try > NUM_OF_ERRORS:  # if there are more than NUM_OF_ERRORS unsuccessful downloads then try running all over again
            break
    print(f"Runtime of all downloads is {(time.time() - start) / 60} minutes")
    return num_already_downloaded


def get_genomes(num_try: int, save_folder, bacteria_input):
    print(f'starting num_try {num_try}')

    ftp = FTP(SERVER_PATH)
    ftp.login()
    print('logged in')

    ftp.cwd(SERVER_BACTERIA_LOCATION)
    print(f'changed dir to {SERVER_BACTERIA_LOCATION}')

    all_ncbi_bacterias = ftp.nlst(SERVER_BACTERIA_LOCATION)  # get folders within the directory
    bacterias_of_interest = [filename for filename in all_ncbi_bacterias if
                              bacteria_input in os.path.basename(filename)]  # get folders of searched bacteria

    local_bacteria_folder = os.path.join(os.path.normpath(save_folder),
                                         os.path.basename(os.path.normpath(bacteria_input)))
    local_stats_folder = os.path.join(local_bacteria_folder, os.path.normpath("stats"))
    local_fasta_folder = os.path.join(local_bacteria_folder, os.path.normpath("fasta"))

    # see if there is a way to write this in one line
    for dir_path in [local_bacteria_folder, local_stats_folder, local_fasta_folder]:
        os.makedirs(dir_path, exist_ok=True)
        print("created dir %s" % dir_path)

    num_stats_downloaded = download_stats(ftp, bacterias_of_interest, local_stats_folder, bacteria_input)
    print(f'num_try {num_try}, already_downloaded {num_stats_downloaded} | len {len(bacterias_of_interest)}')

    accession_ids_to_download = list(get_best_assemblies_per_org_df(local_stats_folder)["last RefSeq accession ID"])
    num_fasta_downloaded = download_fastas(ftp, bacterias_of_interest, local_fasta_folder, accession_ids_to_download, bacteria_input)
    print(f'num_try {num_try}, already_downloaded {num_fasta_downloaded} | len {len(bacterias_of_interest)}')

    ftp.quit()  # This is the “polite” way to close a connection
    print("Successfully quited ftp connection")


def run(save_folder, bacteria_input):
    for i in range(NUMBER_OF_TRIES):
        try:
            get_genomes(i, save_folder, bacteria_input)
        except Exception as e:
            print(e)
            print(f'finished try number: {i} out of {NUMBER_OF_TRIES}')
        else:
            print("Exiting")
            break


def main():
    args = sys.argv[1:]
    run(args[0], args[1])


if __name__ == "__main__":
    main()
