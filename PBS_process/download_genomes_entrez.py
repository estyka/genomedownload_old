from Bio import Entrez
import os
import urllib

"""
TODO:
- test more to see that it downloads like assembly ncbi
- downloading only latest_refseq
- test Rickettsiales
"""

def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record

def get_assemblies(term, download=True, path='assemblies'):
    """Download genbank assemblies for a given search term.
    Args:
        term: search term, usually organism name
        download: whether to download the results
        path: folder to save to
    """

    #provide your own mail here
    Entrez.email = "estykatzeff@mail.tau.ac.il"
    handle = Entrez.esearch(db="assembly", term=term, retmax='200')
    record = Entrez.read(handle)
    ids = record['IdList']
    print (f'found {len(ids)} ids')
    latest_refseq_assemblies_count = 0
    links = []
    for id in ids:
        summary = get_assembly_summary(id)
        is_refseq = 'latest_refseq' in summary['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']
        latest_refseq_assemblies_count += int(is_refseq)
        print(is_refseq)
        url_genbank = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']  #get ftp url
        url_stats_rpt = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_Stats_rpt']  #get ftp url
        url_assembly_rpt = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_Assembly_rpt']  #get ftp url
        # 'FtpPath_GenBank', 'FtpPath_RefSeq', 'FtpPath_Assembly_rpt', 'FtpPath_Stats_rpt'
        if url_genbank == '':
            continue
        label = os.path.basename(url_genbank)
        link = url_genbank + "/" + label + '_genomic.fna.gz'  #get the fasta link - change this to get other formats

        print (link)
        links.append(link)
        dir = r"C:\Users\97252\Documents\year_4\bio_project_data\download_results\results"
        if download == True:
            urllib.request.urlretrieve(link, os.path.join(dir, f'{label}.fna.gz'))  #download and save link
            # TODO: add a verification here that the download happened (check that file in location isfile()).
            #urllib.request.urlretrieve(url_stats_rpt)
    print(f'latest_refseq_assemblies_count: {latest_refseq_assemblies_count}')
    return links

bacteria2search = "Xanthomonas campestris"
links = get_assemblies(bacteria2search, download=True)
