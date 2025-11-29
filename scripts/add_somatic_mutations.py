import sys
from pathlib import Path

annotation_file = sys.argv[1]

def read_maf(file):
    somatic_mutations = []
    with open(file, 'rt') as fh:
        for line in fh:
            if line.startswith("Hugo"):
                continue
            line = line.strip().split("\t")
            gene = line[0]
            chrom, start, end = line[4:7]
            effect = line[8]
            if not effect in ["Intron","5'Flank","3'Flank","3'UTR","5'UTR"]:
                hgvs = line[52]
                if effect == "Splice_Site":
                    hgvs = f"{chrom}:{start}-{end}" #add mutation id for splice site effect. it's null.
                somatic_mutations.append(f"{gene}.{hgvs}.{effect}")

    return somatic_mutations


def read_cnv(file):
    ascat_call = dict()
    cnv = []
    with open(file,'rt') as fh:
        for line in fh:
            line = line.strip().split('\t')
            if not ascat_call.get(line[0]) or ascat_call[line[0]] == "NA:NA": # if multiple calls, keep the first one except NA
                ascat_call[line[0]] = f"{line[2]}:{line[1]}" # add absolute copynumber 
    for gene in gene_order:
        cnv.append(ascat_call[gene])

    return "\t".join(cnv)

gene_order = ["BRCA1","BRCA2","ATM","BARD1","BRIP1","CHEK2","PALB2","PTEN","RAD51C","RAD51D","TP53"]
with open(annotation_file, 'rt') as fh:
    for line in fh:
        line = line.strip()
        if line.startswith("cohort"):
            new_header = f"{line}\tsomatic_mutations\t{'\t'.join(gene_order)}"
            print(new_header)
            header = line.split("\t")
            continue
        cols={ i:j for i,j in zip(header,line.split('\t'))}
        cohort = cols['cohort']
        tumour_id = cols['tumour_id']
        donor_id = cols['sample_id']
        cohort = "TCGA-BRCA" if cohort == "TCGA" else cohort  # I named TCGA as TCGA-BCRA at first, so the folder name is different from cohort name.
        cohort = "Nones2019" if cohort == "Familial Breast" else cohort
        cohort = "Qimprove" if cohort == "Q-IMPROvE" else cohort
        if cohort in ["MAGIC","Qimprove"]:
            maf_file = f"../3.collect_snvs/somatic_snvs/{cohort}/{donor_id}.{tumour_id}.maf"
            cnv_file = f"ascat_annotation/{cohort}/{donor_id}.{tumour_id}.cn_genelist.txt"
        else:
            maf_file = f"../3.collect_snvs/somatic_snvs/{cohort}/{donor_id}.maf"
            cnv_file = f"ascat_annotation/{cohort}/{donor_id}.cn_genelist.txt"
        maf_file = Path(maf_file)
        cnv_file = Path(cnv_file)

        print(line,end="\t")
        if maf_file.is_file() and cnv_file.is_file():
            somatic_mutations = read_maf(maf_file)
            print(";".join(somatic_mutations), end="\t")
            cnv_str = read_cnv(cnv_file)
            print(cnv_str)
        else:
            print("NA\tNA")
