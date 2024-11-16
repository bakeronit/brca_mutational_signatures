import subprocess
import sys
import pandas as pd
from collections import defaultdict

cohort_input = {
    "familial": {
        "vcf_file": "vcf_files/all_Simpson_FamilialBreast_GRCh38.germline.August2024_final.vcf",
        "report": "../../3.collect_snvs/Simpson_-_FamilialBreast_Germline_Report_with_Nannotation_May2024.csv",
        "cf_annot": "Familial.tsv",
        "donor_col": "donorPublicationID"
    },
    "qimprove": {
        "vcf_file": "vcf_files/all_qimprove_genomiqa_GRCh38.germline.July2024_final.vcf",
        "report": "../../3.collect_snvs/GenomiQa_-_qimprove_GRCh38_Germline_Report_July2024_WITH_NANNO.tsv",
        "cf_annot": "Qimprove.tsv",
        "donor_col": "patientLabel"
    },
    "tcga": {
        "vcf_file": "vcf_files/all_TCGA_Breast_cancer_GRCh38.germline.August2024_final.vcf",
        "report": "../../3.collect_snvs/TCGA_-_breast_cancer_Germline_Report_with_nannotation.June17updated.csv",
        "cf_annot": "TCGA.tsv",
        "donor_col": "donorLabel",
        "uuid_mapping": "../../0.data_resource/GRCh38/uuid_case.txt"
    },
    "magic": {
        "vcf_file": "vcf_files/all_Waddell_-_Magic_WGS_GRCh38.germline.July2024_final.vcf",
        "report": "../../3.collect_snvs/Waddell_-_Magic_WGS_GRCh38_Germline_Report_with_Nanno_Annotation_May2024.csv",
        "cf_annot": "MAGIC.tsv",
        "donor_col": "donorLabel"
    }
}

all_sample_id = pd.read_csv("../brca_germline_clf_tumor_hrd.v5.tsv", sep="\t")["donor_id"].values

cohort = sys.argv[1].lower()
if not cohort in cohort_input:
    raise ValueError(f"Invalid cohort: {cohort}, only {', '.join(cohort_input.keys())} are allowed")
if cohort == "tcga":
    tcga_mapping = {line.strip().split('\t')[0]: line.strip().split('\t')[1] for line in open(cohort_input[cohort]["uuid_mapping"], 'rt').readlines()}

vcf_file = cohort_input[cohort]["vcf_file"]
report_file = cohort_input[cohort]["report"]

if report_file.endswith("tsv"):
    report = pd.read_csv(report_file, sep="\t")
else:
    report = pd.read_csv(report_file)

cf_annot = cohort_input[cohort]["cf_annot"]

sample_list = subprocess.check_output(f"bcftools query -l {vcf_file}", shell=True).decode().split("\n")[:-1]

def get_sample_gt(sample_list, sv_genotypes):
    donor_gt = {}
    for idx, uuid in enumerate(sample_list):
        if not (donor_id:= report[ report["sampleIRI"].str.contains(uuid) ]).empty:
            donor_id = donor_id[cohort_input[cohort]["donor_col"]].values[0]
            if cohort == "tcga":
                donor_id = tcga_mapping[donor_id]
            if pd.isna(donor_id) or not donor_id in all_sample_id:
                continue
            gt = sv_genotypes[idx]
            donor_gt[donor_id] = gt
    return donor_gt

sv = defaultdict(list)
sv_info = defaultdict(str)
with open(cf_annot, "r") as fh:
    for line in fh:
        line = line.strip()
        sv_id, sv_type, sv_len, bp1_gene, bp2_gene = line.split("\t")[4:9]
        sv_class = line.split("\t")[-1]
        sv_string = f"{sv_class}_{sv_id}_{sv_type}_{sv_len}_{bp1_gene}_{bp2_gene}"
        sv[sv_class].append(sv_id)
        sv_info[sv_id] = sv_string
        if sv_type == "BND":
            sv["BND"].append(sv_id)
            sv_string = f"{sv_id}_{sv_type}_{sv_len}_{bp1_gene}_{bp2_gene}"
            sv_info[sv_id] = sv_string


querys = subprocess.check_output(f"bcftools query -f '%CHROM\t%POS\t%INFO/CHR2\t%INFO/END\t%ID\t%INFO/SVTYPE\t%INFO/INSLEN\t[%GT ]' {vcf_file}", shell=True).decode().split("\n")[:-1]
for query in querys:
    sv_id = query.split("\t")[4]
    if sv_id in sv["A"] or sv_id in sv["D"] or sv_id in sv["E"] or sv_id in sv["BND"]:
        genotypes = query.split("\t")[-1].split(" ")
        donor_gt = get_sample_gt(sample_list, genotypes)
        for donor_id in donor_gt:
            if "1" in donor_gt[donor_id]:
                print(f"{donor_id}\t{sv_info[sv_id]}") 
