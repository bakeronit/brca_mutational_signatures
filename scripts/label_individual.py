import sys, os
from cyvcf2 import VCF
from collections import defaultdict

cristina_clf = "germline_mutation_classification.final.txt"

brca_mutation = defaultdict(list)

def add_brca12(cohort, mutations):
    brca1 = ""
    brca2 = ""
    for mutation in mutations:
        if 'BRCA1' in mutation:
            brca1 = "_brca1"
        elif 'BRCA2' in mutation:
            brca2 = "_brca2"
    return f"{cohort}{brca1}{brca2}"


def get_cohort_type(mutation_group):
    if len(mutation_group['A']) > 0:
        cohort = "exclude"
    elif len(mutation_group['D']) > 0:
        cohort = add_brca12("positive", mutation_group['D'])
    elif len(mutation_group['E']) > 0:
        cohort = add_brca12("additional", mutation_group['E'])
    elif len(mutation_group['B']) > 0:
        cohort = "exclude"
    else:
        cohort = "negative_brca1/2"
    return cohort

def get_mutation_AD(mutation, vcf):
    chrom, pos, ref, alt, gene = mutation.split("_")
    ref_alt_ad = ""
    for v in vcf(f"{chrom}:{pos}-{pos}"):
        if ref != v.REF or not alt in v.ALT:
            print(f'WARNING: {mutation} not match')
            sys.exit()
        ref_alt_ad = f"{v.gt_ref_depths[0]}:{v.gt_alt_depths[0]}|{v.gt_ref_depths[1]}:{v.gt_alt_depths[1]}"  # in the vcf file, the first sample is normal blood and the second one is primary tumour
    
    return ref_alt_ad

with open(cristina_clf,'r') as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        line = line.strip().split()
        mutation = "_".join(line[:5])
        group = line[-1]
        brca_mutation[group].append(mutation)


filename = sys.argv[1]
sample_id = os.path.basename(filename).replace('.final.tsv','')

groups = ['A','B','C','D','E','F']
mutation_group = {group:[] for group in groups}
hgvs_group = {group:[] for group in groups}

with open(filename, 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        mutation = "_".join(line[:4] + [line[8]])
        hgvs = line[18]
        gene = line[8]
        for group in groups:
            if mutation in brca_mutation[group]:
                if group in ["A", "D"]:
                    #vcf_file = filename.replace('../3.collect_snvs/germline_snvs','vcf_files').replace('.final.tsv','.vcf.gz')
                    vcf_file = filename.replace('germline_snvs','../4.brca_classification/vcf_files').replace('.final.tsv','.vcf.gz')
                    vcf = VCF(vcf_file)
                    ref_alt_ad = get_mutation_AD(mutation, vcf)
                    hgvs += f"[{ref_alt_ad}]"  ## for a BRCA1/2 path mutation, we want to know their AD in ref and alt, so can infer more in combine with copynumber
                mutation_group[group].append(mutation)
                hgvs_group[group].append(f"{hgvs}({gene})")

cohort_type = get_cohort_type(mutation_group)

print(f"{sample_id}\t{"\t".join(str(len(value)) for value in mutation_group.values())}\t{cohort_type}\t{";".join(hgvs_group["D"])}\t{";".join(hgvs_group["E"])}\t{";".join(hgvs_group["A"])}\t{";".join(hgvs_group["B"])}")
