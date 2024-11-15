## read `Samples.txt` file from results and sorted as FFPEsig.py required.
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input", help="Input file path")
parser.add_argument("output", help="Output file path")
parser.add_argument("--reverse", action="store_true", help="Reverse the input content")
args = parser.parse_args()

channel6 = ['C>A','C>G','C>T','T>A','T>C','T>G']
channel96 = ['ACA', 'ACC', 'ACG', 'ACT', 'CCA', 'CCC', 'CCG', 'CCT', 'GCA',
            'GCC', 'GCG', 'GCT', 'TCA', 'TCC', 'TCG', 'TCT', 'ACA', 'ACC',
            'ACG', 'ACT', 'CCA', 'CCC', 'CCG', 'CCT', 'GCA', 'GCC', 'GCG',
            'GCT', 'TCA', 'TCC', 'TCG', 'TCT', 'ACA', 'ACC', 'ACG', 'ACT',
            'CCA', 'CCC', 'CCG', 'CCT', 'GCA', 'GCC', 'GCG', 'GCT', 'TCA',
            'TCC', 'TCG', 'TCT', 'ATA', 'ATC', 'ATG', 'ATT', 'CTA', 'CTC',
            'CTG', 'CTT', 'GTA', 'GTC', 'GTG', 'GTT', 'TTA', 'TTC', 'TTG',
            'TTT', 'ATA', 'ATC', 'ATG', 'ATT', 'CTA', 'CTC', 'CTG', 'CTT',
            'GTA', 'GTC', 'GTG', 'GTT', 'TTA', 'TTC', 'TTG', 'TTT', 'ATA',
            'ATC', 'ATG', 'ATT', 'CTA', 'CTC', 'CTG', 'CTT', 'GTA', 'GTC',
            'GTG', 'GTT', 'TTA', 'TTC', 'TTG', 'TTT']

formatted_channel = []

for i,context in enumerate(channel96):
    index_m = int((i)/16)
    mutation = channel6[index_m]
    c = context[:1] + f"[{mutation}]" + context[2:]
    formatted_channel.append(c)

df = pd.read_csv(args.input,sep="\t")

if args.reverse:
    df_og = pd.read_csv("results/GRCh38/MAGIC/COSMIC2/SBS96/Samples.txt",sep="\t")  # re-match the row order of sigprofilerExtractor.
    df2 = pd.read_csv("data/GRCh38/MAGIC.sorted.matrix.txt",sep="\t") # the dataframe with mutation type column
    df.insert(0,'MutationType', df2['MutationType'])
    sorted_df = df.set_index('MutationType').reindex(df_og['MutationType']).reset_index()
else:
    df['MutationType']=pd.Categorical(df['MutationType'],categories=formatted_channel, ordered=True) ## sorted the row to be able to run FFPEsig.py
    sorted_df = df.sort_values(by='MutationType')

sorted_df.to_csv(args.output, index=False, sep="\t")


