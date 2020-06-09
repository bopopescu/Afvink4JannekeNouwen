from flask import Flask, render_template, request
import re
import math
from Bio.Blast import NCBIXML

app = Flask(__name__)


@app.route('/')
def hello_world():
    dna = rna = protein = False
    result = ""
    seq = ""
    blast_result = []
    dna_result = {}
    if request.args.get("seq"):
        seq = request.args.get("seq").upper()
        if re.search('^[ATCG]+$', seq):
            dna = True
            result = "The given sequence is a DNA sequence"
            dna_result['dna'] = seq
            dna_result['rna'] = seq.replace("T", "U")
            dna_result['protein'] = dna_to_protein(seq[:math.floor(len(seq)/3)*3])
        if re.search('^[AUCG]+$', seq):
            rna = True
            if dna:
                result = "The given sequence could be a DNA sequence or an RNA" \
                         " sequence."
            else:
                result = "The given sequence is an RNA sequence."
        if re.search('^[GPAVLIMCFYWHKRQNEDST\\*]+$', seq):
            if not dna and not rna:
                protein = True
                result = "The given sequence is a protein sequence."
                blast(seq)
                blast_result = get_blast_info()
        if not dna and not rna and not protein:
            result = "The given sequence is neiter dna, nor rna nor a protein."
    return render_template("index.html", result=result, protein=protein,
                           seq=seq, blast_result=blast_result,
                           dna_result=dna_result)


def blast(seq):
    blast_result = NCBIWWW.qblast("tblastn", "nr", seq,
                                  format_type="XML", hitlist_size=1)
    blast_data_file = open("html_blast.xml", "w")
    blast_data_file.write(blast_result.read())
    blast_data_file.close()


def get_blast_info():
    with open("html_blast.xml", "r") as out_handle:
        blast_records = NCBIXML.parse(out_handle)
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    #  De titel wordt opgesplitst in acessiecode en desciptie
                    accession = alignment.title.split("|")[3]
                    description = alignment.title.split("|")[4].lstrip()
                    e_value = hsp.expect
                    perc_identity = round(
                        hsp.identities / hsp.align_length * 100, 3)
                    total_score = hsp.score
    return description, accession, e_value, perc_identity, total_score


def dna_to_protein(seq):
    """zet coding stand van DNA om naar aminozuren"""
    dna_codon_to_protein = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X',
        'TGC': 'C', 'TGT': 'C', 'TGA': 'X', 'TGG': 'W'
    }
    protein = []
    for i in range(0, len(seq), 3):
        protein.append(dna_codon_to_protein[seq[i:i+3]])
    return "".join(protein)


if __name__ == '__main__':
    app.run()

# Ontwikkel een webapplicatie die in staat is om van een gegeven sequentie te bepalen of het DNA, RNA, eiwit of geen van drieen is.
# Wanneer de applicatie in staat is om te detecteren dat het om DNA gaat, geeft het de corresponderende RNA en eiwit sequenties.
# Als het een eiwit is geeft het meest waaarschijnlijke gen waar het van afkomstig is.

# vamaeeegsrfpyvfwgsknifglanpddvrnicdakgnsfvdamkacgftlpnapltpr
