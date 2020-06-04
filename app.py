from flask import Flask, render_template, request
import re
from Bio.Blast import NCBIWWW, NCBIXML


app = Flask(__name__)


@app.route('/')
def hello_world():
    dna = rna = protein = False
    result = ""
    seq = ""
    if request.args.get("seq"):
        seq = request.args.get("seq").upper()
        if re.search('^[ATCG]+$', seq):
            dna = True
            result = "The given sequence is a DNA sequence"
        if re.search('^[AUCG]+$', seq):
            rna = True
            if dna:
                result = "The given sequence could be a DNA sequence or an RNA sequence"
            else:
                result = "The given sequence is an RNA sequence"
        if re.search('^[GPAVLIMCFYWHKRQNEDST\\*]+$', seq):
            protein = True
            if dna and rna:
                result = "The given sequence could be a DNA sequence, an RNA sequence, or a protein sequence"
            elif dna:
                result = "The given sequence could be a DNA sequence or a protein sequence"
            elif rna:
                result = "The given sequence could be an RNA sequence or a protein sequence"
            else:
                result = "The given sequence is a protein sequence"
                blast(seq)
                get_blast_info()
        if not dna and not rna and not protein:
            result = "The given sequence is neiter dna, nor rna nor a protein."
    return render_template("index.html", result=result, protein=protein, seq=seq)


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
                    #  De titel wordt opgesplitst in acessiecode, desciptie
                    #  en organisme
                    accession = alignment.title.split("|")[1]
                    description = alignment.title.split("|")[2]
                    organism = (alignment.title.split("[")[1]).split("]")[0]

                    e_value = hsp.expect
                    perc_identity = hsp.identities / hsp.align_length * 100
                    total_score = hsp.score
                    print(accession, description, organism, e_value, perc_identity, total_score)
    return accession, description, organism, e_value, perc_identity, total_score



if __name__ == '__main__':
    app.run()


# Ontwikkel een webapplicatie die in staat is om van een gegeven sequentie te bepalen of het DNA, RNA, eiwit of geen van drieen is.
# Wanneer de applicatie in staat is om te detecteren dat het om DNA gaat, geeft het de corresponderende RNA en eiwit sequenties.
# Als het een eiwit is geeft het meest waaarschijnlijke gen waar het van afkomstig is.

# vamaeeegsrfpyvfwgsknifglanpddvrnicdakgnsfvdamkacgftlpnapltpr