beforeBarcode = "GCAATACCGTCTCACTGAACTGGCCGATAATTGCAGACGGGATCCGTATC"
afterBarcode = "AACGGTGTTGTTATCCGATACAACCGGATATTTTTCTTTTAATGAGTCTA"
beforeMutant = "TTTCTCAAAATAAATCGATACTGCATTTCTAGGCATATCCAGCGAGATCT"
afterMutant = "TAACACAAAGACACCGACAACTTTCTTGTATCGGATCCATCCTAACTCGA"

earlyBackbone = ""
middleBackbone = ""
lateBackBone = ""
from Bio import SeqIO

def main():
    try:
        #change when fastq data is available
        with open('/Users/robert/Developer/research/dictionary/UTF-8QYT8X7_1_linear_pSBPevol_2band.fastq', 'r') as handle:
            records = list(SeqIO.parse(handle, "fastq"))
            for record1 in records:
                print("hi")
                for record2 in records:
                    if record1.seq != record2.seq:
                        backboneFinder(record2.seq,record1.seq)
        print(earlyBackbone)
        print(middleBackbone)
        print(lateBackBone)
    except FileNotFoundError:
        return("file not found")

def findSequenceAndBarcode(read):
    try:
        read.index(beforeBarcode)
        read.index(afterBarcode)
        read.index(beforeMutant)
        read.index(afterMutant)
    except ValueError:
        return ValueError
    barcode = read[read.find(beforeBarcode) + len(beforeBarcode):read.find(afterBarcode)]
    sequence = read[read.find(beforeMutant) + len(beforeMutant):read.find(afterMutant)]
    if(barcode == "" or sequence == ""):
        return ValueError
    return sequence, barcode

def backboneFinder(sequence1,sequence2):
    try:
        barcode1, mutant1 = findSequenceAndBarcode(sequence1)
        barcode2, mutant2 = findSequenceAndBarcode(sequence2)

        if sequence1[0:sequence1.find(barcode1)] == sequence2[0:sequence2.find(barcode2)]:
            if earlyBackbone != sequence1[0:barcode1]:
                print("new earlybackbone")
                earlyBackbone = sequence1[0:barcode1]

        if sequence1[sequence1.find(barcode1) + len(barcode1) : sequence1.find(mutant1)] == sequence2[sequence2.find(barcode2) + len(barcode2) : sequence1.find(mutant2)]:
            if middleBackbone != sequence1[sequence1.find(barcode1) + len(barcode1) : sequence1.find(mutant1)]:
                print("new middlebackbone")
                middleBackbone = sequence1[sequence1.find(barcode1) + len(barcode1) : sequence1.find(mutant1)]
        
        if sequence1[sequence1.find(mutant1) + len(mutant1):] == sequence2[sequence2.find(mutant2) + len(mutant2):]:
            if lateBackBone != sequence1[sequence1.find(mutant1) + len(mutant1):]:
                print("found new late backbone")
                lateBackBone = sequence1[sequence1.find(mutant1) + len(mutant1):]
    except:
        return


if __name__ == "__main__":
    main()
