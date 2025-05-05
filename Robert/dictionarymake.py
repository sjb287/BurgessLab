from Bio import SeqIO
import matplotlib.pyplot as plt
beforeBarcode = "GCAATACCGTCTCACTGAACTGGCCGATAATTGCAGACGGGATCCGTATC"
afterBarcode = "AACGGTGTTGTTATCCGATACAACCGGATATTTTTCTTTTAATGAGTCTA"
beforeMutant = "TTTCTCAAAATAAATCGATACTGCATTTCTAGGCATATCCAGCGAGATCT"
afterMutant = "TAACACAAAGACACCGACAACTTTCTTGTATCGGATCCATCCTAACTCGA"

middleBackbone = ""
earlyBackbone = ""
lateBackbone = ""
def main():
    try:
        #change when fastq data is available
        with open('/Users/robert/Developer/research/dictionary/UTF-8QYT8X7_1_linear_pSBPevol_2band.fastq', 'r') as handle:
            records = SeqIO.parse(handle, "fastq")
            goodreads = (record
                 for record in records
                 if min(record.letter_annotations["phred_quality"]) >= 20)
            lookupTable = {} 
            frequencyTable = {}
            for record in records:  
                try:
                    seq,barcode = findSequenceAndBarcode(record.seq)
                    
                    #checks if backbone between mutant and barcode is not mutated
                    if middleBackbone != record.seq[record.seq.find(barcode) + len(barcode):record.seq.find(seq)]:
                        break
                    #check up until the barcode
                    if earlyBackbone != record.seq[0:record.seq.find(barcode)]:
                        break
                    #check after the mutant
                    if lateBackbone != record.seq[record.seq.find(seq)]:
                        break

                    #by now sequence should not have any mutantions except for the variable region.
                    if seq not in lookupTable:
                        lookupTable[barcode] = seq
                        frequencyTable[barcode] = 1
                    else:
                        frequencyTable[barcode] += 1
                except TypeError:
                    print("couldnt find barcode in reverse or forward")
                except AttributeError:
                    print("mutation in backbone")
            sequences = list(lookupTable.values())
            values = list(frequencyTable.values())
            plt.hist(values)
            plt.show()
    except FileNotFoundError:
        print("file not found")
        return

def findSequenceAndBarcode(read):
    try:
        read.index(beforeBarcode)
        read.index(afterBarcode)
        read.index(beforeMutant)
        read.index(afterMutant)
    except ValueError:
        return findReverse(read)
    barcode = read[read.find(beforeBarcode) + len(beforeBarcode):read.find(afterBarcode)]
    sequence = read[read.find(beforeMutant) + len(beforeMutant):read.find(afterMutant)]
    if(barcode == "" or sequence == ""):
        return findReverse(read)
    return sequence,barcode
    
def findReverse(read):
    try:
        read.index(beforeBarcode)
        read.index(afterBarcode)
        read.index(beforeMutant)
        read.index(afterMutant)
    except ValueError:
        return TypeError
    barcode = read[read.find(beforeBarcode) + len(beforeBarcode):read.find(afterBarcode)]
    sequence = read[read.find(beforeMutant) + len(beforeMutant):read.find(afterMutant)]
    if(barcode == "" or sequence == ""):
        return TypeError
    return sequence,barcode

def reverseCompiment(read):
    newString = ""
    for char in read:
        if char == 'A':
            newString += 'T'
        if char == 'T':
            newString += 'A'
        if char == 'C':
            newString += 'G'
        if char == 'G':
            newString += 'C'
    return newString[::-1] #return the reversed string
if __name__ == "__main__":
    main()
