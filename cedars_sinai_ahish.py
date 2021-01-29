#!/usr/bin/env python3

import argparse, re, csv

#STEP 1:
#Creating dictionary where key:header and value:protein sequence
def retrieveFasta(fasta_file):
    with open(fasta_file) as fh1:
        header_seqs = {}

        for line in fh1:
            line = line.strip() #Removing "\n" if any

            if line[0] == ">":
                words = line.split()
                header = words[0][1:]   #Removing ">" character
                header_seqs[header] = ''

            else:
                header_seqs[header] = header_seqs[header] + line

    return header_seqs

#STEP 2:
#Getting the modified peptide and the position of the modification in the peptide:
def getPTM(peptide_file):
    with open(peptide_file) as fh:
        result_dict = {}

        for seq in fh:
            i = 0
            index = 0
            seq = seq.strip("\n")

            for letter in seq:

                if letter == "(":
                    start_point_site = i   #start_point_site becomes 0
                    start_point = index    #start_point becomes 0

                if letter == ")":
                    end_point = index      #end_point becomes 0

                    if seq in result_dict:  #Next iteration, as the seq is already added, appending the different start_point_site to the list
                        result_dict[seq].append(start_point_site)

                    else:   #Adds respective seq to the dict with the start_point_site of respective PTM
                        result_dict[seq] = [start_point_site]

                    i = i - (end_point - start_point) - 1   #Resetting i for the next iteration

                i += 1
                index += 1

        #Removing duplicates. Could have taken a set too, but the order in a set is undefined hence appending to a list and removing dups is easier:
        for key in result_dict:
            a = 0

            while a < len(result_dict[key]):
              b = a + 1

              while b < len(result_dict[key]):

                 if result_dict[key][a] == result_dict[key][b]:
                    del result_dict[key][b]

                 else:
                    b += 1
              a += 1

    return result_dict

#STEP 3:
#Printing out final result
def getPos(header_seqs, result_dict):
    test_dict = {}  #Use test_dict in order to store positions and update result_dict with positions of PTMs with respect to protein sequence
    header_list = []
    sequence_list = []
    position_list = []
    for key1,value1 in result_dict.items():

        for key2,value2 in header_seqs.items():
            match = re.search(re.sub("\(.*?\)", "", key1), value2)  #Searching for sequence in fasta dict

            if match:
                match_start_position = match.start()
                test_dict[key1] = match_start_position
                result_dict[key1] = [x + test_dict[key1] for x in result_dict[key1]]    #Updating result_dict positions

                header_list.append(key2)
                sequence_list.append(re.sub("\(.*?\)", "", key1))
                position_list.append(";".join(repr(e) for e in result_dict[key1]))

    final_result = zip(header_list, sequence_list, position_list)

    return final_result

#STEP 4:
#Writing out results to csv:
def writeCSV(final_result, output_file):
    with open(output_file+".csv", "a") as csvfile:
        writer = csv.writer(csvfile, delimiter=",", )
        writer.writerow(["ProteinName", "FullPeptideName", "ModifiedAminoAcidPosition"])
        writer.writerows(final_result)

def main():

    #Argparse code:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required = True, help = "Enter complete path of input file containing peptides.")
    parser.add_argument("-f", required = True, help = "Enter complete path of protein FASTA file.")
    parser.add_argument("-o", required = True, help = "Enter your output file name.")
    args = parser.parse_args()

    #Populating variables:
    peptide_file = args.i
    fasta_file = args.f
    output_file = args.o

    #Calling functions:
    header_seqs = retrieveFasta(fasta_file)
    result_dict = getPTM(peptide_file)
    final_result = getPos(header_seqs, result_dict)
    writeCSV(final_result, output_file)

if __name__ == "__main__":
    main()
