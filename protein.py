import protein_tests as test
import json
### PHASE 1 ###

def read_file(filename: str) -> str:
    """ Read the file concatenate strings while stripping newline characters.

    Args:
        filename (str): The file name
    """
    dna = ""
    with open(filename, "r") as file:
        for line in file.readlines():
            dna += line.strip("\n")
        file.close()
    return dna

def t_to_u(codon:str):
    u_base_codon=""
    for j in range(3):
        if codon[j] == "T":
            u_base_codon+="U"
        else:
            u_base_codon+=codon[j]
    return u_base_codon

def dna_to_rna(dna, start_index):
    # One codon, ATG, signals the Start of an RNA strand; three other codons (TAA, TAG, and TGA) signal the end (Stop)
    # stop codons UAA, UAG, and UGA.
    rna = []
    i = start_index
    for i in range(start_index, len(dna), 3):
        current_codon = dna[i:i+3]
        if "T" in current_codon:
            rna.append(t_to_u(current_codon))
        else:
            rna.append(current_codon)
        if rna[len(rna)-1] in ["UAA", "UAG", "UGA"]:
            return rna
    return rna

def make_codon_dictionary(filename):

    jsonfile=open(filename)
    codon_table=json.load(jsonfile)
    jsonfile.close()
    new_codon_table = {}
    for amino_acid in codon_table.keys():
        for codon in codon_table[amino_acid]:
            if "T" in codon:
                new_codon_table[t_to_u(codon)]=amino_acid
            else:
                new_codon_table[codon]=amino_acid
    return new_codon_table

def generate_protein(codons, codon_dict):
    protein = []
    for idx in range(len(codons)):
        if idx == 0:
            protein.append("Start")
        else:
            protein.append(codon_dict[codons[idx]])
    return protein

def synthesize_proteins(dna_filename, codon_filename):
    """dna_filename uses read_file and codon_file uses make codon dictionary

    
    -[x] Implement the function synthesize_proteins(dna_filename, codon_filename) in the starter file. This program should read the DNA from the given filename (using read_file) and produce a codon dictionary by calling make_codon_dictionary.



    -[] The program should then identify all of the RNA strands that can be produced from the DNA by iterating through all the indexes in the DNA string, looking for the start code (ATG) at each point. Note that you’ll need to keep track of a list of proteins and a count variable outside of the loop.

    -[] If you identify an index in the DNA that corresponds to ATG, call dna_to_rna starting from that index to extract the entire RNA sequence, then call generate_protein on the resulting RNA (and codon dictionary) to produce a protein. That protein should be added to an overall protein list. 
    
    -[] Then update the index in the DNA strand to skip past all the already-checked bases (by adding 3 * the length of the RNA strand).
    If you get to an index that does not correspond to ATG, move on to the next base.

    When you finish looping, you’ll have a list of all the proteins synthesized from the DNA. Return that as your final result.

    """
    # read the DNA from the given filename (using read_file)
    # produce a codon dictionary by calling make_codon_dictionary.
    # The program should then identify all of the RNA strands that can be produced from the DNA by iterating through all the indexes in the DNA string, looking for the start code (ATG) at each point. Note that you’ll need to keep track of a list of proteins and a count variable outside of the loop.
    dna_file = read_file(dna_filename)

    codons_mapping = make_codon_dictionary(codon_filename)
    proteins=[]
    idx=0
    while idx < len(dna_file):
    # for idx in range(len(dna_file)):
        if dna_file[idx:idx+3] == "ATG":
            protein = generate_protein(dna_to_rna(dna=dna_file, start_index=idx),codon_dict=codons_mapping)
            idx += 3*len(protein)
            proteins.append(protein)
            print(dna_file[idx:idx+3])
        else:
            idx+=1
    return proteins
    # print("-------------------------------------------------")
    # print("-------------------------------------------------")
    # If you identify an index in the DNA that corresponds to ATG, call dna_to_rna starting from that index to extract the entire RNA sequence, then call generate_protein on the resulting RNA (and codon dictionary) to produce a protein. That protein should be added to an overall protein list. 


    # Then update the index in the DNA strand to skip past all the already-checked bases (by adding 3 * the length of the RNA strand).
    # If you get to an index that does not correspond to ATG, move on to the next base.

### PHASE 2 ###

def common_proteins(protein_list1, protein_list2):
    unique_proteins=[]
    for protein_of1 in protein_list1:
        if protein_of1 in protein_list2 and protein_of1 not in unique_proteins:
            unique_proteins.append(protein_of1)
    return unique_proteins

def combine_proteins(protein_list):
    all_proteins=[]
    for sublist in protein_list:
        all_proteins+=sublist
    return all_proteins

def amino_acid_dictionary(aa_list):
    amino_acid_count={}
    for amino_acid in aa_list:
        # if amino_acid
        if amino_acid_count.get(amino_acid):
            amino_acid_count[amino_acid]+=1
        else:
            amino_acid_count[amino_acid]=1
    return amino_acid_count

"""12. [15 pts] Find Amino Acid Differences
Now that we know how common each amino acid is, we can start comparing amino acids between genes of different lengths. Because the genes have different lengths, we can’t just compare the counts of amino acids. Instead, we’ll compare the frequencies of amino acids- in other words, how frequently it occurs in the gene. We can determine the frequency of an amino acid by finding its count (from amino_acid_dictionary) and dividing it by the total number of amino acids in the gene.
Implement the function find_amino_acid_differences(protein_list1, protein_list2, cutoff) in the starter file. This takes two protein lists and a float cutoff and returns a list of three-element lists, where the first element in the list is an amino acid, the second element is the frequency of that amino acid in protein_list1, and the third element is the frequency of that amino acid in protein_list2. You should only include amino acids in this returned list if the difference between their frequencies is greater than the provided cutoff. This cutoff is given as a decimal- in other words, 0.02 is 2%. You should also not include the Start and Stop amino acids in the list, as they are not interesting for this analysis (though they should still count towards the overall length of the gene).
To generate this list, you should first use your combine_proteins and amino_acid_dictionary functions to generate data about amino acid frequencies for each protein list. Then go through each amino acid in the lists, and add each amino acid if and only if the two frequencies are sufficiently different between the two genes. If an amino acid does not occur in one of the two lists, its frequency is 0.
"""

def find_amino_acid_differences(protein_list1, protein_list2, cutoff):
    print("-------------------------------------------------")
    print("-------------------------------------------------")
    flat_protein_list1 = combine_proteins(protein_list1)
    list1_length = len(flat_protein_list1)
    print(len(flat_protein_list1))
    amino_acid_dict1 = amino_acid_dictionary(flat_protein_list1)
    print(amino_acid_dict1)
    print("-------------------------------------------------")
    print("-------------------------------------------------")
    flat_protein_list2 = combine_proteins(protein_list2)
    list2_length = len(flat_protein_list2)
    print(len(flat_protein_list2))
    amino_acid_dict2 = amino_acid_dictionary(flat_protein_list2)
    print(amino_acid_dict2)
    print("-------------------------------------------------")
    print("-------------------------------------------------")
    return

def display_commonalities(commonalities):
    return

def display_differences(differences):
    return

### RUN CODE ###

if __name__ == "__main__":
    # test.test_read_file()
    # test.test_dna_to_rna()
    # test.test_make_codon_dictionary()
    # test.test_generate_protein()
    # test.test_synthesize_proteins()
    # test.test_common_proteins()
    # test.test_combine_proteins()
    # test.test_amino_acid_dictionary()
    test.test_find_amino_acid_differences()
    # test.test_all()
    # test.run()