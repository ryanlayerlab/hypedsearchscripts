
from numpy.lib.function_base import append
import pandas as pd
from pyteomics import fasta
from collections import namedtuple
import json

Amino_Acid = namedtuple('Amino_Acid', ['Id', 'Weight' ],defaults=['',0.0])
Protein = namedtuple('Protein',['Name', 'Amino_Acids'],defaults=['', []])
Protein_Match = namedtuple('Protein_Match',['Protein_Name', 'Weight', 'Start_Index', 'End_Index'], defaults=['',0,0,0])
Weight_Protein_Matches = namedtuple('Weight_Protein_Matches',['Weight','Protein_Match'],defaults=[0,[]])

def extract_protein_name(description):
    return description.split('|')[-1].split()[0]

def generate_amino_acid(character):
    weights = {
    "A": 71.037114,"R": 156.101111,"N": 114.042927,"D": 115.026943,"C": 103.009185,
    "E": 129.042593,"Q": 128.058578,"G": 57.021464,"H": 137.058912,"I": 113.084064,
    "L": 113.084064,"K": 128.094963,"M": 131.040485,"F": 147.068414,"P": 97.052764,
    "S": 87.032028,"T": 101.047679,"U": 150.95363,"W": 186.079313,"Y": 163.06332,
    "V": 99.068414,"X": 0, "B": 113.084064, "Z": 0 }
    weight = weights.get(character)
    return Amino_Acid(character,weight)

def extract_amino_acids(sequence):
    amino_acids = []
    for character in sequence:
        amino_acid = generate_amino_acid(character)
        amino_acids.append(amino_acid)
    return amino_acids      

def generate_proteins(fasta_parser):
    proteins = []
    for entry in fasta_parser:
        protein_name = extract_protein_name(entry.description)
        amino_acids = extract_amino_acids(entry.sequence)
        protein = Protein(protein_name,amino_acids)
        proteins.append(protein)
    return proteins

def get_cumulative_weights(amino_acids, kmer_length):
    df_all = pd.DataFrame(amino_acids)
    df_weights = df_all.loc[:,'Weight']
    windows = df_weights.rolling(kmer_length).sum()
    no_nan_windows = windows.fillna(0)
    rounded_windows = no_nan_windows.apply(lambda x: round(x,2))
    return rounded_windows

def generate_protein_match(protein, data_tuple, kmer_length):
    protein_name = protein.Name
    (start_index, weight) = data_tuple
    end_index = start_index + kmer_length
    return Protein_Match(protein_name,weight, start_index,end_index)

def get_protein_matches(protein, kmer_length):
    protein_matches = []
    cumulative_weights = get_cumulative_weights(protein.Amino_Acids,kmer_length)
    indexes = cumulative_weights.index.tolist()
    values = cumulative_weights.values.tolist()
    data_tuples = list(zip(indexes,values))
    for data_tuple in data_tuples:
        protein_match = generate_protein_match(protein, data_tuple, kmer_length)
        protein_matches.append(protein_match)
    return protein_matches

def generate_proteins_matches(proteins,kmer_length):
    proteins_matches = []
    for protein in proteins:
        protein_matches = get_protein_matches(protein,kmer_length)
        proteins_matches = proteins_matches + protein_matches
    return proteins_matches

def handle_weight_protein_matches(weight_protein_matches, all_weight_protein_matches):
    exists = False
    for item in all_weight_protein_matches:
        if item.Weight == weight_protein_matches.Weight:
            item.Protein_Match.append(weight_protein_matches)
            exists = True
            break
    if exists == False:
        all_weight_protein_matches.append(weight_protein_matches)

def generate_all_weight_protein_matches(protein_matches):
    all_weight_protein_matches = []
    for protein_match in protein_matches:
        weight = protein_match.Weight
        weight_protein_matches = Weight_Protein_Matches(weight, [protein_match])
        handle_weight_protein_matches(weight_protein_matches,all_weight_protein_matches)
    return all_weight_protein_matches

def create_json_file(kmer_length):
    file_path = '/Users/jamesdixon/Documents/Hypedsearch_All/hypedsearch/data/database/sample_database.fasta'
    fasta_parser = fasta.read(file_path) 
    proteins = generate_proteins(fasta_parser)
    proteins_matches = generate_proteins_matches(proteins,kmer_length) 
    all_weight_protein_matches = generate_all_weight_protein_matches(proteins_matches) 
    serialzed = json.dumps(all_weight_protein_matches)
    text_file = open("all_weight_protein_matches_kmer" + str(kmer_length) + ".json", "w")
    n = text_file.write(serialzed)
    text_file.close()
    print("Done")

create_json_file(9)
#n = 279 for fasta
#n=103163 for kmer=2 combined to n=176 for kmer=2 and round to 2 decimal places
