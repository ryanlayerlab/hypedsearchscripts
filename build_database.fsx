
#r "nuget: Newtonsoft.Json"

open System
open System.IO
open Newtonsoft.Json

//FASTA PARSER EQUIV
type FastaEntry = {Description:String; Sequence:String}

let chunk lines =
    let step (blocks, currentBlock) s =
        match s with
        | "" -> (List.rev currentBlock :: blocks), []
        | _ -> blocks, s :: currentBlock
    let (blocks, lastBlock) = Array.fold step ([], []) lines
    List.rev (lastBlock :: blocks)

let generateFastaEntry (chunk:String List) =
    match chunk |> Seq.length with
    | 0 -> None
    | _ ->
        let description = chunk |> Seq.head
        let sequence = chunk |> Seq.tail |> Seq.reduce (fun acc x -> acc + x)
        Some {Description=description; Sequence=sequence}

let parseFastaFileContents fileContents = 
    chunk fileContents
    |> Seq.map(fun c -> generateFastaEntry c)
    |> Seq.filter(fun fe -> fe.IsSome)
    |> Seq.map(fun fe -> fe.Value)

//BUILD DATABASE

type AminoAcid = {Id:String; Weight:float}
type Protein = {Name:string; AminoAcids: AminoAcid List}
type ProteinMatch = {ProteinName:string; Weight:float; StartIndex:int; EndIndex:int}
type WeightProteinMatches = {Weight:float; ProteinMatchs:ProteinMatch list}

let extractProteinName (description:string) =
    description.Split('|')
    |> Seq.skip 1
    |> Seq.head

let generateAminoAcid (character:char) =
    match character with
    | 'A' -> {Id="A"; Weight=71.037114}| 'R' -> {Id="R"; Weight=156.101111} 
    | 'N' -> {Id="N"; Weight=114.042927} | 'D' -> {Id="D"; Weight=115.026943} 
    | 'C' -> {Id="C"; Weight=103.009185} | 'E' -> {Id="E"; Weight=129.042593} 
    | 'Q' -> {Id="Q"; Weight=128.058578} | 'G' -> {Id="G"; Weight=57.021464} 
    | 'H' -> {Id="H"; Weight=137.058912} | 'I' -> {Id=":I"; Weight=113.084064} 
    | 'L' -> {Id="L"; Weight=113.084064} | 'K' -> {Id="K"; Weight=128.094963} 
    | 'M' -> {Id="M"; Weight=131.040485} | 'F' -> {Id="F"; Weight=147.068414} 
    | 'P' -> {Id="P"; Weight=97.052764} | 'S' -> {Id="S"; Weight=87.032028} 
    | 'T' -> {Id="T"; Weight=101.047679}| 'U' -> {Id="U"; Weight=150.95363} 
    | 'W' -> {Id="W"; Weight=186.079313} | 'Y' -> {Id="Y"; Weight=163.06332} 
    | 'V' -> {Id="V"; Weight=99.068414} | 'X' -> {Id="X"; Weight=0.0} 
    | 'B' -> {Id="B"; Weight=113.084064} | 'Z' -> {Id="Z"; Weight=0.0}
    | _ -> {Id="Z"; Weight=0.0}

let extractAminoAcids(sequence:string) =
    sequence.ToUpper().ToCharArray()
    |> Seq.map(fun c -> generateAminoAcid c)
    |> Seq.toList

let generateProtein(fastaEntry:FastaEntry)=
    let name = extractProteinName fastaEntry.Description
    let aminoAcids = extractAminoAcids fastaEntry.Sequence
    {Name=name;AminoAcids=aminoAcids}

let generateProteins parsedFileContents =
    parsedFileContents
    |> Seq.map(fun fc -> generateProtein fc)
    |> Seq.toList

let generateProteinMatch (protein: Protein) (index:int) (aminoAcids:AminoAcid array) (kmerLength:int) =
    let name = protein.Name
    let weight = 
        aminoAcids 
        |> Array.map(fun aa -> aa.Weight)
        |> Array.reduce(fun acc x -> acc + x)
    let startIndex = index * kmerLength 
    let endIndex = index * kmerLength + kmerLength      
    {ProteinName = name; Weight = weight; StartIndex = startIndex; EndIndex = endIndex}

let generateProteinMatches (protein: Protein) (kmerLength:int) =
    protein.AminoAcids
    |> Seq.windowed(kmerLength)
    |> Seq.mapi(fun idx aa -> generateProteinMatch protein idx aa kmerLength)

let generateAllProteinMatches (proteins: Protein list) (kmerLength:int) =
    proteins 
    |> Seq.map(fun p -> generateProteinMatches p kmerLength)
    |> Seq.reduce(fun acc pm -> Seq.append acc pm)

let generateWeightProteinMatches (proteinMatches: ProteinMatch list)=
    proteinMatches
    |> Seq.map(fun pm -> {ProteinName=pm.ProteinName;StartIndex=pm.StartIndex;EndIndex=pm.EndIndex;Weight=System.Math.Round(pm.Weight,2)})
    |> Seq.groupBy(fun pm -> pm.Weight)
    |> Seq.map(fun (w,pm) -> {Weight= w; ProteinMatchs = pm |> Seq.toList })

let createJsonFile kmerLength =
    let inputFilePath = @"/Users/jamesdixon/Documents/Hypedsearch_All/hypedsearch/data/database/sample_database.fasta"
    let fileContents = File.ReadLines(inputFilePath) |> Seq.cast |> Seq.toArray
    let parsedFileContents = parseFastaFileContents fileContents
    let proteins = generateProteins parsedFileContents
    let proteinMatches = generateAllProteinMatches proteins kmerLength
    let weightProteinMatches = generateWeightProteinMatches (proteinMatches |> Seq.toList)
    let json = JsonConvert.SerializeObject(weightProteinMatches);
    let outputFilePath = "/Users/jamesdixon/Desktop/Hypedsearchservice/weight_protein_matches/all_weight_protein_matches_kmer" + kmerLength.ToString() + ".json"
    File.WriteAllText(outputFilePath,json)
    

[2..9]
|> Seq.iter(fun x -> createJsonFile x)
 
