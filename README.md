 
NetSyn is a tool to detect conserved genomic contexts (i.e. synteny conservation) among a list of protein targets. Synteny are computed using a re-implementation of the method described in Boyer *et al.* article ([https://doi.org/10.1093/bioinformatics/bti711](https://doi.org/10.1093/bioinformatics/bti711)).

  # Installation
Prerequisites: python 3.6.5 or 3.7 and MMseqs2.

Netsyn installation: `python3 install setup.py`
 
For the installation of MMseqs2 please refer to https://github.com/soedinglab/MMseqs2/wiki#installation
 
 Python library list:

-   pyyaml
    
-   python-igraph
    
-   jsonschema
    
-   networkx
    
-   markov_clustering
    
-   urllib3
    
-   biopython
    
-   scipy==1.1.0
 
 
# Basic Usage

NetSyn can be used with 2 different input file formats. One is a file containing a list of UniProt accessions (`-u`), while the other one is a file of correspondences (`-c`). The two types of file are described in the Input Data part. It is possible to start an analysis with both input file formats. It leads to 3 NetSyn basic usage callings:

-   with the UniProt accessions list:
 
	`netsyn -u <UniProtAC.list> -o <OutputDirName>`

-   with the correspondences file:
    
	`netsyn -c <CorrespondencesFileName> -o <OutputDirName>`

-   with both entries:
    
	`netsyn -u <UniProtAC.list> -c <CorrespondencesFileName> -o <OutputDirName>`

  

Whatever the type of input, it is necessary to provide an output directory (-o). NetSyn creates the directory and stores all the results in it.

# Settings

### General Settings

-   `-o/--OutputDirName`: Name of the directory which is going to be created into the working directory from where the analysis starts. This directory must not exist already. It is the only one required parameter of this section.
    
-   `-u/--UniProtACList`: Name of the file containing the Uniprot accession list (full description in the Input Data part).
    
-   `-c/--CorrespondencesFile`: File of correspondences between protein accession, nucleic accession and the location of the genomic files (full description in the Input Data part)
    
-   `-md/--MetaDataFile`: file containing metadata information (full description in the Input Data part)
    
-   `-np/--newProject`: force the creation of a new project. If a directory or NetSyn project already exists under the same name as provided with --OutputDirName parameter, this one is overwritten
    
-   `--conserveDownloadedINSDC`: conserve the downloaded INSDC files. By default, the downloaded INSDC files are removed. This parameter can be only used with the --UniProtACList option
    

### Clustering Settings

These parameters control the call to MMseqs2. The call to MMseqs2 is detailed in the Dependencies part.

-   `-id/--Identity`: minimum sequence identity (default value: 0.3).
    
-   `-cov/--Coverage`: minimum sequence coverage (default value: 0.8).


----------


    

### Synteny Settings

-   `-ws/--WindowSize`: size of the genomic contexts to compare. The value must be an odd number between 3 to 11. The target is considered in the middle of the genomic context. If the target is close to a border of the contig, the larger existing genomic context is taken. By default, the larger genomic context is taken (11)
    
-   `-sg/--SyntenyGap`: minimal number of genes without homologue in the second synteny between two genes with homologue genes into the second synteny to be considered as part of the synteny. The higher the value is , the less stringent the definition of conserved synteny is. Target genes do not enter into consideration in this setting. The default value is equal to 3.
    

![](https://lh4.googleusercontent.com/v836YNzK32giBIglaOIw78Dm0eD_Vwj_nvKBvHYhupDE38kOrb-Cas2GaIMO9muNWpYUuss8Lk61gFQ_XG-3wP2s78QMcdAZIpozCP30c-5AQqa0xjrUhZaRr78gM-nz7vcJeYtb)

-   `-ssc/--SyntenyScoreCutoff`: score lower threshold between two synteny mandatory to create an edge between the two target gene in the graph. By default, the minimum threshold is equal to 3
    

  

### Graph reduction settings

In order to reduce the graph complexity target nodes identified as belonging to the same synteny cluster (several clustering methods are available) and sharing a given property (indicated by the `-gl` option) will be merged into a unique node. The merging can be applied on only one property at a time. It is not allowed to specify a taxonomic rank (with the `-gt` option) and a metadata label (with the `-gl` option) in a single analysis. In order to define which synteny cluster repartition to use to compute the redundancy removal, a clustering method (`-cm`) must be provided.

-   `-cm/--ClusteringMethod`: clustering method that will provide the synteny clusters repartition in order to determine the nodes to merge according to a given property. Clustering methods are: {MCL, Infomap, Louvain, WalkTrap}
    
-   `-gt/--GroupingOnTaxonomy`: taxonomic rank on which the grouping clustering will be computed. A choice is given on the list of taxonomic ranks retrieved from NCBI taxonomy request: {superkingdom, phylum, class, order, family, genus, species} (see more on Web Requests part). Only one rank can be specified at a time. This option is not compatible with the `--GroupingOnLabel` option but requires the `--ClusteringMethod` option to make the functionality enabled
    
-   `-gl/--GroupingOnLabel`: label taken from the list of metadata labels provided by the user. It requires the use of the K option. The name must be equal to the header of the provided metadata file. Names “accession_type” and “accession” are not available for the option. Only one label can be specified at a time. This option is not compatible with the `--GroupingOnTaxonomy` option but requires the `--ClusteringMethod` option to make it enabled
    

  

### Advanced Settings

Some additional settings for MMseqs or graph clustering methods can be specified.

These settings are transmitted through two YAML files as follows.

#### MMseqs advanced settings

-   `-mas/--MMseqsAdvancedSettings`: YAML file name.
    
-   Example of YAML file with MMSeqs default advanced settings:
    
```yaml
MMseqs advanced settings:
	MMseqs_cov-mode: 1
	MMseqs_cluster-mode: 1
	MMseqs_kmer-per-seq: 80
	MMseqs_max-seqs: 300
	MMseqs_single-step-clustering: false
	MMseqs_threads: 4
```

#### Graph clustering methods advanced settings

-   `-asc/--AdvancedSettingsClustering`: YAML file name.
    
-   Example of YAML file with clustering method default advanced settings:
    
```yaml
MCL advanced settings:
	MCL_inflation: 2
	MCL_expansion: 2
	MCL_iterations: 1000
	
WalkTrap advanced settings:
	walktrap_step: 4
	
Infomap advanced settings:
	infomap_trials: 10
```

  

# Input Data

This part is dedicated to the description of the files used as input for NetSyn. They are provided by the user. NetSyn can take two kind input files a user can provide with the `--UniProtACList` or K option and one metadata file with the `--MetaDataFile` option.

## UniProt Accessions list

This file must be written as a one column-like. The head column name must be UniProt_AC. The next lines contain only one value, a UniProt accession.
```
UniProt_AC
A0A090WTF7
D7NAF6
K1S2E2
A0A1Q5ZYL9
A0A1J5SHC1
...
A0A1I1VZJ3
A0A0Q0SHG0
A0A1C5UBS1
A0A1H4FNN6
```
  

Starting from the UniProt AC in the list, NetSyn sends a request to the UniProt website (see more on Web Requests part) to get the EMBL protein_id and the EMBL nucleic_id accessions for every UniProt accession. With these identifiants NetSyn is able to download the INSDC file where the genomic context for each UniProt_AC is computed.

  

Some target (sequences given by the user) may be loose at this stage. If there is no correspondence between the UniProt_AC and a INSDC file, the target sequence will not be taken into account for the rest of the analysis and will not be retrieved into the final graph.


## File of Correspondences

The file of correspondences is created by NetSyn when using a UniProt accessions list as input. However, the user may have his own data, which are not stored on UniProt. This is the reason why the user has the possibility to start an analysis with different kind of inputs, but must create the correspondence file by his own.

  

The correspondences file requires 5 columns separated by tabulations (UniProt_AC can be considered as an optional column):

  
```
UniProt_AC	protein_AC	protein_AC_field	nucleic_AC	 nucleic_File_Format	nucleic_File_Path
B3QMV9	ACF11262.1	protein_id	CP001099	embl	data/CP001099.embl
A2VZ04	EAY64950.1	protein_id	CH482378	embl	data/CH482378.embl
B0MVB0	EDS03995.1	protein_id	ABFK02000017	embl	data/ABFK02.embl
A8F3R2	ABV32796.1	protein_id	CP000812	embl	data/CP000812.embl
A3P412	ABN94303.1	protein_id	CP000573	embl	data/CP000573.embl
A3P387	ABN95365.1	protein_id	CP000573	embl	data/CP000573.embl
A6G0K8	EDM80654.1	protein_id	ABCS01000009	embl	data/ABCS01.embl
A3P8Z7	ABN93085.1	protein_id	CP000573	embl	data/CP000573.embl
C5AU40	ACS42739.1	protein_id	CP001510	embl	data/CP001510.embl
```
  -   protein_AC: protein accession that allows NetSyn to identify the target protein
    
-   protein_AC_field: name of the field where the protein_AC is provided for every protein. Syntax to used if field is a dbxref: “dbxref:MaGe” (MaGe is the desired dbxref name)
    
-   nucleic_AC: identifier of the contig that contains the protein_AC, unless identifier of the genome in the file contains the whole assembled genome
    
-   nucleic_File_Format: format of INSDC file. NetSyn2 support file formats that the BioPython library (see more in Dependences part) is able to parse: .embl (embl), .gbff et .gbk (genbank or gb)
    
-   nucleic_File_Path: relative or absolute path where NetSyn2 will find the INSDC file to parse
    
-   UniProt_AC: UniProt accession of the protein. This column is optional and can be filled with “NA” values. If the UniProt accession stored in the INSDC file differs from the one provided by the user (unless ‘NAs’), the UniProt accessions of the user have the priority and a Warning message is printed
    

## Metadata file

The user has the possibility to add information on the sequences computed. These informations can be map on the final graph.

The user adds information on target proteins via the metadata file.
  
Various information may be added for every target proteins, or only for a subset of these target proteins.

Several metadata may be specified to a subset of target proteins.

It consists of 2 required columns in order to specify the concerned target protein and as many columns as the existing metadata fields.
```
accession_type	accession	metada_1	medata_n
UniProt_AC	A0A090WTF7	2	classification_4
UniProt_AC	D7NAF6	1	classification_1
UniProt_AC	K1S2E2	2	classification_2
UniProt_AC	A0A1Q5ZYL9	NA	classification_3
UniProt_AC	A0A1H4FNN6	1	classification_1
UniProt_AC	A0A1V9GD33	2	classification_4
```
-   accession_type: according to the origin of the target protein ("UniProt_AC" if contained in the input file `--UniProtACList` or "Protein_AC" if contained in the input file `--CorrespondencesFile`).
    
-   accession: accession of the protein used to identify the target protein.
    
-   Any other column useful for characterizing a target protein by a metadata. In this example, two metadata ("metadata_1" and "metadata_2") have been used. "NA" is the default value if the metadata value unknown for one protein target.
    

# Data and multiple analysis inside a project

A NetSyn project is defined by the `-o/--outputDirName` option (output directory) and the input file(s) (`-u/--UniProtList` and/or `-c/--CorrespondencesFiles`).

  

A NetSyn analyse corresponds to a NetSyn run with specific parameters. One NetSyn parameter change will launch a new analyse from the corresponding step in the change. The results of each analysis will be in a directory inside the output directory. Each of these analysis directory will contain:

-   A result file in graphML format
    
-   A file to open with your web browser to explore the results
    
-   A resume file of the analyse
    
-   A summary file of the parameters used
    
-   A directory with synthesis files by NetSyn clustering
    

Besides the results files, NetSyn2 creates some intermediate files. NetSyn2 might be separated into 5 steps: 1) GetINSDCFiles, 2) ParseINSDCFiles_GetTaxonomy, 3) ClusteringIntoFamilies, 4) SyntenyFinder and then 5) DataExport. At the end of each part, a check on the generated files is done. It is possible to launch each of these steps independently. Below the details of the input file by step:

### GetINSDCFiles step
    

-   The UniProt accession list, same as netsyn  `--UniProtList` option.
    

### ParseINSDCFiles_GetTaxonomy step
    

-   The file of correspondences, same as netsyn  `--CorrespondencesFile` option.
    

### ClusteringIntoFamilies step
    

-   The file of protein sequences in fasta format.
    
	```
	>828 // protein unique identifier from proteins_parsingStep.json
	MNDQLFKKVLGYIESESYLMAYRELHKLADEYMPLATRMDFDALHSSLSIIIGERSGYPDIADQLADTAGFYERLAYLLTKKLLGDDEAGEKADTLMLCVVAFGNHRRN
	```
  

-   The file of protein data in json format (see … for more information).
	    
	```json
	[
		{
			"id": "828", // protein unique identifier
			"protein_AC": "ACF11257.1",
			"begin": 916518,
			"end": 916847,
			"strand": "1",
			"products": "conserved hypothetical protein",
			"ec_numbers": "NA",
			"UniProt_AC": "B3QMV4",
			"gene_names": "NA",
			"locus_tag": "Cpar_0841",
			"targets": ["833"], // protein unique identifier list from this file
			"targets_idx": ["5"] // protein index list from this file
		}
	]
	```
### Synteny Finder step
    

-   The file of protein data in json format (see … for more information).
    
	```json
	[
		{
			"id": "828", // protein unique identifier
			"protein_AC": "ACF11257.1",
			"begin": 916518,
			"end": 916847,
			"strand": "1",
			"products": "conserved hypothetical protein",
			"ec_numbers": "NA",
			"UniProt_AC": "B3QMV4",
			"gene_names": "NA",
			"locus_tag": "Cpar_0841",
			"targets": ["833"], // protein unique identifier list from this file
			"targets_idx": ["5"], // protein index list from this file
			"family": 453  // family identifier
		}
	]
	```
-   The file of target data in json format (see … for information).
    
```json
{
	"5": { // protein index from proteins_parsingsStep.json
		"id": "833",
		"organism_id": 1,
		"context": ["828","829 ", "830", "831", "832", "833", "834", "835", "836", "837", "838"], // proteins unique identifier list from proteins_parsingsStep.json
		"context_idx": ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"], // proteins index list from proteins_parsingsStep.json
		"UniProt_AC": "B3QMV9",
		"protein_AC": "ACF11262.1",
		"organism_idx": 0  // organism index from organisms_taxonomyStep.json
	}
}
```
  

### DataExport step
    

-   The file of node data in json format (see … for more information).
    
	```json
	[
		{
			"target_idx": "5", // protein index from proteins_familiesStep.json
			"id": "1200", // protein unique identifier from proteins_familiesStep.json
			"UniProt_AC": "A2VZ04",
			"protein_AC": "EAY64950.1",
			"context": ["1195", "1196", "1197", "1198", "1199", "1200", "1201", "1202", "1203", "1204", "1205"],// proteins unique identifier list from proteins_familiesStep.json
			"context_idx": ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"], // proteins index list from proteins_familiesStep.json
			"organism_id": 2,
			"organism_idx": 1,
			"clusterings": { "WalkTrap": 0, "Louvain": 0, "Infomap": 2, "MCL": 0 }, // cluster identifier by clustering methods
			"families": [1064, 2154, 171, 237, 47, 1809, 1458, 883, 756, 565, 2230], // families identifier list same proteins_familiesStep.json
			"Size": 1
		}
	]
	```
-   The file of edge data in json format (see … for more information).
    
	```json
	[
		{
			"source": "0", // nodes index from nodes_list.json
			"target": "1", // nodes index from nodes_list.json
			"proteins_idx_source": ["6", "10", "9", "3", "4", "5"], // proteins index list from nodes_list.json
			"proteins_idx_target": ["48", "47", "53", "51", "50", "49"], // proteins index list from nodes_list.json
			"weight": 4.800000000000001
		}
	]
	```
-   The file of protein data in json format (see … for more information).
    
```json
[
	{
		"id": "1195", // protein unique identifier
		"protein_AC": "EAY64945.1",
		"begin": 452015,
		"end": 452842,
		"strand": "-1",
		"products": "Arginine 3rd transport system periplasmic binding protein",
		"ec_numbers": "NA",
		"UniProt_AC": "A2VYZ9",
		"gene_names": "NA",
		"locus_tag": "BCPG_03285",
		"targets": ["1200"], // protein unique identifier list from this file
		"targets_idx": ["5"], // protein index list from this file
		"family": 2230  // family identifier
	}
]
```
-   The file of organism data in json format (see … for more information).
    
	```json
	[
		{
			"id": 1, // unique identifier
			"name": "Chlorobaculum parvum NCIB 8327",
			"strain": "NCIB 8327",
			"taxon_id": "517417",
			"targets_idx": ["5"],
			"lineage": [
				{
					"rank": "superkingdom",
					"scientificName": "Bacteria",
					"tax_id": "2",
					"level": 1
				},
				{
					"rank": "phylum",
					"scientificName": "Chlorobi",
					"tax_id": "1090",
					"level": 5
				},
				{
					"rank": "class",
					"scientificName": "Chlorobia",
					"tax_id": "191410",
					"level": 8
				},
				{
					"rank": "order",
					"scientificName": "Chlorobiales",
					"tax_id": "191411",
					"level": 13
				},
				{
					"rank": "family",
					"scientificName": "Chlorobiaceae",
					"tax_id": "191412",T
					"level": 17
				},
				{
					"rank": "genus",
					"scientificName": "Chlorobaculum",
					"tax_id": "256319",
					"level": 21
				},
				{
					"rank": "species",
					"scientificName": "Chlorobaculum parvum",
					"tax_id": "274539",
					"level": 26
				}
			]
		}
	]
	```
# Web Requests

UniProt: allows to recover the protein accession and nucleic accession from a UniProt accession (into the GetINSDCFiles part).

  

EBI-ENA: allows to recover the INSDC file (embl format) from a nucleic accession (into the GetINSDCFiles part).

  

NCBI-taxonomy: allows to recover the lineage taxonomic from a toxon identifier (into the ParseINSDCFles_GetTaxonomy part).

