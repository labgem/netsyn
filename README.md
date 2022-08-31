# NetSyn: detect synteny conservation among a list of protein targets

NetSyn is a tool to detect conserved genomic contexts (i.e. synteny conservation) among a set of protein (call protein targets). Synteny are computed using a re-implementation of the method described in Boyer *et al.* article ([https://doi.org/10.1093/bioinformatics/bti711](https://doi.org/10.1093/bioinformatics/bti711)).

## Installation

Prerequisites: python 3.8 and MMseqs2.

Netsyn installation: `python3 setup.py install`

For the installation of MMseqs2 please refer to https://github.com/soedinglab/MMseqs2/wiki#installation

 Python library list:

  - pyyaml

  - python-igraph

  - jsonschema

  - networkx (>= 2.8)

  - markov_clustering

  - urllib3

  - biopython

  - requests

You can easely install NetSyn using an virtual environment with pip (command lines example below):

  - virtualenv venv_netsyn

  - source venv_netsyn/bin/activate

  - pip install pyyaml

  - pip install python-igraph

  - pip install jsonschema

  - pip install networkx

  - pip install markov_clustering

  - pip install urllib3

  - git clone https://github.com/labgem/netsyn

  - python3 setup.py install

## Basic Usage

NetSyn can be used with 2 different input file formats. One is a file containing a list of UniProt accessions (`-u` option), while the other one is a correspondences file (`-c` option). The two types of file are described in the Input Data part. It is possible to start an analysis with both input file formats. It leads to 3 NetSyn basic usage callings:

  - with the UniProt accessions list:
    `netsyn -u <UniProtAC.list> -o <OutputDirName>`

  - with the correspondences file:
    `netsyn -c <CorrespondencesFileName> -o <OutputDirName>`

  - with both entries:
    `netsyn -u <UniProtAC.list> -c <CorrespondencesFileName> -o <OutputDirName>`

Whatever the type of input, it is necessary to provide an output directory name (-o) which will be created and where NetSyn will stores all the results.

## Settings

### General Settings

  - `-o/--OutputDirName`: Name of the directory which is going to be created into the working directory and where the analysis  will start. This directory must not already exist. This option is required.

  - `-u/--UniProtACList`: Name of the file containing the Uniprot accession list (full description in the Input Data part).

  - `-c/--CorrespondencesFile`: File of correspondences between protein accession, nucleic accession and the location of the genomic files (full description in the Input Data part)

  - `-md/--MetaDataFile`: file containing metadata information (full description in the Input Data part)

  - `-np/--newProject`: force the creation of a new project. If a directory or NetSyn project already exists under the same name as provided with --OutputDirName parameter, this one is overwritten

  - `--conserveDownloadedINSDC`: conserve the downloaded INSDC files. By default, the downloaded INSDC files are removed. This parameter can be only used with the --UniProtACList option


### Clustering Settings

These parameters control the MMseqs2 call (more details in the Dependencies part).

  - `-id/--Identity`: minimum sequence identity (default value: 0.3).

  - `-cov/--Coverage`: minimum sequence coverage (default value: 0.8).

### Synteny Settings

  - `-ws/--WindowSize`: Number of genes considered into the comparing genomic contexts. The value must be an odd number between 3 to 11. The target protein is considered at the middle of the genomic context. If the target is close to a border of a contig (i.e. near a end of a INSDC file), the larger existing genomic context is taken. By default, the larger genomic context is taken (11)

<p align="center"><img src="/images/windows_size.jpg" width="100%"></p>

  - `-sg/--SyntenyGap`: maximal number of genes without homologue in the second synteny between two genes with homologue genes into the second synteny to be considered as part of the synteny. The definition of conserved synteny is less stringent when this value is higher. The default value is equal to 3.


<p align="center"><img src="/images/gap.jpg" width="100%"></p>

  - `-ssc/--SyntenyScoreCutoff`: score lower threshold between two synteny mandatory to create an edge between the two target gene in the graph. By default, the minimum threshold is equal to 3

### Graph reduction settings

In order to reduce the graph complexity, target nodes identified as belonging to the same synteny cluster (several clustering methods are available) and sharing a given property (indicated by the `-gl` option) will be merged into a unique node. The merging can be applied on only one property at a time. It is not allowed to specify a taxonomic rank (with the `-gt` option) and a metadata label (with the `-gl` option) in a single analysis. In order to define which synteny cluster repartition to use to compute the redundancy removal, a clustering method (`-cm`) must be provided.

  - `-cm/--ClusteringMethod`: clustering method used to group syntenies sharing hight similarity. Several clustering methods are available: {MCL, Infomap, Louvain, WalkTrap, All}. By default this option is set to all. In order to reduce the graph, only one method must be chosen.

  - `-gt/--GroupingOnTaxonomy`: taxonomic rank use to reduce the graph. Nodes belonging to a same cluster and a same taxonomic rank will be merged. A choice is given on the list of taxonomic ranks retrieved from NCBI taxonomy request: {superkingdom, phylum, class, order, family, genus, species} (see more on Web Requests part). Only one rank can be specified at a time. This option is not compatible with the `--GroupingOnLabel` option but requires the `--ClusteringMethod` option to make the functionality enabled

  - `-gl/--GroupingOnLabel`: label taken from the list of metadata labels provided by the user and use to merge nodes in a same cluster. The given name must the same to the header of the provided metadata file. Names “accession_type” and “accession” are not available for this option. Only one label can be specified at the same time. This option is not compatible with the `--GroupingOnTaxonomy` option but requires the `--ClusteringMethod` option.

### Advanced Settings

Some additional settings for MMseqs or graph clustering methods can be specified.

These settings are transmitted through two YAML files as follows.

#### MMseqs advanced settings

  - `-mas/--MMseqsAdvancedSettings`: YAML file name.

  - Example of YAML file with MMSeqs default advanced settings:

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

  - `-asc/--AdvancedSettingsClustering`: YAML file name.

  - Example of YAML file with clustering method default advanced settings:

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

## Input Data

This part is dedicated to the description of the files used as input for NetSyn. They are provided by the user. NetSyn can take two kind input files a user can provide with the `--UniProtACList` or the `--CorrespondencesFile` option and one metadata file with the `--MetaDataFile` option.

### UniProt Accessions list

This file must contain only one column labeled as "UniProt_AC" with one UniProt accession per line.
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

Starting from the UniProt AC in the list, NetSyn sends a request to the UniProt website (see more on Web Requests part) to get the corresponding EMBL protein_id and the EMBL nucleic_id accessions for every UniProt accession. With these identifiants NetSyn is able to download the INSDC file where the genomic context for each UniProt_AC can be retrieved.

Some target (sequences given by the user) may be loose at this stage. If there is no correspondence between the UniProt_AC and a INSDC file, the target sequence will not be taken into account for the rest of the analysis and will not be retrieved into the final graph.

### File of Correspondences

The file of correspondences is created by NetSyn when using a UniProt accessions list as input. However, the user may have his own data, which are not stored on UniProt. This is the reason why the user has the possibility to start an analysis with different kind of inputs, but must create the correspondence file by his own.

The correspondences file requires 6 columns separated by tabulations:

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

  - protein_AC: protein accession that allows NetSyn to identify the target protein

  - protein_AC_field: name of the field where the protein_AC is provided for every protein. Syntax to used if field is a dbxref: “dbxref:MaGe” (MaGe is the desired dbxref name)

  - nucleic_AC: identifier of the contig that contains the protein_AC, or identifier of the genome if the file contains the whole assembled genome

  - nucleic_File_Format: format of INSDC file. NetSyn support file formats that the BioPython library (see more in Dependences part) is able to parse: .embl (embl), .gbff et .gbk (genbank or gb)

  - nucleic_File_Path: relative or absolute path where NetSyn can find the INSDC file to parse

  - UniProt_AC: UniProt accession of the protein. This column is optional and can be filled with “NA” values. If the UniProt accession stored in the INSDC file differs from the one provided by the user (unless ‘NAs’), the UniProt accessions of the user have the priority and a Warning message is printed

### Metadata file

Various information may be added for every target proteins, or only for a subset of them, with a metadata file. These informations will be map on the final graph.

The metadata file consists of 2 required columns in order to specify the concerned target protein and as many columns as metadata fields.
```
accession_type	accession	metada_1	medata_n
UniProt_AC	A0A090WTF7	2	classification_4
UniProt_AC	D7NAF6	1	classification_1
UniProt_AC	K1S2E2	2	classification_2
UniProt_AC	A0A1Q5ZYL9	NA	classification_3
UniProt_AC	A0A1H4FNN6	1	classification_1
UniProt_AC	A0A1V9GD33	2	classification_4
```
  - accession_type: according to the origin of the target protein ("UniProt_AC" if contained in the input file `--UniProtACList` or "Protein_AC" if contained in the input file `--CorrespondencesFile`).

  - accession: accession of the protein used to identify the target protein.

  - Any other column useful for characterizing a target protein by a metadata. In this example, two metadata ("metadata_1" and "metadata_2") have been used. "NA" is the default value if the metadata value unknown for one protein target.



## Output format

NetSyn will create different output files:

  - A .txt: NetSyn log file

  - A .graphML file: final graph which can be read by graph visualization tools like Gephi or cytoscape

  - A .yaml file: summary file of the parameters used

  - A .html file : final graph which can be opened into your web browser to explore the results

<p align="center"><img src="/images/netsyn_app.png" width="100%"></p>

The html file is generated with the D3.js library.The NetSyn web interface is divided into 4 panels. The final graph is displayed into the central panel. The nodes can be colored according to the target proteins attributes like the cluster to which belongs the protein, taxonomy information, metadata given by the users. The upper left panel is the legend of the graph. When a user clicks on a colored spot of the legend, two other  panels appear, one on the right of the graph and the second below. The right panel shows a schematic view of the context of the select nodes. Each context is centered on the target protein with the five genes before and after it. Genes belonging to the same MMSEQ family have the same color. If the users pass his cursor above the one the gene, a pop up shows some information about this gene like the protein_ac, the organism where the gene comes from and annotation found into the INSDC file. The last panel below the graph gives data for families generated by MMSEQ as the number of synteny the family is involved in, the number of protein belonging to this family, the number of species, of strain and genome with protein belonging to this family. These datas concern only the protein of the selected graph cluster.


NetSyn will also create a directory where it will store results on protein families defined by MMseqs2. For each network clustering methods, 2 tabulate files will be created:

(Clustering_Methods_name)_interCluster_family_netsyn.tsv: for each protein family, it give the number of network cluster and their identification number it has been found, if a target protein belongs to, the number of syntenies it is involved, the number of species it has been found, the number of strain it has been found, the number of organisms it has been found and the Product, genes names, EC number, locus tag, protein acession, metadata associated to proteins in this family

(Clustering_Methods_name)_intraCluster_family_netsyn.tsv: for each couple of network cluster and mmseq family, it give the following information: if this family contain a target protein, the number of syntenies this family is involve into this cluster, the number of strain, the number of organism, the list of the organisms, the number of protein in synteny, the product, gene names, EC numbers, locus tag, proteins accession number and metadata of the proteins in this family in this cluster

## Data and multiple analysis inside a project

A NetSyn analyse corresponds to a NetSyn run with specific parameters. If the user want to lauch netsyn with the input, NetSyn will launch a new analyse using, when possible, the already computed results and start from the corresponding step in the change. The results of each analysis will be stored in a new directory created into the same output directory.

Besides the results files, NetSyn creates some intermediate files. NetSyn2 might be separated into 5 steps: 1) GetINSDCFiles, 2) ParseINSDCFiles_GetTaxonomy, 3) ClusteringIntoFamilies, 4) SyntenyFinder and then 5) DataExport. At the end of each part, a check on the generated files is done. It is possible to launch each of these steps independently. Below the details of the input file by step:

### GetINSDCFiles step

At This step, Netsyn dowload the INSDC file from the EBI server for each UniProt accession given as input.

### ParseINSDCFiles_GetTaxonomy step

NetSyn parse the INSDC file dowloaded in the previous step or given by the with a correspondence file. It will retrieve all the protein sequences of the genomic context of the protein target and taxonomic informations in the INSDC file. These taxonomic informations can be display on the final network. All protein sequences are writen in a fasta file and json format file.

### ClusteringIntoFamilies step

The fasta file is given to MMseqs2 which will group proteins into family with the parameters given by the user. Proteins belonging to a same MMseqs2 family are considered as homologous. This homologous relation is used by the synteny finder step

  - The file of protein sequences in fasta format.
    ```
    >828 // protein unique identifier from proteins_parsingStep.json
    MNDQLFKKVLGYIESESYLMAYRELHKLADEYMPLATRMDFDALHSSLSIIIGERSGYPDIADQLADTAGFYERLAYLLTKKLLGDDEAGEKADTLMLCVVAFGNHRRN
    ```

  - The file of protein data in json format.
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

At the end of this step a new json file of protein data is created where the family identifier is added (see example below)

- The file of protein data in json format.
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


### Synteny Finder step

At this step, NetSyn compute synteny between each pair of protein target. The definition of a synteny between genomic contexts of two target proteins is computed by the exact graph-theoretical approach which has been described by Boyer et al 2015. At this step the  This step create 2 files :

 - The file of node data in json format.
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

- The file of edge data in json format.
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




### DataExport step

This step export the data on the network. To the 3 previous file, NetSyn add a file with the organism data. 

  - The file of organism data in json format.
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

## Web Requests

UniProt: allows to recover the protein accession and nucleic accession from a UniProt accession (into the GetINSDCFiles part).

EBI-ENA: allows to recover the INSDC file (embl format) from a nucleic accession (into the GetINSDCFiles part).

NCBI-taxonomy: allows to recover the lineage taxonomic from a taxon identifier (into the ParseINSDCFles_GetTaxonomy part).

##. CONTRIBUTORS
Celine CHEVALIER
Jordan LANGLOIS
Mark STAM

