# sell-gene-explore

To quickly look at the transcriptomic expression of a gene across different external datasets and internally generated data.

Make sure you have access to Sellgren lab's KI server before running the script.

## Datasets Used (external)

1. **BrainSpan Data:** 8 pcw - 40 y bulk RNAseq [[www.brainspan.org](http://www.brainspan.org)] Spatio-temporal transcriptome of the human brain. Nature. 2011 Oct 26; 478(7370):483-9. PMID: 22031440
2. **Human brain 1st trimester:** Braun et al. DOI: 10.1126/science.adf1226. Epub 2023 Oct 13. PMID: 37824650
3. **Human brain (2nd trimester-Adult):** Velmeshev et al. DOI: 10.1126/science.adf0834. Epub 2023 Oct 13. PMID: 37824647; PMCID: PMC11005279

## Running the Script
You will first need to access the Sellgren lab folder at least once before starting. If you are working from home or remotely, make sure your VPN is connected and you have opened the server at least once. 

### Installing R (skip if already installed)
To run certain analyses, you need to have R installed. Download and install R from:
- [CRAN (Comprehensive R Archive Network)](https://cran.r-project.org/)
- Follow the installation instructions for your operating system.

#### macOS Installation
1. Download the latest R package for macOS from the CRAN website.
2. Open the downloaded `.pkg` file and follow the on-screen instructions to install R.

#### Windows Installation
1. Download the latest R installer for Windows from the CRAN website.
2. Run the `.exe` installer and follow the setup instructions.


### Start Here: Cloning the Repository
Open your Terminal on Mac (or command line in windows), copy-paste these lines in your terminal:
```bash
git clone https://github.com/SellgrenLab/sell-gene-explore/
cd sell-gene-explore
```

To analyze gene expression, run (copy-paste) the following command in your command terminal:

```bash
Rscript GeneExplorer.r --gene <GENE_NAME> 
```
Replace ```<GENE_NAME>``` with the gene symbol you want to analyze. Check for gene aliases. 

The default color palette is set to YlGnBu. If you want to change it, you can specific any palette code from the Rcolorbrewer palettes with a -p. Tip: Select either a sequential or diverging color palette. 

```bash
Rscript GeneExplorer.r --gene <GENE_NAME> -p <color_palette>
```
Replace ```<color_palette>``` with the palette. For example RdBu, YlOrRd, Reds, Blues, PuBuGn, RdYlGn, RdYlBu, Spectral, RdGy or PrGn, etc. 

## Output

Three files with figures corresponding to the above datasets are generated in the ```results/``` folder.

### Happy Exploring

##### Note 
The code is adapted to specific files that are present in my (SM) folder on the server. If the path is modified, the script needs to be updated accordingly. The external datasets are used as published by the authors and not modified. If there are issues/errors while running or if you have any questions, feel free to get in touch with me. 





