# TODO

## Directory content

There are two datasets:
1. RNA-Seq
2. Cell wall measurements

Cell-wall-measurement.tsv is a table of the average thickness (µm) of the algal cell wall at the different experimental time points (and for the control + cold stress conditions). The algae were grown in control (20°C) and cold stress (5°C) conditions. 

## Aim

Look at the changes in the overall transcriptome over time, with the intention of identifying candidate genes that might control/explain the algae’s resistance to cold stress. I measured the cell wall of the algae at the different time points and we see that the cell wall becomes thicker when the algae is under cold stress as compared to the control (this starts to happen after 24hours it seems). We want to match these measurements with the transcriptomic data in order to understand what happens to the cell wall when the algae is exposed to cold stress.

## Tasks

### Results

1. PCA on vst,blind=FALSE
    1. Probably done in BioQA - check src/R (use the right sample file)
    2. Redo with vst aware data
    3. ellipse?
    4. Rm S/N (independent filtering) - DONE (unnecessary as filtering by taxonomy reduced the number of genes to ~50k)
2. Heatmap
    1. Order by time
    2. Order genes
3. DE more/ less cold exposure, GO, KEGG
    1. Venn Diagram (acute (1) early(4-12), late(24-120))
    1. Check what was done in src/R
    2. Identify cold response genes, but be wary of cell cycle!!! - cold acclimation; response to cold; cellular response to cold (also Venn)
    3. plot that cold activation is time sensitive
    4. cell wall genes 4 vs. 0 and 120 vs. 0 (volcano plot)
4. Generate upset plot (probably exists as Venns) - no we want Venn
5. Lipid / carbs (also Venn)
    1. As for sucrose synthase (that is done, check for more? Ask Olivia for lipids)
6. Mfuzz on cold (done, what figure?)
7. cold shock proteins (probably check the annotation / literature?)
8. cold response genes (probably check the annotation / literature?) that are irresponsive in that dataset
9. cell wall genes (done)
10. pathways with enzymes involved in cell wall thickening (ask Olivia for pathways?)
11. cell wall synthesis pathway (as for lignin)
