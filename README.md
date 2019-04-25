# single-cell-RNAseq-HSV1
Analysis files accompanying eLife paper entitle "Single cell analysis of HSV-1 infection reveals anti-viral and developmental programs are activated in distinct sub-populations with opposite outcomes"

This repository contains the Matlab code used to analyze and generate the images for the single-cell anlyses in the paper.
Files contained are:
1. scRNAseq_analysis_final_042519_for_Github.m : this is a matlab code file that will run all the analyses and procude the figures from the paper. It was constructed and used in MATLAB2018a version.
2. Files_for_analysis_intitation.mat : this matlab files contains the raw data from sequencing of the three samples (mock_orig, wt_orig ad do_orig), the gene names (names_mock, names_wt and names_d0) and "colors" - used for coloring the plots. the mock_orig, wt_orig and d0_orig are DGE (digital gene expression) matrices. each row is a gene and is coloum is a cell. the names_mock, names_wt and names_d0 are arrays holding the names of the genes (=rows) for each corresponding matrix. mock are mock-infected cells, wt are wild-type HSV-1 infected cells and d0 are deltaICP0 HSV-1 infected cells.
3. a folder named distributionPlot. This is a function for creating violin plots. Please note it was authored by Jonas and was obatained from: https://www.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions-distributionplot-m. 
4. a folder named kakearney-boundedline-pkg-8179f9a. This is a function to create shaded lines. Please note is was authored by Kelly Kearney and was obtained from https://www.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m.

Major steps in the analysis:
1. Pre-processing: Fileting out cells with low and high UMI and high mitochondrial gene content.
2. Analysis of viral genes heterogeneity: analyses the distrubiton of total HSV-1 genes among indidual cells, as well as the distributuion of the four viral gene classes (IE - immediate-early, E - early, L1 - late1 and L2 - late2).
3. Analysis of differnitally expressed host genes: regression-out of UMI and cell-cycle effects using a linear model. Differential gene identification by Wilcoxson ranksum test, followed by Benjamini and Hochberg false detection rate (FDR) correction for multiple comparisions.
4. Analysis of ISGs (Interferon stimulated genes) in cells expressing high vs. low levels of viral genes

Please contact the first author of the paper, Nir Drayman (nirdra@uchicago.edu), with any comments, concerns or requests. We'll be happy to help you adapt our pipeline for your analyses!
