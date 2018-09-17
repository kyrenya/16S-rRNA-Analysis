In this repository you will find:

  Folder containing: 
  
      - RMarkdown files: 
      
            Rmd: 16_Q1_Q2_Dada2.Rmd
            
            Interactive Report: 
            
            html: 16_Q1_Q2_Dada2.html  
            
      - data-set: 16SDATA folder
      
      - three folders containing the files corresponding to the analysis through: 
      
          - QIIMEI
          
          - QIIME2
          
          - DADA2
          
      - folders needed for further analysis of the results: 
      
          - COMPARING DM
          
          - OTU_SUMMARIES
          
      - cache folders for faster knitting of the RMarkdown file
      
          - 16S_Q1_Q2_Dada2_cache
          
          - 16S_Q1_Q2_Dada2_files
          
      - R environment: workspace.RData
      
      

![alt text](https://raw.githubusercontent.com/kyrenya/16S-rRNA-Analysis/master/workflow.jpg)

QIIME, QIIME2 and DADA2 workflows overview. Simplified scheme of processes (blue boxes) and objects (blue ellipses) yielded in each pipeline. QIIME1 and QIIME2 allow the complete analysis of sequencing data, providing tools for carrying out statistics and plotting and visualizing data. On the contrary, DADA2 package pipeline ends in the production of the feature table, which shall be used with other tools to complete the analysis. Adapted from Navas-Molina et al., 2013 and https://qiime2.org



Tools needed for Running the Pipelines in the 16_Q1_Q2_Dada2.Rmd: 
  Follow the instructions for installing QIIME using miniconda:

QIIMEI
http://qiime.org/install/install.html


Dependency versions
===================
          QIIME library version:	1.9.1
           QIIME script version:	1.9.1
    qiime-default-reference version:	0.1.3
                  NumPy version:	1.15.1
                  SciPy version:	1.1.0
                 pandas version:	0.23.1
             matplotlib version:	1.4.3
            biom-format version:	2.1.6
                   h5py version:	2.8.0 (HDF5 version: 1.10.2)
                   qcli version:	0.1.1
                   pyqi version:	0.3.2
             scikit-bio version:	0.2.3
                 PyNAST version:	1.2.2
                Emperor version:	0.9.51
                burrito version:	0.9.1
       burrito-fillings version:	0.1.1
              sortmerna version:	SortMeRNA version 2.0, 29/11/2014
              sumaclust version:	SUMACLUST Version 1.0.00
                  swarm version:	Swarm 1.2.19 [Mar  1 2016 23:41:10]
                          gdata:	Installed.


QIIME2
https://docs.qiime2.org/2018.8/install/native/#install-qiime-2-within-a-conda-environment
(version 2018.8.0)


DADA2 R PACKAGE
https://bioconductor.org/packages/release/bioc/html/dada2.html


The R packages versions used for this work are listed in the Session_Info section in the report 16S_Q1_Q2_Dada2.html



[Click here to pre-view the 16SAnalysis Report document (16S_Q1_Q2_Dada2.html)]
(http://htmlpreview.github.io/?https://github.com/kyrenya/master/16S-rRNA-Analysis/16SAnalysis/16S_Q1_Q2_Dada2.html 
http://htmlpreview.github.io/?https://github.com/kyrenya/16S-rRNA-Analysis/blob/master/16SAnalysis/16S_Q1_Q2_Dada2.html
