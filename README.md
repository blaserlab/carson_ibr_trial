# carson_brooke
Analysis of mdscs in cancer patients treated with ibrutinib

Initial data report uploaded 4/9/2020.  

An easy way to provide comments is to click on issues at the top of the page.  Another way is to click on "N commits" (N = 2 currently) just above the blue line, then click on the most recent commit, then scroll down and add a comment.

--------------------------------------------
Update April 11 2020

1. The top markers function was used on the mdsc population to identify all of the genes that are significantly different between c1 dminus7 and c1 dplus1.  The raw output is in "data_out/TM_mdsc_cycle_day.csv".  This file has to be interpreted carefully.  Monocle uses several metrics to select these top markers.  First it calculates the mean expression for all genes in all of the cells in the dataset.  From this it calculates the Jensen-Shannon divergence for all genes.  This is a metric to determine whether the data distributions in the two groups of cells are different.  It then selects the top N markers (I chose 50) based on Jensen-Shannon divergence and performs logistic regression to determine the association of each gene with a group of cells.  This test is more computationally intensive so it only performed on the subset of genes selected by Jensen-Shannon. From this it calculates a p value using the likelihood ratio test and a q value using Bonferroni correction.  

    This is great, but still it is not in a useful format if you want to know which genes are up or down in your cells of interest.  So what we can also do is select the genes that are significant with q value < 0.05 and plot their expression in both groups.  These are in the new folder "plots_out/mdsc_single_violins/".  In these plots, the p value is recalculated using the wilcox test.  This is the best choice because it is non-parametric and the values we are using in the plots are <img src="https://render.githubusercontent.com/render/math?math=\log_{10}(expression + 1)>.  This is a standard trick:  adding a pseudocount of 1 to all of the values allows you to plot the zeros on a log scale since $\log_{10}0 = -\infty$ and therefore $\log_{10}(0+pseudocount) = 0$.

    I gave you all of the genes with q<0.05 in this folder.  MT-genes are mitochondrial gene RNA and should be ignored.  Probably also actin and gapdh even if they look significant.  Most of the other stuff looks interesting.  In addition to the plots, for each gene there is a small text file with a stat report for that gene.

2.  "data_out/TM_TNK_reclust.csv" and "plots_out/TNK_violin.pdf" were calculated in the same way as I just described.  However the violin plots shown are the top 4 in each group by q value.  I made single violin plots and stat reports for genes with q<0.05 like I did for MDSC ("plots_out/TNK_single_violins/").

3.  I don't think we should compare MDSC by response, since it is so heavily skewed to pt 11.  But T cells by response is probably ok.  I made a new plot with the T cells binned by Tc/NK or Mem/Naive phenotype and response: "plots_out/TNK_by_response.pdf".  I also generated the stats ("data_out/TM_TNK_by_response.csv") and the single plots and stat reports as before ("plots_out/TNK_by_response_single_violins/").

4.  To determine the TCR diversity in each sample, I calculated the Shannon diversity index as $$-1*\sum_{i=1}^{R}p_{i}\ln p_{i}$$

    where $p_{i}$ is the proportional abundance of the $i^{th}$ clonotype within the sample and $R$ is richness or total number of observed clonotypes in the sample.  This is one of the standard ways to quantify T cell diversity.  It is basically a weighted average based on the abundance of each T cell clone.  
    
    It doesn't look like there is a consistent difference here between pre- and post-treatment.  See "data_out/clonotype_shannon.csv".  Patients 11 and 22 go down, 15 goes up and 17 stays about the same.  
    
Happy to discuss further.  If there are other analyses you want, please let me know.  If there are plots you like I can make cosmetic adjustments to fit the rest of your figures at the last minute prior to submission. 

Thanks!
