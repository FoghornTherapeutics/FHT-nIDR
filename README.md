# Motivations

The Irreproducible Discovery Rate (IDR) is a method to combine combining replicates to only get the highly reproducible peaks and consistency between replicates in high-throughput experiments. It is used to measure the consistency of results (such as the identification of transcription factor binding sites, histone marks, or gene expression levels) across different experimental replicates. The standard method involves comparing the rank of results across different replicates to estimate the proportion of findings that are reproducible (consistent across replicates) versus those that are irreproducible (inconsistent or likely to be noise).  One way to assess concordance of peak calls between replicates is to implement a statistical procedure. A popular method is the [IDR framework developed by Qunhua Li and Peter Bickel’s group](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-5/issue-3/Measuring-reproducibility-of-high-throughput-experiments/10.1214/11-AOAS466.full). It compares a pair of ranked lists of regions/peaks and assigns values that reflect its reproducibility. The core IDR R package can be downloaded from the [IDR download page](http://cran.r-project.org/web/packages/idr/index.html).

This method can only compare a pair of replicates. When you have n replicates, you have to compare $\binom{N}{2}$ pair of replicates. If N=4, you then have to compare by hand 6 pair of replicates. We propose another method to compare all the replicates at once and give the best set of peaks that are consistent accross them. We call it nIDR where n is the number of replicates.



# FHT-nIDR data flow


To illustrate the mechanism, we take a group of 3 replicates: A1, A2 and A3.  

**Step 1**: Combine and merge three replicates peak id, peak origin and logFC on the same narrowPeak file. For this step, we use bedtools merge to merge the peaks. In consequence, we can have several logFC from the same replicate. <br/>
**Step 2**: Reformat bed file with a logFC per replicate and by peak id. If there is logFC from a replicate, the value is 0. When there are several logFC from the same replicates, we compute the average between them. <br/>
**Step 3**: Generate a shuffled distribution of the logFC over the given list of peak ID for each replicate.<br/>
**Step 4**: Compute the min percent rank of the shuffled distribution for each peak over the three replicates. In other words, we compute the percent rank of each replicates, i.e., we have three percent rank for each replicate. Then we only keep the minimum of the percent rank for each peak id. <br/>
**Step 5**: Compute the ECDF and select the min percentage rank that corresponds to keeping above 90% of the reads.  In this example, the min percent rank is around 0.53. <br/>
**Step 6**: Filter all peak id that have a min percent rank lower than the selected threshold from the null distribution. <br/>
**Step 7**: One of the nIDR outputs is the Empirical Cumulative Distribution Function (ECDF) plots of the consistency across replicates of all peaks found in a group of replicates. The x-axis is the "min percent rank" which indicates the consistency of a peak across the replicates, a higher value corresponds to a higher consistency of peaks accross replicates. The y-axis is the fraction peaks that have at least that value. The red curve is plotting the actual data and the blue curve is the simulated null distribution. The green dashed line indicates a p-value of 0.1 based on the blue null curve and determines the consistency score threshold to use for keeping the real peaks. In this case, the green dashed line indicates the threshold where 90% of the null peaks have a consistency score below 0.53. Therefore this sets a threshold for choosing peaks with a p value < 0.1. 0.53 is considered to be a relatively high min percent rank. The null distribution is above the true values, meaning that replicates of the same group show more consistency and they are not random.


CHANGE ECDF WITH THRESHOLD OF .53!!!!
<img src="readme_figures/data_flow.JPG" alt="image" style="width:1000px;height:auto;">



# Example 1

### Data
The standard pipeline was run on publicly available data from paper "[Chromatin accessibility underlies synthetic lethality of SWI/SNF subunits in ARID1A-mutant cancers](https://elifesciences.org/articles/30506#content)" looking for potential PD markers as well as what an ATAC-seq profile looks like. This paper has ATACseq results of ARID1A-/- cancer cell lines (native or CRISPR knockout) with ARID1B knockdown. 

**Data from GEO series**:  [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101975)

Overview of experiment:
* Biological context (N=2):  TOV21G, HCT116
* wild type and modified with stable ARID1A knockout
* Perturbagens (N=1):  shRNA knockdown of ARID1B
* Negative control (N=1):  wild type / untreated
* Replicates:  N=2

#### HCT116 (ACH-000971)
* WT: SRR5876158 & SRR5876159
* ARID1B knockdown: SRR5876160 & SRR5876161
* ARID1A knockout: SRR5876162 & SRR5876163
* ARID1A knockout ARID1B knockdown:SRR5876164 & SRR5876165

#### TOV21G (ACH-000885)
* WT: SRR5876661 & SRR5876662
* ARID1B knockdown: SRR5876663 & SRR5876664

### Results

In this example, there are only two replicates by group. Therefore, there are only $\binom{N}{2} = \binom{2}{2}= 1$ comparison by group. This data is a good example to first verify that the peaks from the nIDR computation overlap with the standard method. 

First of all, the output from the standard IDR show that the replicates have high consitency betwen groups.
 
<img src="readme_figures/ARID_paper_standard_IDR.JPG" alt="image" style="width:900px;height:auto;">

In the same way, the results from the ECDF when using the nIDR method shows that we have a really high min percent rank around 0.65. 

XXXXXXXXXXXXXXXX ADD OUR nIDR ECDF XXXXXXXXXXXXXXXX

To compare the two IDR narrowPeak files, we study their overlap with a Ven Diagram. And we can see that the standard IDR narrowPeak is almost included in the nIDR narrowPeak. The additional peaks from nIDR could be explained by the fact that we are adding more information ?? XXXXXXXXXXXXXXXX

<img src="readme_figures/ARID_paper_nIDR.JPG" alt="image" style="width:900px;height:auto;">




# Example 2

### Data

In this example, we compare some samples from random data with 3 replicates by group. Therefore, there are only $\binom{N}{2} = \binom{3}{2} = 3$ comparison by group.

Overview of experiment:
* Biological context (N=2): Cell line 1
* Perturbagens (N=1): treatment
* Time point (N=2): 24h and 72h
* Negative control (N=1): wild type / untreated
* Replicates: N=3


First of all, the output from the standard IDR show that the replicates have high consitency betwen groups.

XXXXXXXXXXXXXXXX ADD OUR STANDARD OUPTUT PNG XXXXXXXXXXXXXXXX
for 24h: chosen DMSO rIDR is A2_A3.
for 24h: chosen FHT rIDR is A5_A6.


The ECDF shows the null distribution above the true values with a wide gap, meaning that replicates of the same group show more consistency and they are not random. The green dashed line indicates the threshold where 90% of the null peaks have a consistency score below 0.53 which is considered to be a relatively high min percent rank.

XXXXXXXXXXXXXXXX ADD OUR nIDR ECDF XXXXXXXXXXXXXXXX

When comparing the [Differential Peak Area (DPA) ](https://github.com/FoghornTherapeutics/FHT-ATACseq-pipeline/blob/main/README.md#2-differential-peak-area-dpa) using the regular narrowPeak IDR (rIDR) and our method nIDR, the DPA logFC are overall well correlated. This scatter plot compares the rIDR on the x-axis and the nIDR on the y-axis for each contrast (DMSO vs treatment) at 24h and 72h. The points on the x-axis or on the y-axis are ?? XXXXXXXXXXXXXXXX

XXXXXXXXXXXXXXXX ADD SCATTER PLOT logFC XXXXXXXXXXXXXXXX


Finally, we compare the overlap of peaks between the nIDR and the individual samples for the replicates at 24h. The first one is for DMSO and the second is for treatment. Most of the peaks are ?? XXXXXXXXXXXXXXXX. The second row compares the  nIDR narrowPeak and each of the pairwise regular IDR narrowPeak. Less that half of the peaks added in nDR are not part of the pairwise rIDR narrowPeak.

XXXXXXXXXXXXXXXX ADD Ven DIAGRAM XXXXXXXXXXXXXXXX


add NS-22.0044:
- nIDR ECDF
- normal IDR png
- scatter plot of logFC
- Ven diagram






# Example 3: Groups with outlier(s)

### Data

One of the challenge of this method is to compare all N replicates of a group regardless of their consistency. However, in the ATACseq pipeline described [here](https://github.com/FoghornTherapeutics/FHT-ATACseq-pipeline/blob/main/README.md), we first compute a series of QC measures (multiqc, FRiP, PCA, sample correlations) to identify any potential outliers in a group of replicates. You can find an example of a QC result that identifies an outlier [here](https://github.com/FoghornTherapeutics/FHT-ATACseq-pipeline/blob/main/QC_example_with_outlier.md). <br/>

This next example compares the nIDR narrowPeak with all 3 replicates of a group with and then without the outlier in the treated group and the negative control group. 

Overview of experiment:
* Biological context (N=1): Cell line 1
* Perturbagens (N=1): treatment
* Time point (N=1): 24h 
* Negative control (N=1): wild type / untreated
* Replicates: N=3


Samples A1, A2 and A3 are for the negative control. A1 was identified as an outlier (low FRiP scores,  did not cluster with any other replicates in PCA plot and lower sample-to-sample correlation).
Samples B1, B2 and B3 are for the negative control. B3 was also identified as an outlier (lower insret size, did not cluster with any other replicates in PCA plot and lower sample-to-sample correlation).


#### Negative control group:

First looking at the standard IDR output, we can see that A2 and A3 show more consistency than A1. In the same way, B1 and B2 have more overlapping peaks.

Even if we decide to ignore the fact that A1 and B3 are outliers from our QC measures and the standard IDR ouptut, once we compute the nIDR narrowPeak, the ECDF plot raises an extra flag. 
Visually, the min rank corresponding to 10% of kept reads computes a threshold that is not where the null distribution is above the real distribution. Quantitatively, the computed threshold is really low at about 0.35 for each group and should be above 0.5. 

XXXXXXXXXXXXXXXX ADD OUR STANDARD OUPTUT PNG XXXXXXXXXXXXXXXX

XXXXXXXXXXXXXXXX ADD OUR nIDR ECDF XXXXXXXXXXXXXXXX

XXXXXXXXXXXXXXXX ADD SCATTER PLOT logFC XXXXXXXXXXXXXXXX




add NS-22.0049:
- nIDR ECDF
- normal IDR png
- scatter plot of logFC
- Ven diagram



