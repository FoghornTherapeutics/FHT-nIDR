# Motivations

The Irreproducible Discovery Rate (IDR) is a method to combine combining replicates to only get the highly reproducible peaks and consistency between replicates in high-throughput experiments. It is used to measure the consistency of results (such as the identification of transcription factor binding sites, histone marks, or gene expression levels) across different experimental replicates. The standard method involves comparing the rank of results across different replicates to estimate the proportion of findings that are reproducible (consistent across replicates) versus those that are irreproducible (inconsistent or likely to be noise).  One way to assess concordance of peak calls between replicates is to implement a statistical procedure. A popular method is the [IDR framework developed by Qunhua Li and Peter Bickel’s group](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-5/issue-3/Measuring-reproducibility-of-high-throughput-experiments/10.1214/11-AOAS466.full). It compares a pair of ranked lists of regions/peaks and assigns values that reflect its reproducibility. The core IDR R package can be downloaded from the [IDR download page](http://cran.r-project.org/web/packages/idr/index.html).

This method can only compare a pair of replicates. When you have n>2 replicates, one method would be to use the above to do pairwise comparisons cand choose the best pair.  However, in that case you have to compare $\binom{N}{2}$ pair of replicates. If N=4, you then have to compare by hand 6 pair of replicates.  In addition you are then losing the statistical power / signal that comes from having greater numbers of replicates. We propose another method to compare all the replicates at once and give the best set of peaks that are consistent accross all of them. We call it nIDR where n indicates it is for arbitrary number of replicates.


# Mathematical derivation of the new algorithm

Starting from equation (2.1) in Li et al:

$$\psi(t) = \frac{1}{N} \sum_{i=1}^{N}{1 \left( X1_i>x_{X1}(t), X2_i>x_{X2}(t) \right) } $$

* $N$ number of features i.e. peaks
* $1(...)$ is a function that returns $1$ if all the arguments are true otherwise returns 0
* $X1_i$ values observed in replicate X1 for $i$ th feature
* $X2_i$ values observed in replicate X2 for $i$ th feature
* $t$ is the fractional rank i.e. number from 0 to 1
* $x_{X1}(t)$ percentile value of replicate X1.  i.e. $x_{X1}(0.5)$ is the 50th percentile aka median of the data for replicate X1
* $x_{X2}(t)$ percentile value of replicate X2.  i.e. $x_{X2}(0.5)$ is the 50th percentile aka median of the data for replicate X2

To help us understand how this works we use this illustrative example for $t=0.5$ -  define:

* $median_{X1} = x_{X1}(t=0.5)$ the median of the data measured for replicate X1
* $median_{X2} = x_{X2}(t=0.5)$ the median of the data measured for replicate X2

then:

$$\psi(t=0.5) = \frac{1}{N} \sum_{i=1}^{N}{1 \left( X1_i>median_{X1},X2_i>median_{X2} \right) } $$

Restating the above:  the operand of the sum $1(...)$ is 1 if for feature $i$ the value for replicate $X1_i > median_{X1}$ and the value for replicate $X2_i > median_{X2}$, otherwise it is 0.  Therefore, $\psi(t=0.5)$ counts up the features features where both X1 and X2 have values greater than their respective medians.

With the above understanding, we can do some mathematical rearrangements to develop an algorithm that allows us to more easily extend this to multiple replicates.  Note that Li et al also address this in section 4 of their supplemental "Extension of our model to the case of m > 2"; our goal here is rearrange for our understanding and to make reasonably efficient computer code.

define fractional rank of X1 - $frac\textunderscore rank\textunderscore X1$ as a vector where each entry is the fractional rank of the corresponding $i$ th value in X1

$$frac\textunderscore rank\textunderscore X1_i = \frac{index\textunderscore of\textunderscore sorted\textunderscore X1(i)}{N}$$

then this is true:

$$ X1_i>x_{X1}(t) == frac\textunderscore rank\textunderscore X1_i>t $$

then the equation for $\psi$ can be rewritten as:

$$\psi(t) = \frac{1}{N} \sum_{i=1}^{N}{1(frac\textunderscore rank\textunderscore X1_i>t, frac\textunderscore rank\textunderscore X2_i>t)} $$

the sum's operand $1(...)$ can be rewritten as:
  
$$1\Bigl(min[frac\textunderscore rank\textunderscore X1, frac\textunderscore rank\textunderscore X2]>t\Bigr)$$

yielding:

$$\psi(t) = \frac{1}{N} \sum_{i=1}^{N}{1 \Bigl(min(frac\textunderscore rank\textunderscore X1_i,frac\textunderscore rank\textunderscore X2_i)>t \Bigr)} $$

Instead of calcuating the sum for each value of $t$ we can generate an ECDF (empirical cumulative distribution).

First, create the vector of values $min\textunderscore rank$ such that:

$$min\textunderscore rank_i = min(frac\textunderscore rank\textunderscore X1_i, frac\textunderscore rank\textunderscore X2_i))$$

Sort the above to create $sort\textunderscore min\textunderscore rank$ such that:

$$sort\textunderscore min\textunderscore rank_i <= sort\textunderscore min\textunderscore rank_{i+1}$$

starting at beginning of vector $sort\textunderscore min\textunderscore rank$, find $sort\textunderscore min\textunderscore rank\textunderscore index(t)$ which is first occurrence of $t$ in $sort\textunderscore min\textunderscore rank$

$$sort\textunderscore min\textunderscore rank\bigl[sort\textunderscore min\textunderscore rank\textunderscore index(t)\bigr] = t$$

Then:

$$\psi(t) = \frac{sort\textunderscore min\textunderscore rank\textunderscore index(t)}{N}$$

The above equation then leads to a fairly straightforward algorithm to implement:
1. calculate the fractional rank of values within each replicate
2. For each feature in the replicates, calculate the minimum fractional rank
3. calculate an ECDF of these minimum fractional ranks - this is now $psi(t)$
This methodology has allowed us to extend IDR to arbitrary numbers of replicates, and to extend the number of features into the hundreds of millions (and we expect billions).

# FHT-nIDR data flow


To illustrate the mechanism, we take a group of 3 replicates: A1, A2 and A3.  

**Step 1**: Combine and merge three replicates peak id, peak origin and logFC on the same narrowPeak file. For this step, we use bedtools merge to merge the peaks. In consequence, we can have several logFC from the same replicate. <br/>
**Step 2**: Reformat bed file with a logFC per replicate and by peak id. If there is logFC from a replicate, the value is 0. When there are several logFC from the same replicates, we compute the average between them. <br/>
**Step 3**: Generate a shuffled distribution of the logFC over the given list of peak ID for each replicate.<br/>
**Step 4**: Compute the min percent rank of the shuffled distribution for each peak over the three replicates. In other words, we compute the percent rank of each replicates, i.e., we have three percent rank for each replicate. Then we only keep the minimum of the percent rank for each peak id. <br/>
**Step 5**: Compute the ECDF and select the min percentage rank that corresponds to keeping above 90% of the reads.  In this example, the min percent rank is around 0.53. <br/>
**Step 6**: Filter all peak id that have a min percent rank lower than the selected threshold from the null distribution. <br/>
**Step 7**: One of the nIDR outputs is the Empirical Cumulative Distribution Function (ECDF) plots of the consistency across replicates of all peaks found in a group of replicates. The x-axis is the "min percent rank" which indicates the consistency of a peak across the replicates, a higher value corresponds to a higher consistency of peaks accross replicates. The y-axis is the fraction peaks that have at least that value. The red curve is plotting the actual data and the blue curve is the simulated null distribution. The green dashed line indicates a p-value of 0.1 based on the blue null curve and determines the consistency score threshold to use for keeping the real peaks. In this case, the green dashed line indicates the threshold where 90% of the null peaks have a consistency score below 0.53. Therefore this sets a threshold for choosing peaks with a p value < 0.1. 0.53 is considered to be a relatively high min percent rank. The null distribution is above the true values, meaning that replicates of the same group show more consistency and they are not random.


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

<img src="readme_figures/HCT116_nIDR.JPG" alt="image" style="width:600px;height:auto;">
<img src="readme_figures/TOV21G_nIDR.JPG" alt="image" style="width:600px;height:auto;">


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

<img src="readme_figures/Example2_ECDF.JPG" alt="image" style="width:600px;height:auto;">

When comparing the [Differential Peak Area (DPA) ](https://github.com/FoghornTherapeutics/FHT-ATACseq-pipeline/blob/main/README.md#2-differential-peak-area-dpa) using the regular narrowPeak IDR (rIDR) and our method nIDR, the DPA logFC are overall well correlated. This scatter plot compares the rIDR on the x-axis and the nIDR on the y-axis for each contrast (DMSO vs treatment) at 24h and 72h. The points on the x-axis or on the y-axis are ?? XXXXXXXXXXXXXXXX

<img src="readme_figures/Example2_logFC_scatter_plot.JPG" alt="image" style="width:600px;height:auto;">


Finally, we compare the overlap of peaks between the nIDR and the individual samples for the replicates at 24h. The first one is for DMSO and the second is for treatment. Most of the peaks are ?? XXXXXXXXXXXXXXXX. The second row compares the  nIDR narrowPeak and each of the pairwise regular IDR narrowPeak. Less that half of the peaks added in nDR are not part of the pairwise rIDR narrowPeak.

<img src="readme_figures/Example2_VenDiagram.JPG" alt="image" style="width:600px;height:auto;">






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

First, we compare the standard output from rIDR, we can identify that A2 and A3 show more consistency than A1. In the same way, B1 and B2 have more overlapping peaks overall.

<img src="readme_figures/Example3_rIDR.JPG" alt="image" style="width:600px;height:auto;">


Ignoring the fact the A1 and B3 are outliers, we apply the nIDR method on all three replicates in both groups. 

Visually, the min rank corresponding to 10% of kept reads computes a threshold that is not where the null distribution is above the real distribution. Quantitatively, the computed threshold is really low at about 0.35 for each group and should be above 0.5. <br/>

We then apply the method removing these replicates on the second row. The true distribution is still very close to the random data. However, the min percent rank is much higher now around 0.67.


<img src="readme_figures/Example3_ECDF.JPG" alt="image" style="width:600px;height:auto;">


Finally, when removing the two replicates A1 and B3, the scatter plot of the DPA logFC shows a stronger correlation between the two  

<img src="readme_figures/Example3_logFC_scatter_plot.JPG" alt="image" style="width:600px;height:auto;">





