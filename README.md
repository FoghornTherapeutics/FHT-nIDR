# Motivations

Measures consistency between replicates in high-throughput experiments. Also uses reproducibility in score rankings between peaks in each replicate to determine an optimal cutoff for significance. The core IDR R package can be downloaded from the IDR download page: http://cran.r-project.org/web/packages/idr/index.html



# FHT-nIDR

One way to assess concordance of peak calls between replicates is to implement a statistical procedure. 
A popular method is the IDR framework developed by Qunhua Li and Peter Bickel’s group. 
It compares a pair of ranked lists of regions/peaks and assigns values that reflect its reproducibility.
