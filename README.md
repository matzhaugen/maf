# Maf 

Maximum autocorrelation factors (MAF) are are linear combinations of a set of concurrent multivariate time series that contain the maximum amount of autocorrelation. It can be shown that if you have a data set where each time series contains a couple of smooth underlying signals corrupted by noise, that MAFs achieve the maximum signal-to-noise ratio and highest correlation with the true underlying signals.

# Install

```
library(devtools)
install_github("matzhaugen/maf")
```

# Run the example

```
maf.object = maf(treeringTimeseries)
plot(maf.object)
```