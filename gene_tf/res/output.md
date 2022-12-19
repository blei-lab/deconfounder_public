
We summarize `linear_rmse_sing` and `logistic_rmse_sing`, the RMSEs
from regression against one SNP at a time, for the low SNR setting
with n=5000 and p=100000. (They are often substantially smaller than
`linear_rmse` and `logistic_rmse`, which are the RMSEs from regression
against all SNPs simultaneously.) `res/*_avg.csv` contains all the raw
output files for these results.

RMSEs in `linear_rmse_sing` calculates the RMSE for each SNP (i.e. the
absolute error of the coefficient for this SNP) and then average over
all SNPs, so it calculates $$1/m \sum_{j=1}^m \sqrt{(\beta_j -
\beta_j_hat)^2} = 1/m \sum_{j=1}^m |\beta_j - \beta_j_hat|$$ , which
is the mean absolute error (MAE) of the algorithms as we usually call
it. RMSEs in `linear_rmse` calculates the RMSE for all SNPs jointly
and do not average over all SNPs, so it calculates $$\sqrt{ 1/m
\sum_{j=1}^m (\beta_j - \beta_j_hat)^2 }$$, the root mean squared error
(RMSE) of the algorithms as we usually call it.

### Balding-Nichols Model

||Real-valued outcome RMSE\*10^{-2}|Binary-valued outcome RMSE\*10^{-2}|
|:----|:----|:----|
|No control|4.76|3.64|
|oracle|4.74|3.56|
|PCA|4.73|3.55|
|PPCA|4.66|3.41|
|PF|4.30|3.50|
|LFA|3.68|3.02|
|GMM|4.76|3.63|
|DEF|1.33|2.62|


### PSD Model (alpha=0.01)

||Real-valued outcome RMSE\*10^{-2}|Binary-valued outcome RMSE\*10^{-2}|
|:----|:----|:----|
|No control|3.54|3.18|
|oracle|3.52|3.13|
|PCA|3.52|3.12|
|PPCA|3.47|3.03|
|PF|3.24|3.09|
|LFA|2.22|2.81|
|GMM|3.54|3.18|
|DEF|1.31|2.56|

### PSD Model (alpha=0.1)

||Real-valued outcome RMSE\*10^{-2}|Binary-valued outcome RMSE\*10^{-2}|
|:----|:----|:----|
|No control|2.76|2.97|
|oracle|2.75|2.95|
|PCA|2.75|2.95|
|PPCA|2.72|2.90|
|PF|2.59|2.93|
|LFA|1.93|2.76|
|GMM|2.76|2.97|
|DEF|1.30|2.58|


### PSD Model (alpha=0.5)

||Real-valued outcome RMSE\*10^{-2}|Binary-valued outcome RMSE\*10^{-2}|
|:----|:----|:----|
|No control|2.67|3.00|
|oracle|2.66|2.97|
|PCA|2.66|2.97|
|PPCA|2.63|2.93
|PF|2.49|2.94|
|LFA|1.92|2.78|
|GMM|2.67|3.00|
|DEF|1.32|2.60|



### PSD Model (alpha=1.0)

||Real-valued outcome RMSE\*10^{-2}|Binary-valued outcome RMSE\*10^{-2}|
|:----|:----|:----|
|No control|2.73|3.05|
|oracle|2.72|3.03|
|PCA|2.72|3.03|
|PPCA|2.71|3.00|
|PF|2.57|3.01|
|LFA|1.94|2.83|
|GMM|2.73|3.05|
|DEF|1.35|2.66|


### Spatial Model (a=0.1)

||Real-valued outcome RMSE\*10^{-2}|Binary-valued outcome RMSE\*10^{-2}|
|:----|:----|:----|
|No control|8.05|4.51|
|oracle|7.99|4.34|
|PCA|7.99|4.34|
|PPCA|7.84|4.06|
|PF|4.32|3.45|
|LFA|4.03|3.36|
|GMM|8.05|4.51|
|DEF|1.32|2.57|



### Spatial Model (a=0.25)

||Real-valued outcome RMSE\*10^{-2}|Binary-valued outcome RMSE\*10^{-2}|
|:----|:----|:----|
|No control|3.83|3.20|
|oracle|3.81|3.17|
|PCA|3.80|3.15|
|PPCA|3.76|3.07|
|PF|2.50|2.86|
|LFA|2.22|2.80|
|DEF|1.29|2.54|


### Spatial Model (a=0.50)

||Real-valued outcome RMSE\*10^{-2}|Binary-valued outcome RMSE\*10^{-2}|
|:----|:----|:----|
|No control|5.25|3.59|
|oracle|5.21|3.49|
|PCA|5.21|3.49|
|PPCA|5.18|3.43|
|PF|2.91|2.94|
|LFA|2.70|2.88|
|GMM|5.25|3.59|
|DEF|1.25|2.47|


### Spatial Model (a=1)

||Real-valued outcome RMSE\*10^{-2}|Binary-valued outcome RMSE\*10^{-2}|
|:----|:----|:----|
|No control|4.49|3.39|
|oracle|4.46|3.31|
|PCA|4.46|3.31|
|PPCA|4.43|3.26|
|PF|2.59|2.87|
|LFA|2.41|2.83|
|GMM|4.49|3.38|
|DEF|1.27|2.51|


