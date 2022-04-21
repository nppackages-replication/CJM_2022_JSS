*-------------------------------------------------------------------------------
* lpdensity: Local Polynomial Density Estimation and Inference
* Matias D. Cattaneo, Michael Jansson, and Xinwei Ma
* Replication Code
*-------------------------------------------------------------------------------
* hlp2winpdf, cdn(lpdensity) replace
* hlp2winpdf, cdn(lpbwdensity) replace
*-------------------------------------------------------------------------------
* net install lpdensity, from(https://raw.githubusercontent.com/nppackages/lpdensity/master/stata) replace
*-------------------------------------------------------------------------------
clear all
set more off
mata: mata mlib index
set scheme sj

*-------------------------------------------------------------------------------
* Generate data
*-------------------------------------------------------------------------------
set seed 42
set obs 4000
gen data = rnormal(1, 1)
drop if data <= 0 // normal(1, 1) truncated from below at 0
drop if _n > 2000 // keep only 2000 observations
sum data

*-------------------------------------------------------------------------------
* lpdensity(): Estimation with bandwidth 0.5 on provided grid points
*-------------------------------------------------------------------------------
// generate grid 
gen grid = -0.5 + 0.5 * _n if _n <= 9

lpdensity data, grid(grid) bw(0.5)
lpdensity data, grid(grid) bw(0.5) level(99) // 99% CI

*-------------------------------------------------------------------------------
* lpdensity(): Estimation with bandwidth 0.5 on quantile grid points
*-------------------------------------------------------------------------------
lpdensity data, bw(0.5)
lpdensity data, bw(0.5) q(2) // conventional CI

*-------------------------------------------------------------------------------
* lpdensity(): Plotting
*-------------------------------------------------------------------------------
capture drop grid
gen grid = -0.05 + 0.05 * _n if _n <= 81

lpdensity data, grid(grid) bw(0.5) plot

*-------------------------------------------------------------------------------
* lpdensity(): Saving results to new variables and plotting
*-------------------------------------------------------------------------------
lpdensity data, grid(grid) bw(0.5) genvars(lpdTemp)

twoway (rarea lpdTemp_CI_l lpdTemp_CI_r lpdTemp_grid, sort color(red%30))              ///
	   (line lpdTemp_f_p lpdTemp_grid, sort lcolor(red) lwidth("medthin") lpattern(solid)) ///
	   , ///
	   legend(off) title("lpdensity (p=2, q=3, v=1)", color(gs0)) xtitle("data") ytitle("")

drop lpdTemp_*

*-------------------------------------------------------------------------------
* lpdensity(): Adding histogram; plotting 90% uniform confidence band
*-------------------------------------------------------------------------------
lpdensity data, grid(grid) bw(0.5) plot histogram // add histogram
lpdensity data, grid(grid) bw(0.5) plot ciuniform level(90) // 90% uniform confidence band
		
*-------------------------------------------------------------------------------
* lpbwdensity(): Bandwidth selection on a provided grid
*-------------------------------------------------------------------------------
capture drop grid
gen grid = -0.5 + 0.5 * _n if _n <= 9

lpbwdensity data, grid(grid) // MSE-optimal bandwidth
lpbwdensity data, grid(grid) bwselect(imse-dpi) // IMSE-optimal bandwidth
