/*	See 
		https://bfschaffner.github.io/bookdown-stata/design-weights.html
	for details on using survey weights in Stata
*/

* Initial settings
version 17

clear all
set more off

cd "~/Dropbox/PhD/GTA/2022-23/EC317/LT/Classes/Class 5"

* Load data
use "./data/ps10.dta", clear

* Question 1
reg lhw1 tu1 age1 age1sq ib(1).sex white [pweight = weight], robust
est sto q1

* Question 2
reg dlhw du [pweight = weight], robust
est sto q2

* Question 3
reg dlhw nu un [pweight = weight], robust
est sto q3

* Question 4
test nu = -un
estadd r(p)

qui etable, est(q*) stars(.1 * .05 ** .01 ***, dec)	mstat(N, nformat(%12.0gc))
collect label dim etable_depvar "Dependent variable", modify															// label dimension etable_depvar to have a common column header
#delimit ;
	collect label levels etable_depvar																					// label levels of dimension etable_depvar to use as column headers
									1		"log(hourly wage)"
									2		"Δlog(hourly wage)"
									3		"Δlog(hourly wage)"
									, modify
	;
	collect label levels cmdset																							// label levels of dimension etable_depvar to use as column numbers
									1		"(1)"
									2		"(2)"
									3		"(3)"
									4		"(4)"
									, modify
	;
	collect label levels colname																						// label levels of dimension colname to edit how regressors are labeled in the table
									tu1		"Union member [wave 1]"
									du		"ΔUnion membership [wave 1 -> wave 5]"
									age1	"Age [wave 1]"
									nu		"Union joiner"
									un		"Union leaver"
									age1sq	"Age squared [wave 1]"
									sex		"Female"
									white	"White"
									female	"Dummy for female"
									iv		"Twin-report IV"															// this one is a new level that I will create or fill later, and will indicate the IV specifications
									, modify
	;
#delimit cr
collect label levels result N "Sample size" p "Symmetry test p-value", modify											// re-label the sample size
collect style header etable_depvar, title(label) level(label)															// show column-specific and common column headers of dimension etable_depvar; defined as the corresponding labels
collect style header cmdset, level(label)																				// show column-specific headers of dimension cmdset
collect style cell cell_type[column-header], border(top, pattern(solid))												// add horizontal lines below each row of column headers
collect style cell cmdset#cell_type[column-header], border(top, pattern(nil))											// remove the horizontal line above column headers of dimensionn cmdset (i.e., the column umbers)
collect style autolevels tu1 du nu un 2.sex age1 age1sq white _cons														// define the default levels of dimension colname to display if included in table layout (will define which regressors appear and in what order)
collect title "Wage effects of union membership"																		// add table title
collect style title, font(, bold)																						// modify table title font (make it bold)
collect style header colname[2.sex], level(hide)
qui collect layout (colname[tu1 du nu un 2.sex age1 age1sq white _cons]#result[_r_b _r_se] result[N p]) ((etable_depvar#cmdset)#stars)
collect preview																											// displays the table
