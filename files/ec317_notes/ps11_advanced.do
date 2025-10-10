******************************
*****  Initial Settings  *****
******************************

version 17
set more off
clear all
cd "~/Dropbox/PhD/GTA/2022-23/EC317/LT/Classes/Class 6"

* Load data
use "./data/ps11.dta", clear

*****************************
*****  Q3: First stage  *****
*****************************

* Declare panel
xtset dept year

* Define estimation samples
local sample1 inrange(year,1850,1905)
local sample2 `sample1' & high

* Run regressions
forvalues j = 1/2 {
	qui xtreg lwine phyll i.year if `sample`j'', fe vce(robust)
	est sto q3_`j'
}

* Produce table
qui etable, estimates(q3*) keep(phyll) stars(.1 * .05 ** .01 ***, dec) mstat(N, nformat(%12.0gc))
collect label levels result N "Sample size", modify
collect label levels colname phyll "Phylloxera {c -} full contamination", modify
collect label levels etable_depvar 1 "Whole sample" 2 "Wine-intensive {it:départments}", modify
collect style column, width(equal)
collect label dim etable_depvar "Log (wine production)", modify
collect style header etable_depvar, title(label)
collect style cell cell_type[column-header], border(bottom, pattern(solid))
collect preview

* Drop unnecessary variables
drop _est*

******************************
*****  Q4: Reduced form  *****
******************************

* Define list of outcome variables
#delimit ;
	local crimes
				violent
				property
				minor
				homicide
				theft_all
				theft_church
				theft_road
				theft_dom
				theft_other
	;
#delimit cr

* Run regressions
local sample
local j = 1
cap est drop q4*
foreach crime of local crimes {
	if ("`crime'" == "theft_all") local sample if year <= 1912
	qui xtreg `crime' phyll i.year `sample', fe vce(robust)
	est sto q4_`j'
	local ++j
}

* Produce table
local eqbase violent
local eqrecode : list local crimes - local eqbase
foreach eq of local eqrecode {
	local recode `recode' eqrecode(`eq' = `eqbase')
}
qui etable, est(q4*) keep(phyll) stars(.1 * .05 ** .01 ***, dec) mstat(N, nformat(%12.0gc)) `recode'
collect label levels result N "Sample size", modify
collect label levels colname phyll "Phylloxera", modify
collect remap etable_depvar[1 2 3] = agg
collect remap etable_depvar[4 5 6 7 8 9] = dis
collect remap etable_dvlabel[1 2 3 4] = notheft
collect remap etable_dvlabel[5 6 7 8 9] = theft
collect label dim agg "1826{c -}1936", modify
collect label dim dis "1826{c -}1912", modify
collect label dim theft "Thefts"
collect style header theft, title(label) level(label)
collect style header notheft, title(hide) level(label)
collect style header agg dis, title(label) level(hide)
collect label levels notheft 1 "Violent" 2 "Property" 3 "Minor" 4 "Homicides", modify
collect label levels theft 	 5 "All" 	 6 "Church"   7 "Road" 	8 "Domestic" 9 "Other", modify
collect recode dis 4 = 1 5 = 1 6 = 1 7 = 1 8 = 1 9 = 1
collect style cell cell_type[column-header], border(bottom, pattern(solid))
qui collect layout (coleq#colname[phyll]#result[_r_b _r_se] result[N]) ((agg dis)#(notheft theft)#stars)
collect preview

* Drop unnecessary variables
drop _est*

********************
*****  Q5: IV  *****
********************

* Run regressions
local sample inrange(year,1850,1905)
local j = 1
cap est drop q5*
foreach crime of local crimes {
	qui xtivreg `crime' i.year (lwine = phyll) if `sample', fe vce(robust)
	est sto q5_`j'
	local ++j
}

* Produce table
local eqbase violent
local eqrecode : list local crimes - local eqbase
foreach eq of local eqrecode {
	local recode `recode' eqrecode(`eq' = `eqbase')
}
qui etable, est(q5*) keep(lwine) stars(.1 * .05 ** .01 ***, dec) mstat(N, nformat(%12.0gc)) `recode'
collect label levels result N "Sample size", modify
collect label levels colname lwine "Log (wine prod.)", modify
collect remap etable_dvlabel[1 2 3 4] = notheft
collect remap etable_dvlabel[5 6 7 8 9] = theft
collect label dim theft "Thefts"
collect style header theft, title(label) level(label)
collect style header notheft, title(hide) level(label)
collect label levels notheft 1 "Violent" 2 "Property" 3 "Minor" 4 "Homicides", modify
collect label levels theft 	 5 "All" 	 6 "Church"   7 "Road" 	8 "Domestic" 9 "Other", modify
collect style cell cell_type[column-header], border(bottom, pattern(solid))
qui collect layout (coleq#colname[lwine]#result[_r_b _r_se] result[N]) ((notheft theft)#stars)
collect preview

* Drop unnecessary variables
drop _est*
