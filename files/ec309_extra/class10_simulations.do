***************************
***	0. Initial settings	***
***************************

version 17																												// declare version since I use data frames (available only since Stata 16)

clear all																												// clear memory
set more off																											// just in case it was not set off permanently in current computer

timer on 1																												// start timer to measure running time

set seed 786542319																										// set seed for reproducibility

* Set global for paths
pwd																														// pwd is unnecessary in macOS but required in Windows since `c(pwd)' is empty
if regexm("`c(pwd)'", "Pinjas") {
	global path "~/Library/CloudStorage/Dropbox/PhD/GTA/2022-23/EC309/Classes/Class 10"
}
else if "`c(pwd)'" ==  "H:") {																							// you can add an -else if- clause for your own working directory
	global path "H:/EC309/class 10"
}

* Set current directory
cd "$path"

/*
	I will use many global macros. I would recommend using local macros
	instead for a do file that is meant to be executed in one run. Here,
	global macros allow me to run the code in blocks without need to
	redefine these macros each time.
*/

************************************************************************************************************************

***************************************************************
***	1. Simulation with uncorrelated product characteristics	***
***************************************************************

/*
	Here, we will artificially draw iid uniform(0,1) prices that are
	statistically independent from product characteristics. This is an
	artificial assumption given our model of Bertrand-Nash price
	competition, but serves to illustrate a few points. We will solve
	the model and generate and endogenous vector of equilibrium prices
	in section 2.
	
	We will simulate M independent markets where a random subset of
	the J firms compete.
*/

*************************
*** 1.1 Simulate data ***
*************************

* Parameters
mata:																													// I use Mata since it is easier to work with matrices to simulate this type of data (would need to generate many variables in Stata)
	I = 1000																											// number of consumers (we will not use this parameter in this section but I prefer to define it now and store it together with other parameters)
	st_global("I",strofreal(I))																							// store in global macro for later use
	J = 50																												// number of products
	st_global("J",strofreal(J))																							// store in global macro for later use
	K = 6																												// number of non-price characteristics
	st_global("K",strofreal(K))																							// store in global macro for later use
	M = 1000																											// number of markets
	st_global("M",strofreal(M))																							// store in global macro for later use
	Nmin = 3																											// maximum number of firms per market
	Nmax = min((30,J))																									// maximum number of firms per market
	alpha = -abs(rnormal(1,1,0,1))																						// draw price utitlity coefficient (I use -abs() to make sure it is negative)
	beta = abs(rnormal(K,1,0,1))																						// draw vector of non-price utility coefficients
end

* Primitives
mata:
	X = J(1,K,0) \ runiform(J,K)																						// matrix of product (non-price) characteristics (one row vector per product; first row is outside option)
end

* Simulate data for M independent markets
mata:
	N 		= 	runiformint(M,1,Nmin,Nmax)																				// vector of number of firms per market
	cum 	= 	(0::M) :+ (0 \ runningsum(N))																			// cumulative number of observations while looping over markets (with 0 for market = -1 which is lag == 0 in the loop)
	obs 	= 	sum(N) + M																								// total number of observations (the outside option appears once per market -> M additional observations)
	st_global("obs",strofreal(obs))																						// store number of observations in global to use in Stata
	Mlist	= 	J(obs,1,.)																								// empty vector to collect the market corresponding to each observation
	F 		= 	J(obs,1,.)																								// empty vector to store the identities of firms in each market
	p 		= 	J(obs,1,.)																								// empty vector to collect prices from all markets
	s 		= 	J(obs,1,.)																								// empty vector to collect shares from all markets
end
forvalues m = 1/$M {
	local lead = `m' + 1
	mata {
		st_global("N`m'",strofreal(N[`m']))																				// store number of firms in the market in global macro for later use
		rng 			= cum[`m']+1 \ cum[`lead']																		// range of observations corresponding to market m
		Mlist[|rng|]	= J(N[`m']+1,1,`m')
		F[|rng|]		= J(1,1,0) \ jumble(1::J)[|1\N[`m']|]															// draw identities of firms participating in market m (including the outside option)
		p[|rng|] 		= J(1,1,0) \ runiform(N[`m'],1)																	// draw product price vector (first element is zero, corresponding to the outside option)
		s[|rng|]		=																								// compute vector of market shares
			exp(alpha * p[|rng|] + X[F[|rng|]:+1,.] * beta) :/ sum(exp(alpha * p[|rng|] + X[F[|rng|]:+1,.] * beta))
	}
}

**********************************
*** 1.2 Store parameter values ***
**********************************

* List of parameters
global betas																											// clear macro content (if any)
forvalues k = 1/$K {
	global betas $betas beta`k'																							// list of beta coefficients
}
local param_names I J K M Nmin Nmax alpha ${betas}																		// list of parameter names
local tempnames : list local param_names - global betas
local param_values = subinstr("`tempnames'"," ","\",.) + "\beta"

* Create data frame
frame create parameters																									// data frame to store parameter metadata

* Store parameters
frame parameters {
	
	local obs = wordcount("`param_names'")																				// number of parameters
	set obs `obs'																										// set number of observations (one per parameter)
	
	mata:
		st_addvar(("double"),("value"))																					// create Stata variable (empty)
		st_store(.,("value"),(`param_values'))																			// fill variable with parameter values
	end
	
	gen parameter = ""																									// empty string variable to store parameter names
	local i = 1																											// counter to loop over observations (one for each parameter)
	foreach name of local param_names {																					// loop over list of parameter names
		replace parameter = "`name'" in `i'																				// fill variable -names- with corresponding parameter names
		local ++i																										// increase counter
	}
	
	order parameter value																								// order variables to have parameter names first and values after
	
	*save "./simulated data/parameters.dta", replace																		// save data frame to a .dta file (just in case I am not able to simulate the data during the class)
}

*******************************
*** 1.3 Load simulated data ***
*******************************

* List of product characteristics
global chars																											// clear macro content (if any)
forvalues k = 1/$K {
	global chars = `"${chars}"' + " " + `""x`k'""'																		// list of product characteristic variable names
}
global chars = subinstr(ltrim(`"$chars"')," ",", ",.)																	// add comas to the list (necessary for st_addvar() and st_store() syntax)

* Create data frame
frame create uncorrelated																								// data frame to load simulated data

* Load data
frame uncorrelated {
	set obs ${obs}																										// set number of observations of data frame
	mata:
		st_addvar(("int","double"),("product_id","market_id","share","price",${chars}))									// create Stata variables (empty)
		st_store(
				.,																										// fill all observations
				("product_id","market_id","share","price",${chars}),													// names of variables to fill 
				(F,Mlist,s,p,X[F:+1,.]))																				// fill variables with corresponding matrices
	end
	*save "./simulated data/uncorrelated.dta", replace																	// save data frame to a .dta file (just in case I am not able to simulate the data during the class)
}

********************
*** 1.4 Analysis ***
********************

/*	If all product characteristics are observed, we can
	recover the indirect utility parameters from projection
	of deltas (can recover them from market shares)
*/

* Scatter plots of price and each characteristic
local cols = ceil(${K}/2)																								// number of graph columns
cap graph drop gu*																										// clear graphs from memory
local graphs																											// empty local to store list of subgraphs
forvalues k = 1/$K {
	local ttl  {it:{stSerif:x{sub:`k'}}}
	#delimit ;
		frame uncorrelated: twoway
			(
				scatter price x`k' if product_id != 0,
					scheme(s1mono)
					color(lavender%15)
					msize(tiny)
			)
			(
				lfit 	price x`k' if product_id != 0,
					lwidth(thick)
					lcolor(emerald)
			),
					xtitle(
						Product characteristic,
							size(small)
					)
					ytitle(
						Product price,
							size(small)
					)
					xlabel(
						,
								grid
								glcolor(gs15)
								labsize(vsmall)
					)
					ylabel(
								,
								angle(0)
								grid
								glcolor(gs15)
								labsize(vsmall)
					)
					title(
						`ttl',
							size(medium)
							margin(medsmall)
					)
					legend(off)
					name(gu`k')
					nodraw
		;
	#delimit cr
	local graphs `graphs' gu`k'
}
#delimit ;
	graph combine `graphs',
					cols(`cols')
					scheme(s1mono)
					title(Uncorrelated Product Characteristics)
					name(uncorrelated)
	;
#delimit cr
*graph export "./graphs/pricex_uncorrelated.pdf", as(pdf) replace														// just to add the graph to my slides

* Recover deltas from market shares
qui {																													// just to omit a lot of output from the -gen- and -replace- commands
	frame uncorrelated {
		gen logshare = ln(share)																						// logarithm of market share
		sum logshare if product_id == 0, meanonly																		// obtain logshare of outside option
		scalar logshare0 = r(mean)																						// store log share of outside option
		gen delta = logshare - scalar(logshare0)																		// invert market share equations to recover deltas
		drop if product_id == 0																							// drop observations corresponding to outside option
	}
}

/*
	While we have data from M different markets, notice that the
	inversion of the market share equations implies that parameters
	are identified from the shares in any one market, no matter the
	sample size (provided it is larger than the number of parameters
	including the constant).
	
	For example, we notice that there are only 8 firms competing in
	market 981 (plus the outside option, which we already excluded).
	We can recover the parameters exactly from a regression of delta
	on price and characteristics on these 8 observations.
*/

* Check number of firms in market m = 981
frame uncorrelated: tab market if market == 981

* Recover utility parameters from only one market
cap est drop uc_nu1m																									// clear estimates
frame uncorrelated: reg delta price x* if market == 981																	// OLS regression of delta on price and product characteristics
estadd local fe "no"																									// indicate no FE
est sto uc_nu1m	

* Recover utility parameters using all markets + market fixed effects
cap est drop uc_nuallfe																									// clear estimates
frame uncorrelated: reghdfe delta price x*, absorb(market) nocons														// OLS regression of delta on price and product characteristics
estadd local fe "yes"																									// indicate FE
est sto uc_nuallfe	

* Estimate utility parameters using all markets without market fixed effects
cap est drop uc_nuallnofe																									// clear estimates
frame uncorrelated: reg delta price x*																					// OLS regression of delta on price and product characteristics
estadd local fe "no"																									// indicate no FE
est sto uc_nuallnofe																									// store results

/*	If some characteristics are unobserved but uncorrelated to prices
	(and other observed characteristics), we can consistently estimates
	utility parameters with OLS regression of deltas on price and observed
	characteristics. We will assume that only x1-xL are observed, where
	L = K - 2 -> 2 product characteristics are unobserved.
*/

* Number of observed characteristics
local L = $K - 2																										// 2 product characteristics are unobserved

* Estimate utility parameters including market fixed effects
cap est drop uc_u																										// clear estimates
frame uncorrelated: reghdfe delta price x1-x`L', absorb(market)	vce(robust)												// OLS regression of delta on observed characteristics (incl. price) with market fixed effects
estadd local fe "yes"																									// indicate FE
est sto uc_ufe																											// store results

* Estimate utility parameters without market fixed effects
cap est drop uc_u																										// clear estimates
frame uncorrelated: reg delta price x1-x`L', robust																		// OLS regression of delta on observed characteristics (incl. price)
estadd local fe "no"																									// indicate no FE
est sto uc_unofe	


* Table settings
local params alpha ${betas}																								// list of parameters to tabulate
local k = 0																												// counter for product characteristics
foreach param of local params {																							// this loop is just to include the true parameters (stored in frame parameters) in the table
	frame parameters: sum value if parameter == "`param'", meanonly
	scalar b`k' = r(mean)
	local true `true' `: display %10.7g b`k''
	if "`param'" == "alpha" {
		local varlab price "alpha"
	}
	else {
		local varlab `varlab' x`k' "beta`k'"
	}
	local ++k

}

*Tabulate true estimates and regression results
#delimit ;
	estout uc_nu1m uc_nuallfe uc_nuallnofe uc_unofe uc_ufe,
					labcol2(
						`true',
							title(True value)
					)
					cells(
						b(
							nostar
							fmt(%10.7g)
						)
						se(
							par
							fmt(%10.7g)
						)
					)
					varlabels(
						`varlab',
					)
					lz
					mlabels(
							`"Ideal case: 1 m"'
							`"Ideal case: all"'
							`"Ideal case: all"'
							`"Unobserved chars."'
							`"Unobserved chars."'
					)
					stats(
							N
							fe
							r2,
								fmt(
									%12.0fc
									%s
									%4.3g
								)
								labels(
									"Observations"
									"Market FE"
									"R squared"
								)
					)
					drop(_cons)
					modelwidth(18)
					numbers
					collabels(none)
	;
#delimit cr

************************************************************************************************************************

**************************************************************************+
***	2. Simulation with correlated unobserved product characteristics	***
***************************************************************************

*************************
*** 2.1 Simulate data ***
*************************

/*	
	Now, we simulate the full model (including the supply side)
	to obtain data generated by equilibrium in many markets. From
	the model, we know that equilibrium prices will be a function
	of product characteristics. 
	
	We will work with the previously defined values for J, K, and
	I. Similarly, we will work with the same parameter values for
	alpha and beta. We will also work with the previously drawn
	matrix of product characteristics X, but we will assume that
	only a random subset of firms participate in each market.
	The additional primitives we will draw are:
	
					- vector of marginal costs: c
					- list of participating firms for market m: Fm
					
	We will compute the equilibrium prices by numerically solving
	the Nash equilibrium of the oligopoly game.
	
	Aknowledgement:
	
		This part of the code is based on Matteo Courthoud's logit
		demand simulation, available at:
		
			https://matteocourthoud.github.io/course/empirical-io/11_logit_demand/
*/

* Functions to solve model in each market
mata:
	mata set matastrict off																								// just to make sure that declaration of all variables is not requ
	real colvector function demand(																						// we define the demand vector function, declaring it to be a real column vector
					real colvector	p,																					// argument p is a real column vector
					real matrix		X,																					// argument x is a real matrix
					real scalar		alpha,																				// argument alpha is a real scalar
					real colvector	beta,																				// argument beta is a real column vector
					real scalar		I																					// argument I is a real scalar
			) {
					return(
						I * exp(alpha * p + X * beta):/sum(exp(alpha * p + X * beta))									// the function takes the arguments p, X, alpha, beta, I and returns the vector of demands
					)
			}
	real colvector function profit(																						// and now you understand the syntax for defining the profit vector function
					real colvector	p,
					real colvector	c,																					// marginal cost vector
					real matrix		X,
					real scalar		alpha,
					real colvector	beta,
					real scalar		I
			) {
				return(
					(p - c) :* demand(p, X, alpha, beta, I)																// returns the profit vector where quantities are determined by the demand function
				)	
			}
	real scalar function profit_dev(																					// this function computes the profits of of firm j after a unilateral change in its price
					real scalar		j,																					// identity of deviating firm
					real scalar		p_j,																				// new price set by firm j (the deviation)
					real colvector	p,
					real colvector	c,
					real matrix		X,
					real scalar		alpha,
					real colvector	beta,
					real scalar		I
			) {
					p[j + 1]	=	p_j																					// here I change the (j + 1)th component since the first component corresponds to the outside option (with p_0 = p[1] = 0)
					profits		=	profit(p, c, X, alpha, beta, I)														// compute vector of implied profits for all firms
					return(
						profits[j + 1]																					// returned implied profit of firm j
					)
			}
	void best_response(		real scalar todo,																			// define evaluator for best response function (which we will optimize for each firm)	
							real scalar		p_j,
							real scalar		j,
							real colvector	p,
							real colvector	c,
							real matrix		X,
							real scalar		alpha,
							real colvector 	beta,
							real scalar		I,
							real scalar		prof_j,
							real matrix		g,
							real matrix		H
					) {
						prof_j = profit_dev(j, p_j, p, c, X, alpha, beta, I)
					}
	real colvector function eq_price(																					// this function will compute the Bertrand-Nash equilibrium prices given the primitives c and X
					real colvector	c,
					real matrix		X,
					real scalar		alpha,
					real colvector	beta,
					real scalar		I,
					real scalar		J
			) {
				p = 1 * c																								// initial price vector we could specify the initial price, tolerance, and max iterations as arguments of the function but will not do so for simplicity
				tol = 1e-2																								// tolerance for convergence (we will stop when the vector of prices canges by les than 'tol'; in a real application you would set a far more stringent tolerance in the order of 1e-6, 1e-8 or even 1e-12, but here I am concerned about speed)
				maxiter = 1e3																							// maximum # of iterations (we will stop if we reach the 'maxiter'th iteation without convergence)
				iter = 0																								// initialize iteration counter
				dist = 1																								// initialize distance for convergence criterion
				new_p = p																								// initialize vector to collect best responses (new prices)
				while (dist > tol & iter < maxiter) {																	// iterate while we do not have convergence AND we have not reached the maximum iterations (i.e., stop if we achieve convergence OR exceed max iterations)
					S = optimize_init()																					// define/start optimization problem (this part is common to all firms)
					optimize_init_evaluator(S, &best_response())														// declare evaluator for objective function (this part is common to all firms)
					optimize_init_evaluatortype(S, "d0")																// declare that optimization is fully numerical (we do not pass gradients nor hessians to the optimization routine)
					optimize_init_argument(S, 2, p)																		// declare arguments that are common to all firms
					optimize_init_argument(S, 3, c)
					optimize_init_argument(S, 4, X)
					optimize_init_argument(S, 5, alpha)
					optimize_init_argument(S, 6, beta)
					optimize_init_argument(S, 7, I)
					for (j = 1; j <= J; j++) {																			// iterate over firms
						optimize_init_argument(S, 1, j)																	// declare firm identity argument						
						optimize_init_params(S, p[j+1])																	// choose starting value for p_j
						new_p[j + 1] = optimize(S)																		// update price of firm j with its best response given current prices
					}
					dist = max(abs(new_p - p))																			// update distance
					p = new_p																							// update prices
					++iter																								// update iteration counter
				}
				return(p)																								// return equilibrium price vector
			}
end

* Primitives
mata: 
	p = J(obs,1,.)																										// prices will now be equilibrium outcomes
	s = J(obs,1,.)																										// shares will now be equilibrium outcomes
	c = 0 \ runiform(J,1)																								// draw marginal cost vector (first element is zero, corresponding to the outside option)
end

* Simulate data for M independent markets
forvalues m = 1/$M {
	local lead = `m' + 1
	mata {
		rng = cum[`m']+1 \ cum[`lead']																					// range of observations corresponding to market m
		p[|rng|] = eq_price(c[F[|rng|]:+1,.], X[F[|rng|]:+1,.], alpha, beta, I, N[`m'])									// compute equilibrium price vector
		s[|rng|] = demand(p[|rng|], X[F[|rng|]:+1,.], alpha, beta, I) / I												// compute equilibrium market share vector
	}
	di "*************************************************"
	di " # of firms competing in market `m': " ${N`m'}
	di "*************************************************"
}

*******************************
*** 2.2 Load simulated data ***
*******************************

* Frame to collect data
frame create full_model																									// data frame to load simulated data

* Load data
frame full_model {
	set obs ${obs}																										// set number of observations of data frame
	mata:
		st_addvar(("int","double"),("product_id","market_id","share","price",${chars}))									// create Stata variables (empty)
		st_store(
				.,																										// fill all observations
				("product_id","market_id","share","price",${chars}),													// names of variables to fill 
				(F,Mlist,s,p,X[F:+1,.]))																				// fill variables with corresponding matrices
	end
	*save "./simulated data/full_model.dta", replace																		// save data frame to a .dta file (just in case I am not able to simulate the data during the class)
}

********************
*** 2.3 Analysis ***
********************

* Scatter plots of price and each characteristic
local cols = ceil(${K}/2)																								// number of graph columns
cap graph drop gfm*	full_model																							// clear graphs from memory
local graphs																											// empty local to store list of subgraphs
forvalues k = 1/$K {
	if inlist(`k', 2) {
		local color`k' ebblue%15																						// this is an ad-hoc way of highlighting
	}																													// the strong correlation between price
	else {																												// and (x1,x4), which I am aware of only
		local color`k' lavender%15																						// after having looked at these scatterplots.
	}
	local ttl  {it:{stSerif:x{sub:`k'}}}
	#delimit ;
		frame full_model: twoway
			(
				scatter price x`k' if product_id != 0,
					scheme(s1mono)
					color(`color`k'')
					msize(tiny)
			)
			(
				lfit price x`k' if product_id != 0,
					lwidth(thick)
					lcolor(emerald)
			),
					xtitle(
						Product characteristic,
							size(small)
					)
					ytitle(
						Product price,
							size(small)
					)
					xlabel(
						,
								grid
								glcolor(gs15)
								labsize(vsmall)
					)
					ylabel(
								,
								angle(0)
								grid
								glcolor(gs15)
								labsize(vsmall)
					)
					title(
						`ttl',
							size(medium)
							margin(medsmall)
					)
					legend(off)
					name(gfm`k')
					nodraw
		;
	#delimit cr
	local graphs `graphs' gfm`k'
}
#delimit ;
	graph combine `graphs',
					cols(`cols')
					scheme(s1mono)
					title(Full model, margin(zero))
					name(full_model)
	;
#delimit cr
*graph export "./graphs/pricex_full_model.pdf", as(pdf) replace															// just to add the graph to my slides

/*	Looking at the scatter plots, the correlation between price
	and x2 look reasonably strong. Therefore, we will assume
	that x2 is unobserved, rendering price endogenous.
*/

* Lists of observed and unobserved characteristics
global xunobs x2																										// we specify the unobserved characteristics
frame full_model: qui ds x*																								// then, we get the full list of characteristics
local xall = "`r(varlist)'"																								// and store it in a local macro
global xobs : list local xall - global xunobs																			// then, we obtain the list of observed characteristics

* Recover deltas from market shares
qui {																													// just to omit a lot of output from the -gen- and -replace- commands
	frame full_model {
		gen delta = .																									// empty variable to store deltas from each market
		gen logshare = ln(share)																						// logarithm of market share
		forvalues m = 1/$M {																							// loop over markets
			sum logshare if product_id == 0 & market == `m', meanonly													// obtain logshare of outside option
			scalar logshare0_`m' = r(mean)																				// store log share of outside option
			replace delta = logshare - scalar(logshare0_`m') if market == `m'											// invert market share equations to recover deltas
		}
		drop if product_id == 0
	}
}

* Ideal regression
cap est drop fm_ideal																									// clear estimates
frame full_model: reg delta price x*, nocons																			// OLS regression of delta on all characteristics (incl. price)
est sto fm_ideal																										// store results

/*
	Notice that we can recover the utility parameters exactly from a
	regression of delta on price and all characteristics using data
	from just one market (provided sufficient observations). For
	example, there are only 8 observations from market 981, yet we
	are able to recover the full vector of utility parameters.
*/


* Check number of firms in market m = 981
frame full_model: tab market if market == 981

* Ideal regression from one small market
cap est drop fm_ideal1m																									// clear estimates
frame full_model: reg delta price x* if market == 981, nocons															// OLS regression of delta on all characteristics (incl. price) using only observations from market 948
est sto fm_ideal1m

* OLS regression
cap est drop fm_ols																										// clear estimates
frame full_model: reg delta price ${xobs}, robust																		// OLS regression of delta on observed characteristics (incl. price)
est sto fm_ols																											// store results

* BLP instruments
local indices = subinstr(subinstr("$xobs","x","",.),"-","/",.)															// remove the x's from the list of observed varibles and substitute / for - in order to use it as a numlist
qui {																													// just to omit a lot of output from the -gen- and -replace- commands
	frame full_model {
		foreach k of numlist `indices' {
			egen z`k' = total(x`k'), by(market)
			replace z`k' = z`k' - x`k'
		}
		gen z0 = .
		forvalues m = 1/$M {
			replace z0 = ${N`m'} - 1 if market == `m'
		}
	}
}

* Scatter plots of price and unobserved char. vs BLP instruments
frame full_model: qui ds z*																								// this will generate a list of all the instruments
local nz = wordcount("`r(varlist)'")																					// we store the number ow words in this list (i.e., number of instruments) in a local
local cols = ceil(`nz'/2)																								// number of graph columns
local yttl_price Price
frame full_model {
	foreach var of varlist $xunobs {
		local yttl_`var' = ///
			"Unobserved characteristic {it:{stSerif:" + substr("`var'",1,1) + "{sub:" + substr("`var'",2,.) + "}}}"	
	}
	foreach var of varlist price $xunobs {																				// loop over price and unobserved variables
		cap graph drop `var'z*																							// clear graphs from memory
		local graphs																									// empty local to store list of subgraphs
		foreach inst of varlist z* {																					// loop over BLP instruments
			local varname = substr("`inst'",1,1) + "{sub:" + substr("`inst'",2,.) + "}"									// just to display zk as z_k (with subindices)
			local ttl {it:{stSerif:`varname'}}
			#delimit ;
				frame full_model: twoway
					(
						scatter `var' `inst' if product_id != 0,
							scheme(s1mono)
							color(lavender%15)
							msize(tiny)
					)
					(
						lfit 	`var' `inst' if product_id != 0,
							lwidth(thick)
							lcolor(emerald)
					),
							xtitle(
								Instrument,
									size(small)
							)
							ytitle(
								`yttl_`var'',
									size(small)
							)
							xlabel(
								,
										grid
										glcolor(gs15)
										labsize(vsmall)
							)
							ylabel(
										,
										angle(0)
										grid
										glcolor(gs15)
										labsize(vsmall)
							)
							title(
								`ttl',
									size(medium)
									margin(medsmall)
							)
							legend(off)
							name(`var'`inst')
							nodraw
				;
			#delimit cr
			local graphs `graphs' `var'`inst'
		}
		#delimit ;
			graph combine `graphs',
							cols(`cols')
							scheme(s1mono)
							title(`yttl_`var'' vs BLP Instruments)
							name(`var'z)
			;
		#delimit cr
		*graph export "./graphs/`var'z_full_model.pdf", as(pdf) replace													// just to add the graph to my slides
	}
}

* 2SLS regression
cap est drop fm_2sls																									// clear estimates
frame full_model: ivregress 2sls delta ${xobs} (price = z*), robust														// 2SLS regression of delta on observed characteristics (incl. price)
est sto fm_2sls																											// store results

* GMM regression
cap est drop fm_gmm																										// clear estimates
frame full_model: ivregress gmm delta ${xobs} (price = z*)																// GMM regression of delta on observed characteristics (incl. price)
est sto fm_gmm																											// store results

* Table settings
local params alpha ${betas}																								// list of parameters to tabulate
local k = 0																												// counter for product characteristics
foreach param of local params {																							// this loop is just to include the true parameters (stored in frame parameters) in the table
	frame parameters: sum value if parameter == "`param'", meanonly
	scalar b`k' = r(mean)
	local true `true' `: display %09.7f b`k''
	if "`param'" == "alpha" {
		local varlab price "alpha"
	}
	else {
		local varlab `varlab' x`k' "beta`k'"
	}
	local ++k
}

*Tabulate true estimates and regression results
#delimit ;
	estout fm_ideal fm_ideal1m fm_ols fm_2sls fm_gmm,
					labcol2(
						`true',
							title(True value)
					)
					cells(
						b(
							nostar
							fmt(%010.7g)
						)
						se(
							par
							fmt(%010.7g)
						)
					)
					varlabels(
						`varlab',
					)
					lz
					mlabels(
							`"Ideal Case"'
							`"Ideal w/1 market"'
							`"OLS"'
							`"2SLS"'
							`"GMM"'
					)
					stats(
							N
							r2,
								fmt(
									%12.0fc
									%4.3g
								)
								labels(
									"Observations"
									"R squared"
								)
					)
					keep(price x*)
					modelwidth(18)
					numbers
					collabels(none)
	;
#delimit cr

************************************************************************************************************************

**************************************************************************+
***	3. Simulation of microdata with correlated product characteristics	***
***************************************************************************

*************************
*** 3.1 Simulate data ***
*************************

* Simulate consumer choices
mata:
	mata drop c s																										// we will not need the marginal cost and market share data
	NM = 100																											// number of markets to sample
	Msample = sort(jumble(1::M)[|1\NM|],1)																				// random sample of markets
	obs = NM * I																										// we will sample the same number I of consumers in each market
	st_global("obs",strofreal(obs))																						// store number of observations in global to use in Stata
	choice = J(obs,1,.)																									// empty vector to store consumer choices (1 product per consumer)
	markets = J(${obs},1,.)
	for (n = 1; n <= NM; n++) {																							// loop over sampled markets
		m = Msample[n]
		markets[|(n-1)*I+1\(n-1)*I+I|] = J(I,1,m)
		rng = cum[m]+1 \ cum[m+1]																						// range of observations in vector p corresponding to market m
		firms = F[|rng|]																								// identities of products in the market
		utility = J(I,1,1) # (alpha*p[|rng|] + X[F[|rng|]:+1,.]*beta)' -ln(-ln(runiform(I,N[m]+1)))						// we are computing mean utilities and drawing idiosyncratic shocks in one line
		for (i = 1; i <= I; i++) {																						// loop over consumers
			index = selectindex(utility[i,.]':==rowmax(utility[i,.]))													// search for the row-index of the maximum utility in the column (transposed) vector of utilities of consumer i and store it in row i of the choice vector (I subtract 1 because of the shift in indices due to inclusion of the outside option)
			choice[(n-1)*I+i,1] = 	firms[index]																		// identity of chosen product
		}
	}
end

*******************************
*** 3.2 Load simulated data ***
*******************************

* Create data frames
frame create micro_choices																								// data frame to load consumer choices
frame create micro_choicesets																							// data frame to load choice sets (products available in each market)

* Load consumer data
frame micro_choices {
	set obs ${obs}																										// set number of observations of data frame
	mata:
		st_addvar(("long","int","double"),("consumer_id","market_id","choice"))											// create Stata variables (empty)
		st_store(.,("consumer_id","market_id","choice"),((1::obs),markets,choice))										// fill variables with corresponding matrices
	end
	*save "./simulated data/micro_cohices.dta", replace																	// save data frame to a .dta file (just in case I am not able to simulate the data during the class)
}

* Load choice set data
mata:
	st_global("obs",strofreal(sum(N[Msample]:+1)))
	Nsample = 0 \ N[Msample]
	firms = J(${obs},1,.)
	markets = J(${obs},1,.)
	prices = J(${obs},1,.)
	run = runningsum(Nsample:+1)
	for (n = 1; n <= NM; n++) {
		m = Msample[n]
		rng = run[n] \ run[n] + Nsample[n+1]
		markets[|rng|] = J(Nsample[n+1]+1,1,m)
		firms[|rng|] = F[|cum[m]+1 \ cum[m+1]|]
		prices[|rng|] = p[|cum[m]+1 \ cum[m+1]|]
	}
	
end

frame micro_choicesets {
	set obs ${obs}
	mata:
		st_addvar(("int","int","double"),("market_id","product_id","price",${chars}))									// create Stata variables (empty)
		st_store(.,("market_id","product_id","price",${chars}),(markets,firms,prices,X[firms:+1,.]))					// fill variables with corresponding matrices
	end
	*save "./simulated data/micro_cohicesets.dta", replace																// save data frame to a .dta file (just in case I am not able to simulate the data during the class)
}

********************
*** 3.3 Analysis ***
********************

* Format choice data for conditional logit estimation
frame micro_choicesets: frame put market_id, into(expand)
frame expand {
	bysort market_id: gen cs_size = _N
	by market_id: keep if _n == 1
}
frame micro_choices {
	frlink m:1 market_id, frame(expand)
	frget cs_size, from(expand)
	drop expand
	frame drop expand
	expand cs_size
	drop cs_size
	levelsof market_id, local(mkts)
	gen product_id = .
	sort market_id consumer_id
	mata {
		st_view(X,.,"product_id")
		lst = 0
		for (n=1; n <= NM; n++) {
			fst = lst + 1
			lst = lst + (Nsample[n+1]:+1)*I
			rngX = fst \ lst
			rngf = run[n] \ run[n] + Nsample[n+1]
			X[|rngX|] = J(I,1,1) # firms[|rngf|]
		}
	}
	mata: mata clear
	replace choice = product_id == choice
	cmset consumer_id product_id
	frlink m:1 market_id product_id, frame(micro_choicesets)
	frget price x*, from(micro_choicesets)
	drop micro_choicesets
	sort market_id consumer_id product_id
}

* MLE: Ideal case all markets
cap est drop micro_fullall
frame micro_choices: cmclogit choice price x*, nocons
est sto micro_fullall

* MLE: Ideal case 1 market
cap est drop micro_full1m
frame micro_choices: cmclogit choice price x* if market == 196, nocons
est sto micro_full1m

* Lists of observed and unobserved characteristics
global xunobs x1 x4 x6																										// we specify the unobserved characteristics
frame full_model: qui ds x*																								// then, we get the full list of characteristics
local xall = "`r(varlist)'"																								// and store it in a local macro
global xobs : list local xall - global xunobs																			// then, we obtain the list of observed characteristics

* MLE: unobserved characteristics
cap est drop micro_unobs
frame micro_choices: cmclogit choice price ${xobs}, nocons
est sto micro_unobs

* Two-step procedure (estimate deltas as constant from MLE and then do OLS)
/*frame create micro_deltas product_id market_id delta
frame micro_choicesets: levelsof market_id, local(mkts)
foreach mkt of local mkts {
	frame micro_choices{
		cmclogit choice if market == `mkt'
		levelsof product_id if market == `mkt', local(pdcts)
	}
	foreach pdct of local pdcts {
		cap frame post micro_deltas (`pdct') (`mkt') (e(b)["y1","`pdct':_cons"])
		if scalar(_rc) != 0 {
			frame post micro_deltas (`pdct') (`mkt') (0)
		}
	}
}*/

/*
	I am just loading the data generated by the code commented above
	since it takes more than 15 minutes to run. There are ways to
	make it run faster with parallel computing but I will not delve
	into this due to time constraints.
*/

frame create micro_deltas
frame micro_deltas: use "./simulated data/micro_deltas.dta"

frame micro_deltas {
	/*frlink 1:1 product_id market_id, frame(micro_choicesets)
	frget price x*, from(micro_choicesets)
	drop micro_choicesets
	*save "./simulated data/micro_deltas.dta", replace																	// save data frame to a .dta file (just in case I am not able to simulate the data during the class)
	*/
	* Ideal case
	cap est drop twostep_ideal
	reg delta price x*, robust
	est sto twostep_ideal
	
	* Unobserved characteristics
	cap est drop twostep_unobsols
	reg delta price ${xobs}, robust
	est sto twostep_unobsols
}

* BLP instruments
local indices = subinstr("$xobs","x","",.) 																				// remove the x's from the list of observed varibles and substitute / for - in order to use it as a numlist
qui {																													// just to omit a lot of output from the -gen- and -replace- commands
	frame micro_deltas {
		foreach k of numlist `indices' {
			egen z`k' = total(x`k'), by(market_id)
			replace z`k' = z`k' - x`k'
		}
		egen z0 = total(1), by(market_id)
		replace z0 = z0 - 1
	}
}

* 2SLS with BLP instruments in second step
cap est drop twostep_unobs2sls
frame micro_deltas: ivregress 2sls delta ${xobs} (price = z*), vce(robust)
est sto twostep_unobs2sls

* Table settings
local params alpha ${betas}																								// list of parameters to tabulate
local k = 0																												// counter for product characteristics
foreach param of local params {																							// this loop is just to include the true parameters (stored in frame parameters) in the table
	frame parameters: sum value if parameter == "`param'", meanonly
	scalar b`k' = r(mean)
	local true `true' `: display %09.7f b`k''
	if "`param'" == "alpha" {
		local varlab price "alpha"
	}
	else {
		local varlab `varlab' x`k' "beta`k'"
	}
	local ++k
}

*Tabulate true estimates and regression results
#delimit ;
	estout micro_fullall micro_full1m micro_unobs twostep_ideal twostep_unobsols twostep_unobs2sls,
					labcol2(
						`true',
							title(True value)
					)
					cells(
						b(
							nostar
							fmt(%010.7g)
						)
						se(
							par
							fmt(%010.7g)
						)
					)
					varlabels(
						`varlab',
					)
					lz
					mlabels(
							`"MLE full"'
							`"MLE full 1 market"'
							`"MLE unobs"'
							`"2S full"'
							`"2S unobs OLS"'
							`"2S unobs 2SLS"'
							
					)
					stats(
							N
							r2,
								fmt(
									%12.0fc
									%4.3g
								)
								labels(
									"Observations"
									"R squared"
								)
					)
					keep(price x*)
					modelwidth(18)
					numbers
					collabels(none)
	;
#delimit cr

* Scatter plots of price and each characteristic
local cols = ceil(${K}/2)																								// number of graph columns
cap graph drop micro*																									// clear graphs from memory
local graphs																											// empty local to store list of subgraphs
forvalues k = 1/$K {
	if inlist(`k', 1, 4, 6) {
		local color`k' ebblue%15																						// this is an ad-hoc way of highlighting
	}																													// the strong correlation between price
	else {																												// and (x1,x4), which I am aware of only
		local color`k' lavender%15																						// after having looked at these scatterplots.
	}
	local ttl  {it:{stSerif:x{sub:`k'}}}
	#delimit ;
		frame micro_deltas: twoway
			(
				scatter price x`k' if product_id != 0,
					scheme(s1mono)
					color(`color`k'')
					msize(tiny)
			)
			(
				lfit price x`k' if product_id != 0,
					lwidth(thick)
					lcolor(emerald)
			),
					xtitle(
						Product characteristic,
							size(small)
					)
					ytitle(
						Product price,
							size(small)
					)
					xlabel(
						,
								grid
								glcolor(gs15)
								labsize(vsmall)
					)
					ylabel(
								,
								angle(0)
								grid
								glcolor(gs15)
								labsize(vsmall)
					)
					title(
						`ttl',
							size(medium)
							margin(medsmall)
					)
					legend(off)
					name(micro`k')
					nodraw
		;
	#delimit cr
	local graphs `graphs' micro`k'
}
#delimit ;
	graph combine `graphs',
					cols(`cols')
					scheme(s1mono)
					title(Micro data)
					name(micro)
	;
#delimit cr
*graph export "./graphs/pricex_micro.pdf", as(pdf) replace																// just to add the graph to my slides

* Scatter plots of price and unobserved char. vs BLP instruments
frame micro_deltas: qui ds z*																								// this will generate a list of all the instruments
local nz = wordcount("`r(varlist)'")																					// we store the number ow words in this list (i.e., number of instruments) in a local
local cols = ceil(`nz'/2)																								// number of graph columns
local yttl_price Price
frame micro_deltas {
	foreach var of varlist $xunobs {
		local yttl_`var' = ///
			"Unobserved characteristic {it:{stSerif:" + substr("`var'",1,1) + "{sub:" + substr("`var'",2,.) + "}}}"	
	}
	foreach var of varlist price $xunobs {																				// loop over price and unobserved variables
		cap graph drop micro`var'z*																							// clear graphs from memory
		local graphs																									// empty local to store list of subgraphs
		foreach inst of varlist z* {																					// loop over BLP instruments
			local varname = substr("`inst'",1,1) + "{sub:" + substr("`inst'",2,.) + "}"									// just to display zk as z_k (with subindices)
			local ttl {it:{stSerif:`varname'}}
			#delimit ;
				twoway
					(
						scatter `var' `inst' if product_id != 0,
							scheme(s1mono)
							color(lavender%15)
							msize(tiny)
					)
					(
						lfit 	`var' `inst' if product_id != 0,
							lwidth(thick)
							lcolor(emerald)
					),
							xtitle(
								Instrument,
									size(small)
							)
							ytitle(
								`yttl_`var'',
									size(small)
							)
							xlabel(
								,
										grid
										glcolor(gs15)
										labsize(vsmall)
							)
							ylabel(
										,
										angle(0)
										grid
										glcolor(gs15)
										labsize(vsmall)
							)
							title(
								`ttl',
									size(medium)
									margin(medsmall)
							)
							legend(off)
							name(micro`var'`inst')
							nodraw
				;
			#delimit cr
			local graphs `graphs' micro`var'`inst'
		}
		#delimit ;
			graph combine `graphs',
							cols(`cols')
							scheme(s1mono)
							title(`yttl_`var'' vs BLP Instruments)
							name(micro`var'z)
			;
		#delimit cr
		*graph export "./graphs/`var'z_micro.pdf", as(pdf) replace														// just to add the graph to my slides
	}
}

* Construct composite unobs. char. (error term)
local indices = subinstr("$xunobs","x","",.)
foreach index of local indices {
	frame parameters: sum value if parameter == "beta`index'", meanonly
	scalar b`index' = r(mean)
	local expr   "`expr' + " "scalar(b`index') * x`index'"
}
local expr = subinstr(subinstr(`"`expr'"',"+ ","",1),`"""',"",.)
frame micro_deltas: gen xi =  `expr'

* Scatter plots of price and unobserved char. vs BLP instruments
frame micro_deltas: qui ds z*																								// this will generate a list of all the instruments
local nz = wordcount("`r(varlist)'")																					// we store the number ow words in this list (i.e., number of instruments) in a local
local cols = ceil(`nz'/2)																								// number of graph columns
frame micro_deltas {
	cap graph drop microxiz*																							// clear graphs from memory
	local graphs																										// empty local to store list of subgraphs
	foreach inst of varlist z* {																						// loop over BLP instruments
		local varname = substr("`inst'",1,1) + "{sub:" + substr("`inst'",2,.) + "}"										// just to display zk as z_k (with subindices)
		local ttl {it:{stSerif:`varname'}}
		#delimit ;
			twoway
				(
					scatter xi `inst' if product_id != 0,
						scheme(s1mono)
						color(lavender%15)
						msize(tiny)
				)
				(
					lfit 	xi `inst' if product_id != 0,
						lwidth(thick)
						lcolor(emerald)
				),
						xtitle(
							Instrument,
								size(small)
						)
						ytitle(
							Composite unobs. char. {it:{&xi}{sub:{it:i}}},
								size(small)
						)
						xlabel(
							,
									grid
									glcolor(gs15)
									labsize(vsmall)
						)
						ylabel(
									,
									angle(0)
									grid
									glcolor(gs15)
									labsize(vsmall)
						)
						title(
							`ttl',
								size(medium)
								margin(medsmall)
						)
						legend(off)
						name(microxi`inst')
						nodraw
			;
		#delimit cr
		local graphs `graphs' microxi`inst'
	}
	#delimit ;
		graph combine `graphs',
						cols(`cols')
						scheme(s1mono)
						title(Regression error {it:{&xi}{sub:{it:i}}} vs BLP Instruments)
						name(microxiz)
		;
	#delimit cr
	*graph export "./graphs/xiz_micro.pdf", as(pdf) replace															// just to add the graph to my slides
}

* Stop timer
timer off 1
timer list
di r(t1)/60
