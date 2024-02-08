/*-----------------------------*/
/*ADO file*/
/* logitjack  - logit jack(knife) estimation */
/* implements procedures found in this paper: */
/* logit-jack */ 
/* written by Matt Webb */
/* version 1 */
/* date 01/18/23 */
/*-----------------------------*/


capture program drop logitjack
program define logitjack, rclass
	syntax varlist(fv) [if] [in], CLUSter(varname) [ABSorb(varlist) VARiables(real 1) boots(real 1) LINear SAMple(string) NOnull ]
	
	
	local SMPL "`sample'"
	
	/*add if to sample*/
	if "`sample'" != "" {
	    local SMPL = "if `sample'"
	}
	
	/*flag for linearized or boot*/
	local flagboot = `boots'>1
	local flaglinear = 0 
	if "`linear'"!=""{
		local flaglinear = 1
	}
	local flagboth = `flagboot'+`flaglinear'
	
	
	/*sort the data by the clustering variable*/
		/*make a sequential-numeric cluster indicator*/
		
	tempvar temp_sample
	mark `temp_sample'
	markout `temp_sample' `varlist' `absorb' `cluster',strok
	
	/*impose the "sample option" restriction*/
	if "`sample'" != "" {
		
		tempvar temp_samp
		qui gen `temp_samp' = 0
		qui replace `temp_samp' = 1 `SMPL' 
		qui replace `temp_sample' = `temp_sample'*`temp_samp'
	}
	
			
		
	tempvar temp_indexa	 
	qui egen `temp_indexa' = group(`cluster') if `temp_sample'==1
	qui summ `temp_indexa'
	local G = r(max)
			
	qui sort `temp_indexa'
	qui putmata CVAR=`temp_indexa' if `temp_sample'==1, replace
	
	mata: numobs = rows(CVAR)
	mata: st_numscalar("numobs", numobs)
	mata: st_local("numobs", strofreal(numobs))	
		
	
	mata: info = panelsetup(CVAR,1)
	mata: ng = info[.,2] - info[.,1] :+ 1
	mata: G = rows(ng)
	

	
	
	if `variables' == 1{
		local x = word("`varlist'",2)
		local y = word("`varlist'",1)
	}
	
	*to impose null
	if "`nonull'" == "" {
		local newvarlist = subinstr("`varlist'","`x'","",1)
	}


	/*when absorb specified*/
	if "`absorb'" != ""{
		
		qui  levelsof(`absorb'), local(abslevels)
 
		local i =0
		
		local allti ""
		foreach abv in `abslevels' {
	
			local i = `i' +1
			local ti = "t_`i'"
			tempvar `ti'
			qui gen ``ti'' = `absorb'==`abv' 
			local allti = "`allti' ``ti''"
		}
		
		qui logit `varlist' `allti', nocons cluster(`cluster')
		
		
		local beta = _b[`x']
		mata: beta = `beta'
		
		local se =  _se[`x']
		mata: cv1se = `se'
		
		mata: tcv1 = beta/cv1se
		
		local yhat = "yhat"
		tempvar yhat
		qui predict `yhat' if `temp_sample'==1
		
		qui putmata yhat = `yhat' if `temp_sample'==1, replace
		
		if "`nonull'" == ""{ 
		
			qui logit `newvarlist' `allti', nocons cluster(`cluster')
					
			local yhatr = "yhatr"
			tempvar yhatr
			qui predict `yhatr' if `temp_sample'==1
			
			qui putmata yhatr = `yhatr' if `temp_sample'==1, replace
			
			/*dump everything to mata*/
			qui putmata ALLR = (`newvarlist' `allti') if `temp_sample'==1, replace
			
			mata numall = cols(ALLR)
			mata XR = ALLR[.,2::numall]
			
		}
		
		/*dump everything to mata*/
		qui putmata ALL = (`varlist' `allti') if `temp_sample'==1, replace
		
		mata numall = cols(ALL)
		mata Y = ALL[.,1]
		mata X = ALL[.,2::numall]
	
		
		/*generate the omit one strings*/
		*disp "i is `i' "
		mata: betas= J(`i',1,.)
		
		forvalues j=1/`i' {
			
			local omit_`j' = ""
			forvalues k = 1/10{
				
				if `k'!=`j'{
					local tk = "t_`k'"
					local omit_`j' = "`omit_`j'' ``tk'' "
				}
			}
			
			*disp "omit`j' is `omit_`j''"
		
		}/*end of j*/
		
		/*brute force jackknife*/

		if `flagboth' == 0 {
			local j = 0
			foreach abv in `abslevels' {
				local j = `j' + 1
				
				 qui logit `varlist' `omit_`j'' if `absorb'!=`abv', nocons cluster(`cluster')
				 local betas = _b[`x']
				 mata: betas[`j',1]= `betas'
				 
			} /*end of j*/
			
			mata: cv3 = ((`i'-1)/`i')*variance(betas:-beta)
			mata: cv3se = sqrt(cv3)
			mata: st_numscalar("cv3se", cv3se)
			mata: st_local("cv3se", strofreal(cv3se))
			
			disp "the cv3se is `cv3se'"
		}
		
	}/*end of absorb*/
	
	/*when absorb not specified*/
	if "`absorb'" == ""{
	
		*estimate logit 
		qui logit `varlist', cluster(`cluster')
		
		
		local beta = _b[`x']
		mata: beta = `beta'
		
		local se =  _se[`x']
		mata: cv1se = `se'
		
		mata: tcv1 = beta/cv1se
		
		
		local yhat = "yhat"
		tempvar yhat
		qui predict `yhat' if `temp_sample'==1
		
		qui putmata yhat = `yhat' if `temp_sample'==1, replace

		
		/*dump everything to mata*/
		qui putmata ALL = (`varlist') if `temp_sample'==1, replace
		
		*add a constant to X
		qui mata: ones = J(numobs,1,1)
		
		mata numall = cols(ALL)
		mata Y = ALL[.,1]
		mata X = (ALL[.,2::numall], ones)
		
		if "`nonull'" == "" {
			*estimate logit 
			qui logit `newvarlist', cluster(`cluster')
			
			local yhatr = "yhatr"
			tempvar yhatr
			qui predict `yhatr' if `temp_sample'==1
			
			qui putmata yhatr = `yhatr' if `temp_sample'==1, replace

			/*dump everything to mata*/
			qui putmata ALLR = (`newvarlist') if `temp_sample'==1, replace
			
			*add a constant to X
			qui mata: ones = J(numobs,1,1)
			
			mata numall = cols(ALLR)
			mata XR = (ALLR[.,2::numall], ones)
			
		} /*end of nonull */	
		
	}
	
	
	/*linearize and optionally bootstrap*/
	if `flagboth' != 0 {
		
		*linearize
		
		
		/*
			scores -equation 18
			\bis_g(\bbeta) = \sum_{i=1}^{N_g} \big(y_{gi} - 
			\Lambda(\biX_{gi}\bbeta)\big)\biX_{gi},\;\; g=1,\ldots,G,
		*/
		
		mata: scores = (Y :- yhat):*X
		mata: score_g =  panelsum(scores, info)
		
		if "`nonull'" == "" {
			
			mata: scoresr = (Y :- yhatr):*X
			mata: score_rg =  panelsum(scoresr, info)
		}

		
		/*
			contributions to info 
			\biJ_g(\bbeta) = \sum_{i=1}^{N_g} \Lambda_{gi}\tk(\bbeta)
			\Lambda_{gi}\tk(-\bbeta) \biX_{gi}^\top\biX_{gi}.
		*/
		
		*calculate ob level contribution
		timer clear 1
		timer on 1
		*local numobs = 2
		*local numobs = 25851

		*forvalues i = 1 / `numobs' {
		*	mata: infomat_`i' = yhat[`i',1]*(1 - yhat[`i',1]):*X[`i',.]'X[`i',.]
		*}
		
		*calculate cluster level contribution
		local G = 10 

		forvalues g = 1 / `G' {
			
			mata: firstrow = info[`g',1]
			mata: st_numscalar("firstrow", firstrow)
			mata: st_local("firstrow", strofreal(firstrow))
			mata: lastrow = info[`g',2]
			mata: st_numscalar("lastrow", lastrow)
			mata: st_local("lastrow", strofreal(lastrow))
			
			local f = `firstrow'
			mata: infomat_g`g' = yhat[`f',1]*(1 - yhat[`f',1]):*X[`f',.]'X[`f',.]
			
			if "`nonull'" == "" {
				mata: infomat_rg`g' = yhatr[`f',1]*(1 - yhatr[`f',1]):*X[`f',.]'X[`f',.]

			}
			
			local next = `firstrow'+1
			forvalues i = `next'/`lastrow'{
				
				mata: infomat_g`g' = infomat_g`g' + yhat[`i',1]*(1 - yhat[`i',1]):*X[`i',.]'X[`i',.]
				
				if "`nonull'" == "" {
				
					mata: infomat_rg`g' = infomat_rg`g' + yhatr[`i',1]*(1 - yhatr[`i',1]):*X[`i',.]'X[`i',.]
				
				}
				
			} /*end of i*/
			
			
		}/*end of g*/
		
		timer off 1 
		timer list 1

		local G =10
		*sum the score_g and info_g
		mata: score_all = colsum(score_g)
		mata: infomat_all = infomat_g1
		forvalues g = 2/`G'{
			mata: infomat_all = infomat_all + infomat_g`g'
		}
		
		mata: invinfomat_all = invsym(infomat_all)

		mata: linbeta = invinfomat_all*score_all'
		
		
		
		if "`nonull'" != ""{ 
			mata: score_us = J(`G',cols(X),.)
		}
		if "`nonull'" == ""{ 
			mata: score_rall = colsum(score_rg)
			mata: infomat_rall = infomat_rg1
			forvalues g = 2/`G'{
				mata: infomat_rall = infomat_rall + infomat_rg`g'
			}
			
			mata: invinfomat_rall = invsym(infomat_rall)

			mata: linbetar = invinfomat_rall*score_rall'
		
			
			mata: score_rs = J(`G',cols(X),.)
		
		} /*end of nonull*/
		
		
		*cv3l
		local G =10
		mata: cv3lsum = J(rows(linbeta),rows(linbeta),0)
		forvalues g = 1 / 10 {
			
			*unrestricted
			mata: beta_o`g' = invsym(infomat_all - infomat_g`g')*(score_all'-score_g[`g',.]')
			
			mata: cv3lsum = cv3lsum + (beta_o`g':-linbeta)*(beta_o`g'':-linbeta')
			
		
			
			*transform the scores
			if "`nonull'" != ""{ 
				mata: score_us[`g',.] = score_g[`g',.] - (infomat_g`g'* beta_o`g')' 
			}
			if "`nonull'" == ""{ 
					*restricted
				mata: beta_ro`g' = invsym(infomat_rall - infomat_rg`g')*(score_rall'-score_rg[`g',.]')
				mata: score_rs[`g',.] = score_g[`g',.] - (infomat_g`g'* beta_ro`g')' 
			}
			
			
		} /*end of g*/
		
		mata: cv3l = ((`G'-1)/`G'):*cv3lsum
		
		mata: sqrt(cv3l[1,1])
		
	} /*end of flag*/
	
	*bootstrap 
	if `boots' > 1 {
		
		timer  clear 2
		timer  on  2
		
		mata: w_sum = J(cols(score_g),cols(score_g),0)
		forvalues g = 1 / `G'{
				mata: w_g = score_g[`g',.]' - infomat_g`g' * linbeta
				mata: w_sum = w_sum + w_g*w_g'
		}
			
		*variance
		mata: var_cv1 = ((`G')/(`G'-1)):* invinfomat_all * w_sum * invinfomat_all
		
		mata: t = linbeta[1,1]/sqrt(var_cv1[1,1])
		
		
		
		mata: tboots = J(`boots',1,.)
		mata: tbootss = J(`boots',1,.)
		
		forvalues b = 1 / `boots' {
			
			*draw a G-vector of rademacher weights
			mata  e = runiform(`G',1)

			mata v = 2:* (e :> 0.5) :- 1
			
			/*unrestricted*/
			if "`nonull'" != "" {
				
				/*classic - C*/
				
				*transform the scores, then sum them			
				mata: sum_score_b = colsum(v:*score_g)
				
				*bootstrap beta
				mata beta_boot = invinfomat_all * sum_score_b'
				
				*bootstrap empirical scores
				mata: w_sum = J(cols(score_g),cols(score_g),0)
				forvalues g = 1 / `G'{
					mata: w_g = v[`g',1]:*score_g[`g',.]' - infomat_g`g' * beta_boot
					mata: w_sum = w_sum + w_g*w_g'
				}
				
				*bootstrap variance
				 mata: var_boot = ((`G')/(`G'-1)):* invinfomat_all * w_sum * invinfomat_all
				
				*bootstrap t 
				mata: tboots[`b',.] = beta_boot[1,1]/sqrt(var_boot[1,1])
				
				/*scores - S*/
				*transform the scores, then sum them			
				mata: sum_score_sb = colsum(v:*score_us)
				
				*bootstrap beta
				mata beta_boot = invinfomat_all * sum_score_sb'
				
				*bootstrap empirical scores
				mata: w_sum = J(cols(score_g),cols(score_g),0)
				forvalues g = 1 / `G'{
					mata: w_g = v[`g',1]:*score_us[`g',.]' - infomat_g`g' * beta_boot
					mata: w_sum = w_sum + w_g*w_g'
				}
				
				*bootstrap variance
				 mata: var_boot = ((`G')/(`G'-1)):* invinfomat_all * w_sum * invinfomat_all
				
				*bootstrap t 
				mata: tbootss[`b',.] = beta_boot[1,1]/sqrt(var_boot[1,1])
				
			
			} /* end of unrestricted */
			
			
			/*restricted*/
			if "`nonull'" == "" {
				
				/*classic - C*/
				*transform the scores, then sum them			
				mata: sum_score_b = colsum(v:*score_rg)
				
				*bootstrap beta
				mata beta_boot = invinfomat_rall * sum_score_b'
				
				*bootstrap empirical scores
				mata: w_sum = J(cols(score_g),cols(score_g),0)
				forvalues g = 1 / `G'{
					mata: w_g = v[`g',1]:*score_rg[`g',.]' - infomat_rg`g' * beta_boot
					mata: w_sum = w_sum + w_g*w_g'
				}
				
				*bootstrap variance
				 mata: var_boot = ((`G')/(`G'-1)):* invinfomat_rall * w_sum * invinfomat_rall
				
				*bootstrap t 
				mata: tboots[`b',.] = beta_boot[1,1]/sqrt(var_boot[1,1])
				
				/*score - S*/
				*transform the scores, then sum them			
				mata: sum_score_sb = colsum(v:*score_rs)
				
				*bootstrap beta
				mata beta_boot = invinfomat_rall * sum_score_sb'
				
				*bootstrap empirical scores
				mata: w_sum = J(cols(score_g),cols(score_g),0)
				forvalues g = 1 / `G'{
					mata: w_g = v[`g',1]:*score_rs[`g',.]' - infomat_rg`g' * beta_boot
					mata: w_sum = w_sum + w_g*w_g'
				}
				
				*bootstrap variance
				 mata: var_boot = ((`G')/(`G'-1)):* invinfomat_rall * w_sum * invinfomat_rall
				
				*bootstrap t 
				mata: tbootss[`b',.] = beta_boot[1,1]/sqrt(var_boot[1,1])
			
			} /* end of unrestricted */

			
			
			
		} /*end of b*/
		
		*bootstrap p
			
			*classic
			mata: bootc_p = mean(abs(tcv1[1,1]):<abs(tboots))
			
			mata: st_numscalar("bootc_p", bootc_p)
			mata: st_local("bootc_p", strofreal(bootc_p))
			
			disp "classic boot p is `bootc_p'"
			
			*scores
			mata: boots_p = mean(abs(tcv1[1,1]):<abs(tbootss))
			
			mata: st_numscalar("boots_p", boots_p)
			mata: st_local("boots_p", strofreal(boots_p))
			
			disp "score boot p is `boots_p'"
			
			timer  off 2 
			timer  list 2
		
	} /*end of boots if*/
	
	
	
	
	*jackknife: logit `varlist' i.`absorb', cluster(`cluster')
	
end