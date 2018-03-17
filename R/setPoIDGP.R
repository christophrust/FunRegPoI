
setPoIDGP <-
function(DGP_name){
    ## ########################
    ## A Function definining given a DGP_name the PoI Setup und the Beta(t) Setup
    ## Input:
    ##  - DGP_name: the Choice of PoI Setup
    ## - "Easy": 2 PoI at 0.3, 0.6, Betas -3, 3; beta(t) = -(t-1)^2+2 
    ## - "Complicated": 3 PoI at 0.3, 0.4, 0.6, Betas -3,-2, 3; beta(t) = 5*(t-0.5)^3-(t-0.5)+0.5
    ## - "Cubic": 3 PoI at 0.3, 0.4, 0.6, Betas -3,-2, 3; beta(t) = 5*(t-0.5)^3-(t-0.5)+0.5
    ## - "Sin": 3 PoI at 0.3, 0.4, 0.6, Betas -3, 2, 3; beta(t) = 0.3*sin(5*t)+0.5
    ## - "OnlyPoI": 3 PoI at 0.3, 0.4, 0.6, Betas -3, 2, 3; beta(t) = 0
    ## - "NoPoI": No PoIs, but beta(t)  = -(t-1)^2+2 
    ## Output:
    ##  - tau_ind_true: Indices of true tau 
    ##  - tau: Point at the intervall 
    ##  - beta_tau: Corresponding impact values
    ##  - fct_text: Beta(t) Function as text in t
	
    if( DGP_name == "Easy"){
        tau      <- c(0.3, 0.6) # PoI-locations tau_j
        beta_tau <- c(-3,3)    # PoI-parameters
        fct_text <- "-(t-1)^2+2"
        
    } else if( DGP_name == "Complicated"){
        tau      <- c(0.3, 0.4, 0.6) # PoI-locations tau_j
        beta_tau <- c(-3,3,3)    # PoI-parameters
        fct_text <- "5*(t-0.5)^3-(t-0.5)+0.5"
    } else if( DGP_name == "NoPoI"){
        tau      <- 0 # PoI-locations tau_j
        beta_tau <- 0  # PoI-parameters
        fct_text <- "-(t-1)^2+2"
    } else if( DGP_name == "OnlyPoI"){
        tau      <- c(0.3, 0.6) # PoI-locations tau_j
        beta_tau <- c(-3,3)    # PoI-parameters
        fct_text <- "0*t"	    
    }  else {
        stop(paste0("DGP_name: ", DGP_name, "could not be found")) 
    }
    
    ## Test for correct specification
    if(length(tau) != length(beta_tau)){
        stop("Taus do not match length of corresponding betas")
    }   

    list(tau = tau, beta_tau = beta_tau, fct_text=fct_text)
}
