### This file is just for testing behavioral stopping - i.e., before environment re-evaluation. 13 Sept 2023. 

### Setup
pacman::p_load(tidyverse, mrgsolve, patchwork, latex2exp) # You may need to install some of these packages from GitHub, and may require RTools. Refer to the documentation for remotes::install_github() and https://cran.r-project.org/bin/windows/Rtools/, respectively. 
## Package versions: 
# tidyverse = 2.0.0
# mrgsolve = 1.0.9
# patchwork = 1.1.2
# latex2exp = 0.9.6

# Change colour palette for graphs
myblues <- c("#08306B", "#2171B5", "#6BAED6", "#9ECAE1")

### -----

# # Create function to generate utility curve for prospect theory. x = vector of x-axis values. alpha = curvature of gains. lambda = loss aversion parameter. beta = curvature of losses. 
# prospect_utility <- function(x, alpha=0.2, lambda=1, beta=0.5) {
#   if (x >= 0) {
#     return(x^alpha)
#   } else {
#     return(-lambda * abs(x)^beta)
#   }
# }

# Load in C++ model code
cpp_code <- "
  $PARAM // Parameters for simulation
    // Clearance rates for compartments - reset by opponentprocess() function call
    k_Dose = 0,
    k_apk = 0,
    k_bpk = 0,
    k_apd = 0,
    k_bpd = 0,
    k_H = 0,
    
    // Pharmacodynamic constants - reset by opponentprocess() function call
    E0_a = 0,
    Emax_a = 0,
    EC50_a = 0,
    gamma_a = 0,
    E0_b = 0,
    Emax_b = 0,
    EC50_b = 0,
    gamma_b = 0, 
    
    // Infusion duration
    infuse = 1
  
  $CMT // Model compartments
    Dose, // Hormonal concentration following Digital Behavior
    apk, // a-process pharmacokinetics
    apd, // a-process pharmacodynamics
    bpk, // b-process pharmacokinetics
    bpd, // b-process pharmacodynamics
    H // Overall hedonic outcomes
  
  $MAIN // Set additional relationships
    D_Dose = infuse; // Sets the infusion duration for digital behavior compartment
  
  $ODE // Ordinary Differential Equations
    dxdt_Dose = - k_Dose * Dose;
    dxdt_apk = k_Dose * Dose - k_apk * apk;
    dxdt_bpk = k_apk * apk - k_bpk * bpk;
    dxdt_apd = Emax_a * pow(apk, gamma_a) - k_apd * apd;
    dxdt_bpd = Emax_b * pow(bpk, gamma_b) - k_bpd * bpd;
    dxdt_H = apd - bpd - k_H * H;
  "

# Compile C++ code
mod <- mcode('Cppcode', cpp_code)

### opponentprocess() takes arguments for PK/PD models, creates an mrgsolve compartmental model, and plots the output.
opponentprocess <- function(
    ii=100000, # Dosing interval
    sim_length=4000, # Time length of PKPD simulation, in minutes
    addl=10000, # Number of additional doses to deliver - essentially infinite.
    plot_utility=FALSE, # Whether to calculate the graphs for biophase or not.
    
    # Set PK/PD constants for C++ code
    k_Dose=10,
    k_apk=0.02,
    k_bpk=0.004,
    k_apd=1,
    k_bpd=1,
    k_H=1,
    E0_a=0,
    Emax_a=1,
    EC50_a=1,
    gamma_a=2,
    E0_b=0,
    Emax_b=1,
    EC50_b=3,
    gamma_b=2,
    
    # Set infusion duration for drug input
    infuse=1,
    
    ## Set values for a-process, b-process, and double gamma plot datasets to return. These values represent the number of doses that fall within the allocated timeframe 
    plot_2=c(0.0015, 0.006) # PKPD models for dose frequencies listed in plot_2 are plotted along with their associated Bode plot.
) { 
  
  # Create data frame of parameters to pass to simulation.
  idataset=data.frame(
    k_Dose=k_Dose,
    k_apk=k_apk,
    k_bpk=k_bpk,
    k_apd=k_apd,
    k_bpd=k_bpd,
    k_H=k_H,
    E0_a=E0_a,
    Emax_a=Emax_a,
    EC50_a=EC50_a,
    gamma_a=gamma_a,
    E0_b=E0_b,
    Emax_b=Emax_b,
    EC50_b=EC50_b,
    gamma_b=gamma_b,
    infuse=infuse
  ) %>% 
    rowid_to_column("ID") # Add column of IDs to start of data frame
  
  if (nrow(idataset) > 4) stop('Number of simulations must be 4 or less. Check idataset') # Stop if number of simulations > 4 
  
  # Print out simulation parameters once
  if (plot_utility==TRUE) {
    cat('\nSimulation parameters =\n\n') 
    print(idataset)
  }  
  
  # Create a list of events
  events <- ev(amt = 1, # Dose amount
               rate = -2, # Signals that duration of infusion is modeled
               ii = ii, # Dosing interval
               ID = 1:nrow(idataset), # Add number of simulations being run
               addl = addl) # No. of additional doses to administer
  
  # Run model
  out <- mrgsim(mod, events, idataset, end=sim_length, maxsteps=50000)
  
  # Calculate integral (AUC, area under curve) of H (hedonic) compartment
  AUC_H <- out@data %>% 
    group_by(ID) %>% 
    summarise(AUC=sum(H))
  
  # Calculate dose frequency
  freq <- 1/ii 
  AUC_H$freq <- freq
  
  # Print results
  cat(paste('Integral of hedonic graph for simulation', AUC_H$ID, '=', AUC_H$AUC, '\n'))
  cat(paste('Dose frequency =', freq, 'per min\n\n'))
  
  # If rounded dose frequency value falls within plot_2 list, then return plot of H compartment
  if (isTRUE(all.equal(freq, plot_2[[1]])) | isTRUE(all.equal(freq, plot_2[[2]]))) { # Use all.equal() to check equivalence of floating point numbers
    
    # Plot of H compartment
    cat('Saving plots for dose frequency above.................\n\n')
    plot_2_freq <- out@data %>%
      ggplot(aes(x=time, y=H, colour=factor(ID))) +
      geom_hline(yintercept=0, linetype='dashed', color='black') +
      geom_line() +
      scale_color_manual(values=myblues) +
      ggtitle(bquote(paste('Dose frequency = ', .(freq)) ~ min^-1)) +
      xlab('Time, t [min]') + {
        if (isTRUE(all.equal(freq, plot_2[[1]]))) ylab(bquote(paste('Hedonic scale, H'[a*','*b])))# Only create y label if first plot
      } +
      theme_light() + {
        if (isTRUE(all.equal(freq, plot_2[[1]]))) { # Only create y label if first plot
          theme(plot.title=element_text(size=9, hjust=0.5, margin=margin(t=0, b=0)),
                legend.position='none')
        } else {
          theme(plot.title=element_text(size=9, hjust=0.5, margin=margin(t=0, b=0)),
                legend.position='none',
                axis.title.y=element_blank())
        }
      }
  } else {
    plot_2_freq <- NULL
  }
  
  ## Create plots for biophase curves for PK -> PD conversion, using biophase equations
  
  if (plot_utility == TRUE) {
    # Set x axis length with dose_seq, then calculate biophase curves
    utility_data <- tibble(x=seq(-20, 20, 0.5))
    
    for (i in 1:nrow(idataset)) { # Calculate biophase curve for each set of parameters
      # Create column name for biophase curve based on ID number
      apd_colname <- paste0('apd', i)
      
      # Calculate biophase curves
      utility_data <- utility_data %>% 
        mutate({{apd_colname}} := case_when(
          x >= 0 ~ idataset$Emax_a[i] * x ^ idataset$gamma_a[i],
          x < 0 ~ -idataset$Emax_b[i] * abs(x) ^ idataset$gamma_b[i])
        )
    }

    # Plots for biophase curves
    utility_graph <- utility_data %>% 
      pivot_longer(cols=starts_with('apd'), names_to='ID', values_to='Values') %>% 
      ggplot(aes(x=x, y=Values, colour=ID)) +
      geom_line() + 
      geom_hline(yintercept=0, linetype='dashed') +
      geom_vline(xintercept=0, linetype='dashed') +
      scale_color_manual(values=myblues) +
      theme_light() +
      theme(plot.title=element_text(size=9, hjust=0.5),
            legend.position='none') +
      ylab('Value') +
      xlab('Outcome (losses vs gains)')
  }
  
  ## ------------------------------------------------------------------------------
  
  # Return necessary objects
  ifelse(plot_utility == TRUE,
         return(list(AUC_H, freq, plot_2_freq, utility_graph)),
         return(list(AUC_H, freq, plot_2_freq)))
}

## Example run of function above:
# opponentprocess(EC50_b=c(8, 6, 4.5, 3))

### bode_plot() takes opponentprocess() and runs it across a range of dose frequencies, thus allowing us to plot the relationship between dose frequency and the integral of hedonic outcomes, and to determine whether this relationship is hormetic. 
bode_plot <- function(
    # Pass on arguments to opponentprocess()
  ..., 
  
  # Set x values for biophase graphs
  seq_1=0.00025,
  seq_2=0.006,
  
  # Set y limit for hormesis graph (integer). If NA, ylim is automatically set
  gg_ylim=NA
) {
  
  # List of dose intervals to pass to opponentprocess()
  dose_interval <- c(0,  seq(seq_1, seq_2, seq_1)^-1) 
  
  # Run loop to calculate Bode magnitude plot across range of frequencies
  H_list <- list() # Create list to house graphs of hedonic outcomes vs time
  for (i in 1:length(dose_interval)) {
    
    # If first dose interval, then set up bode_data data frame
    if (i == 1) {
      
      loop_list <- opponentprocess(ii=dose_interval[2], 
                                   ...)
      
      # Create data frame to store wellbeing scores in, based on number of simulations performed
      bode_data <- tibble(
        'ID' = 1:nrow(loop_list[[1]]), # ID of each mrgsolve simulation
        'AUC' = rep(0, nrow(loop_list[[1]])), # AUC scores for H compartment graphs
        'freq' = rep(0, nrow(loop_list[[1]])) # Dose frequency
      )
    } else if (i == 2) {
      # If second dose interval, then calculate biophase graphs. Otherwise just calculate loop_list to append to bode_data
      loop_list <- opponentprocess(ii=dose_interval[i],
                                   plot_utility=TRUE,
                                   ...)
      utility_plot <- loop_list[[4]]
    } else {
      loop_list <- opponentprocess(ii=dose_interval[i],
                                   ...)
    }
    
    # Append AUC scores (hedonic outcomes) and dose frequencies, and store H compartment graphs
    bode_data <- bode_data %>% 
      rbind(loop_list[[1]])
    H_list[[length(H_list) + 1]] <- loop_list[[3]] # H compartment simulations
  }
  
  # Create Bode magnitude plot
  bode_graph <- bode_data %>% 
    ggplot(aes(x=freq, y=AUC, colour=factor(ID))) +
    geom_hline(yintercept=0, linetype='dashed', color='black') + 
    geom_line() +
    scale_color_manual(values=myblues) + {
      if (!is.na(gg_ylim)) {
        coord_cartesian(ylim=c(gg_ylim, NA))
      }
    } +
    xlab(bquote('Dose frequency, f' ~ '[' * min^-1 * ']')) +
    ylab(bquote(integral(H[a*','*b](t)[total]*dt, 0, t[sim]))) +
    theme_light() +
    theme(legend.position='none')
  
  # Patch pharmacodynamic, temporal, and Bode plots  together with patchwork package, and print
  tryCatch(bode_patch <- utility_plot / (first(H_list) | last(H_list)) / bode_graph,
           error=function(e) {
             warning(e)
             stop("Error: bode_patch didn't patch together correctly. Double-check that the values in plot_2 are exact multiples of seq_1.")
           }
  )
  suppressWarnings(print(bode_patch))
}

###
### EXAMPLES IN ARTICLE:

# bode_plot() # Default hormesis curve by itself.
# bode_plot(EC50_b=c(1.8, 2.4, 3, 3.6), seq_1=0.00025, sim_length=4000, gg_ylim=-285) # Example of modifying EC50_b. Figure 4 in article.
# bode_plot(seq_1=0.005, seq_2=0.2, plot_2=c(0.005, 0.2), EC50_a=3, EC50_b=35, k_apk=c(0.02, 0.04), k_bpk = c(0.005, 0.004), addl=30, sim_length=3000) # Replication of Solomon & Corbit's affective dynamics graph. Figure 5 in article. Note how multiple parameter combinations can be run simultaneously.
# 
# bode_plot(seq_1=0.005, seq_2=0.2, plot_2=c(0.005, 0.2), EC50_a=3, EC50_b=35, k_apk=c(0.02, 0.04), k_bpk = c(0.005, 0.004), addl=60, sim_length=10000) # Stopping behaviors before re-evaluation occurs

bode_plot(gamma_a=0.2, gamma_b=0.5, Emax_a = 1, Emax_b = 1, k_apk = 0.005, k_bpk = 0.004, seq_1 = 0.0006, seq_2=0.006, plot_2=c(0.0006, 0.006))









### Test some parameter variations

# bode_plot() # Default hormesis curve by itself
# bode_plot(seq_1=0.0002) # Figure 3

# bode_plot(EC50_b=c(1.8, 2.4, 3, 3.6), seq_1=0.00025, sim_length=4000, gg_ylim=-285) # Modifying EC50_b. Figure 4

# bode_plot(seq_1=0.002, seq_2=0.1, plot_2=c(0.002, 0.1), EC50_a=c(1.5,3.5), EC50_b=20, addl=200) # Try instead by moderating EC50A to reproduce S&C's tolerance development. Figure 5 (deprecated)

# bode_plot(seq_1=0.005, seq_2=0.1, plot_2=c(0.005, 0.1), EC50_a=c(2.2, 3.2), EC50_b=c(15, 12), k_bpk = c(0.006, 0.003), addl=200) # S&C replication attempt 1
# bode_plot(seq_1=0.005, seq_2=0.1, plot_2=c(0.005, 0.1), EC50_a=c(2, 3.4), EC50_b=17, k_bpk = c(0.006, 0.004), addl=150, sim_length=3000) # S&C replication attempt 1
# bode_plot(seq_1=0.005, seq_2=0.2, plot_2=c(0.005, 0.2), EC50_a=2, EC50_b=25, k_apk=c(0.015, 0.04), k_bpk = c(0.0045, 0.003), addl=300, sim_length=3000) # S&C replication attempt 1
# bode_plot(seq_1=0.005, seq_2=0.2, plot_2=c(0.005, 0.2), EC50_a=3, EC50_b=35, k_apk=c(0.015, 0.04), k_bpk = c(0.0045, 0.003), addl=300, sim_length=3000) # S&C replication attempt 1
# bode_plot(seq_1=0.005, seq_2=0.2, plot_2=c(0.005, 0.2), EC50_a=3, EC50_b=40, k_apk=c(0.015, 0.04), k_bpk = c(0.004, 0.003), addl=300, sim_length=3000) # S&C replication attempt 1
# bode_plot(seq_1=0.005, seq_2=0.2, plot_2=c(0.005, 0.2), EC50_a=3, EC50_b=35, k_apk=c(0.02, 0.04), k_bpk = c(0.005, 0.0035), addl=300, sim_length=3000) # Replication attempt of Solomon & Corbit affective dynamics graph. Figure 5.

# bode_plot(EC50_a=c(1, 1.4), seq_1=0.005, seq_2=0.7, plot_2=c(0.02, 0.08), addl=50) # Drugs cut off after 100 doses
# bode_plot(EC50_b=c(6, 4), seq_1=0.005, seq_2=1, plot_2=c(0.02, 0.5), addl=750, infuse=1) # Infusion time set to virtually nil, and EC50_b modified
# bode_plot(EC50_a=c(1, 1.4), seq_1=0.005, seq_2=1, plot_2=c(0.01, 0.5), addl=750, infuse=1) # Infusion time set to virtually nil, and EC50_a modified
# 
# ## Trying multiple parameters at once
# 
# bode_plot(k_apk=c(0.1, 0.05, 0.005), k_bpk=c(0.004, 0.002, 0.0002), seq_1=0.05, seq_2=1, sim_length=25000, plot_2=c(0.05, 0.85), addl=750, infuse=1) # Infusion time set to virtually nil, and k_apk and k_bpk modified, with longer simulation length
# bode_plot(k_apk=c(0.1, 0.02), Emax_b=c(0, 14.5), k_bpk=c(0.007, 0.007), seq_1=0.02, seq_2=1, plot_2=c(0.02, 1), addl=400, infuse=1, sim_length=1250) # Infusion time set to virtually nil, and modifying k_apk and Emax_b such that we can replicate Paalzow & Sjalander-Brynne, 1998 - dopamine model
# bode_plot(infuse=c(20, 40, 80)) # Note also that increasing infusion time decreases hormesis!!!
# 
# # Increase dose frequency and remove dose to generate withdrawal
# 
# bode_plot(seq_1=0.1, seq_2=0.15, plot_2=c(0.1, 0.15), EC50_b=20) # Have to modify EC50_b to match opponent process
# bode_plot(seq_1=0.005, seq_2=0.1, plot_2=c(0.005, 0.1), EC50_b=c(3,20), addl=200) # Matches opponent process figure in Solomon & Corbit's model perfectly!
# bode_plot(seq_1=0.005, seq_2=0.1, plot_2=c(0.005, 0.1), EC50_a=c(1,3.5), EC50_b=20, addl=200) # Try instead by moderating EC50A to reproduce S&C's tolerance development.
# bode_plot(seq_1=0.5, seq_2=1, plot_2=c(0.5, 1), Emax_a=c(3,1.5), EC50_b=20, addl=2000) # Try focus on getting infinite dose
# bode_plot(seq_1=0.1, seq_2=1, plot_2=c(0.1, 1), EC50_b=20, addl=2000) # At even higher concentrations, a-process and b-process reach equilibrium at hedonic scale = 0. 
# 
# ## Let's try increasing the timescale of these plots
# 
# bode_plot(seq_2=0.15, plot_2=c(0.005, 0.15)) # Crazy - turns out that hormetic relationship is actually a triphasic relationship! Beyond a certain dose frequency, the A process dominates if it isn't given time to recover. E.g. dopamine - acute exposure is positive, moderately high frequency of exposure is negative, but continuous exposure is positive again!!!
# bode_plot(seq_1=0.05, seq_2=1, plot_2=c(0.05, 1)) # If you extend the timescale far enough, you get convergence to some value of integral(H(t))
# bode_plot(seq_1=0.05, seq_2=1, plot_2=c(0.05, 1), k_DB=c(0.2, 1)) # If you decrease the dopamine decay constant (i.e. longer half life), then you will find that dopamine equilibrium is reached much faster!
# bode_plot(seq_1=0.05, seq_2=1, plot_2=c(0.05, 1), k_DB=5) # But what if you INCREASE the dopamine decay constant? You get negative hormesis...
# bode_plot(seq_1=0.002, seq_2=0.05, plot_2=c(0.002, 0.05), k_DB=5) # This explains the negative curve I saw before.
# 
# 
# # Hill equation with gamma parameter of 2 represents a departure from Michaelis-menten kinetics, with positive cooperativity since gamma > 1. What if we try changing cooperativity?
# bode_plot(seq_2=0.05, plot_2=c(0.005, 0.05), gamma_a=1, gamma_b=1) # No triphasic relationship
# bode_plot(seq_2=0.05, plot_2=c(0.005, 0.05), gamma_a=1.5, gamma_b=1.5) # Closer to triphasic relationship but still monopositive
# bode_plot(seq_2=0.05, plot_2=c(0.005, 0.05), gamma_a=1.5, gamma_b=1.5, EC50_b=3) # Try reducing EC50_b also to increase b-process: get negative hormesis.
# bode_plot(seq_2=0.05, plot_2=c(0.005, 0.05), gamma_a=1, gamma_b=2) # Try only setting gamma_b to 2 - monopositive increase.
# bode_plot(seq_2=0.05, plot_2=c(0.005, 0.05), gamma_a=1, gamma_b=2, EC50_b=2) # Back to triphasic hormesis
# bode_plot(seq_1=0.0001, seq_2=0.05, plot_2=c(0.005, 0.05), gamma_a=0.5, gamma_b=1, EC50_b=2) # Try negative cooperativity: monopositive trend.
# 
# ## So ultimately for hormesis to work, you need a high gamma (close to 2 or above) for the b-process.

## Unusual modifications: [NOTE: need to add k_apd, k_bpd, k_H]

# bode_plot(seq_1=0.005, seq_2=0.2, plot_2=c(0.005, 0.2), EC50_a=3, EC50_b=35, k_apk=c(0.02, 0.025), k_bpk = c(0.005, 0.004), k_apd=c(1, 1.1), k_bpd=c(1, 0.9), addl=300, sim_length=3000)
# bode_plot(seq_1=0.005, seq_2=0.2, plot_2=c(0.005, 0.2), EC50_a=3, EC50_b=35, k_apk=c(0.02, 0.04), k_bpk = c(0.005, 0.004), addl=300, -sim_length=3000)
# bode_plot(seq_1=0.005, seq_2=0.2, plot_2=c(0.005, 0.2), EC50_a=3, EC50_b=35, k_apk=c(0.02, 0.04), k_bpk=c(0.005, 0.003), k_apd=c(1, 1.5), k_bpd=c(1, 1.8), addl=300, sim_length=3000)

## Working on Behavioral Posology for AI paper: trying some variations
# bode_plot(EC50_b=c(3, 5), k_apk=1000000, k_bpk=c(0.004, 0.0016), seq_1=0.00025, sim_length=4000, gg_ylim=-285) # Remove a-process, and just look at how modifying half-life vs modifying intensity affects the hormetic curve. Looks like half-life has a much greater effect, which makes sense. 

# Testing out JUST PK alterations. Can STILL ACHIEVE HORMESIS??
# bode_plot(EC50_b = 1, k_bpk = 0.02, seq_2 = 0.039, plot_2 = c(0.0015, 0.039))
