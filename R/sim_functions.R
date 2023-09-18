### This file contains functions for simulating repeated opponent processes using the behavioral posology paradigm, and for using Bode plots to analyze the healthy hormetic limits of behaviors. 

### ----

#' opponentprocess() simulates the opponent process model. It can be run by itself for diagnostic purposes (for example, calculating the area under the hedonic curve, or plotting either the opponent processes or utility function), but in general use this should be called through bode_plot(). 
#'
#' @param ii Dosing interval.
#' @param sim_length Time length of PKPD simulation, in minutes.
#' @param addl Number of additional doses to deliver - essentially infinite.
#' @param plot_utility Whether to calculate the graphs for biophase or not.
#' @param plot_op Whether to plot the graph of opponent processes (Hedonic state, H compartment)
#' @param join_plots Whether to join plots as subplots or plot them separately.
#'
#' @param k_Dose Clearance rate for compartments.
#' @param k_apk a-process pharmacokinetic clearance rate.
#' @param k_bpk b-process pharmacokinetic clearance rate.
#' @param k_apd a-process pharmacodynamic clearance rate.
#' @param k_bpd b-process pharmacodynamic clearance rate.
#' @param k_H Clearance rate for hedonic pharmacodynamic compartment.
#' @param lambda_a Pharmacodynamic constant for a-process magnitude, based on Kahneman/Tversky's utility function.
#' @param gamma_a Pharmacodynamic constant for a-process curvature, based on Kahneman/Tversky's utility function.
#' @param lambda_b Pharmacodynamic constant for b-process magnitude, based on Kahneman/Tversky's utility function.
#' @param gamma_b Pharmacodynamic constant for b-process curvature, based on Kahneman/Tversky's utility function.
#' @param infuse Infusion duration.
#' @param plot_frequencies # If the dose frequency (or frequencies) are in plot_frequencies, then create plots for those frequencies. Should be a vector of length 2 - will always plot
#' @param verbose # Print output to console.
#'
#' @return A list containing the integral of hedonic outcomes, dose frequency, and optional plots.
#'
#' @examples
#' opponentprocess() # Default parameters
#' opponentprocess(ii=10, sim_length=4000, addl=10000, plot_utility=FALSE, join_plots=TRUE, k_Dose=1, k_apk=0.01, k_bpk=0.01, k_apd=1, k_bpd=0.01, k_H=1, lambda_a=1, gamma_a=0.5, lambda_b=1, gamma_b=0.7, infuse=1, plot_frequencies=c(0.002), verbose=TRUE) # Default parameters
#' opponentprocess(plot_utility=TRUE) # Plots utility function
#' opponentprocess(plot_op=TRUE) # Plots hedonic compartment values (opponent processes)
#' 
#' @export
opponentprocess <- function(
    ii=10000, # Dosing interval
    sim_length=4000, # Time length of PKPD simulation, in minutes
    addl=10000, # Number of additional doses to deliver - essentially infinite.
    plot_utility=FALSE, # Whether to calculate the graphs for biophase or not.
    plot_op=FALSE, # Whether to plot the graph of opponent processes
    join_plots=TRUE, # Whether to join plots as subplots or plot them separately
    
    # Set PK/PD constants for C++ code
    k_Dose=1, # Higher k_Dose leads to more instantaneous dose release; smaller k_Dose smooths out the curves
    k_apk=0.01,
    k_bpk=0.01,
    k_apd=1,
    k_bpd=1,
    k_H=1,
    lambda_a=1,
    gamma_a=0.5,
    lambda_b=1,
    gamma_b=0.7,
    
    # Set infusion duration for drug input
    infuse=1,
    
    # If the dose frequency (or frequencies) are in plot_frequencies, then create plots for those frequencies. 
    plot_frequencies=c(0.002, 0.008), 
    
    # If TRUE, print output on state of simulation
    verbose=TRUE
) { 
  # Convert plot_frequencies to vector if it is numeric
  if (length(plot_frequencies) == 1) {
    plot_frequencies <- rep(plot_frequencies, 2)
    plot_second <- FALSE # Set flag to prevent plotting H compartment twice later on
  }
  
  # Create data frame of parameters to pass to simulation.
  idataset=data.frame(
    k_Dose=k_Dose,
    k_apk=k_apk,
    k_bpk=k_bpk,
    k_apd=k_apd,
    k_bpd=k_bpd,
    k_H=k_H,
    lambda_a=lambda_a,
    gamma_a=gamma_a,
    lambda_b=lambda_b,
    gamma_b=gamma_b,
    infuse=infuse
  ) %>% 
    rowid_to_column("ID") # Add column of IDs to start of data frame
  
  if (nrow(idataset) > 4) stop('Number of simulation variants must be 4 or less.') # Stop if number of simulations > 4 
  
  # Print out simulation parameters once
  if (plot_utility & verbose) {
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
  if (verbose) {
    cat(paste('Integral of hedonic compartment values during simulation', AUC_H$ID, '=', AUC_H$AUC, '\n'))
    cat(paste('Dose frequency =', freq, 'per min\n\n'))
  }

  # If rounded dose frequency value falls within plot_frequencies list, then return plot of H compartment
  if (plot_op | isTRUE(all.equal(freq, plot_frequencies[[1]])) | isTRUE(all.equal(freq, plot_frequencies[[2]]))) { # Check if plot frequency should be plotted
    
    # Plot of H compartment
    if (verbose) cat('Saving plots for dose frequency above.................\n\n')
    plot_hedonic <- out@data %>%
      ggplot(aes(x=time, y=H, colour=factor(ID))) +
      geom_hline(yintercept=0, linetype='dashed', color='black') +
      geom_line() +
      scale_color_manual(values=mycolors) +
      ggtitle(bquote(paste('Dose frequency = ', .(freq)) ~ min^-1)) +
      xlab('Time, t [min]') + {
        if (isTRUE(all.equal(freq, plot_frequencies[[1]])) | join_plots == FALSE) ylab(bquote(paste('Hedonic state, H'[a*','*b]*(t), ' [arb. units]'))) # Only create y label if first plot
      } +
      theme_light() + {
        if (isTRUE(all.equal(freq, plot_frequencies[[1]])) | join_plots == FALSE) { # Only create y label if first plot
          theme(plot.title=element_text(size=9, hjust=0.5, margin=margin(t=0, b=0)),
                legend.position='none')
        } else {
          theme(plot.title=element_text(size=9, hjust=0.5, margin=margin(t=0, b=0)),
                legend.position='none',
                axis.title.y=element_blank())
        }
      }
    
    # Plot the graph of opponent processes (hedonic state, H compartment)
    if (plot_op) print(plot_hedonic)
  } else {
    plot_hedonic <- NULL
  }
  
  ## Create plots for biophase curves for PK -> PD conversion, using biophase equations
  
  if (plot_utility) {
    # Set x axis length with dose_seq, then calculate biophase curves
    utility_data <- tibble(x=seq(-20, 20, 0.1))
    
    for (i in 1:nrow(idataset)) { # Calculate biophase curve for each set of parameters
      # Create column name for biophase curve based on ID number
      apd_colname <- paste0('apd', i)
      
      # Calculate biophase curves
      utility_data <- utility_data %>% 
        mutate({{apd_colname}} := case_when(
          x >= 0 ~ idataset$lambda_a[i] * x ^ idataset$gamma_a[i],
          x < 0 ~ -idataset$lambda_b[i] * abs(x) ^ idataset$gamma_b[i])
        )
    }

    # Plots for biophase curves
    utility_graph <- utility_data %>% 
      pivot_longer(cols=starts_with('apd'), names_to='ID', values_to='Values') %>% 
      ggplot(aes(x=x, y=Values, colour=ID)) +
      geom_line() + 
      geom_hline(yintercept=0, linetype='dashed') +
      geom_vline(xintercept=0, linetype='dashed') +
      scale_color_manual(values=mycolors) +
      theme_light() +
      theme(plot.title=element_text(size=9, hjust=0.5),
            legend.position='none') +
      geom_text(aes(x = 17, y = -2, label = "Gain"), hjust = 0) +
      geom_text(aes(x = -17, y = -2, label = "Loss"), hjust = 1) +
      ylab('Value [arb. units]') +
      xlab('Outcome [arb. units]')
  }
  
  ## ------------------------------------------------------------------------------
  
  # Return necessary objects
  ifelse(plot_utility,
         return(list(AUC_H, freq, plot_hedonic, utility_graph)),
         return(list(AUC_H, freq, plot_hedonic)))
}

### ----

#' bode_plot() creates a Bode plot using opponentprocess() for a range of dose frequencies, allowing us to plot the relationship between dose frequency and the integral of hedonic outcomes, and to determine whether this relationship is hormetic. 
#'
#' @param freq_interval The interval between dose frequencies.
#' @param multiply A multiplier for the dose frequencies.
#' @param gg_ylim Y-axis limit for the Bode graph (optional).
#' @param join_plots Whether to join graphs as subplots or plot them separately.
#'
#' @param ii Dosing interval for opponentprocess().
#' @param sim_length Time length of PKPD simulation for opponentprocess(), in minutes.
#' @param addl Number of additional doses for opponentprocess() - essentially infinite.
#' @param plot_utility Whether to calculate the graphs for biophase for opponentprocess() or not.
#' @param join_plots_opponent Whether to join plots as subplots or plot them separately in opponentprocess().
#' @param plot_frequencies # If the dose frequency (or frequencies) are in plot_frequencies, then create plots for those frequencies.
#'
#' @param k_Dose Parameter for opponentprocess(): Clearance rates for compartments.
#' @param k_apk Parameter for opponentprocess(): a-process pharmacokinetic clearance rate.
#' @param k_bpk Parameter for opponentprocess(): b-process pharmacokinetic clearance rate.
#' @param k_apd Parameter for opponentprocess(): a-process pharmacodynamic clearance rate.
#' @param k_bpd Parameter for opponentprocess(): b-process pharmacodynamic clearance rate.
#' @param k_H Parameter for opponentprocess(): Clearance rate for hedonic pharmacodynamic compartment.
#' @param lambda_a Parameter for opponentprocess(): Pharmacodynamic constant for a-process magnitude, based on Kahneman/Tversky's utility function.
#' @param gamma_a Parameter for opponentprocess(): Pharmacodynamic constant for a-process curvature, based on Kahneman/Tversky's utility function.
#' @param lambda_b Parameter for opponentprocess(): Pharmacodynamic constant for b-process magnitude, based on Kahneman/Tversky's utility function.
#' @param gamma_b Parameter for opponentprocess(): Pharmacodynamic constant for b-process curvature, based on Kahneman/Tversky's utility function.
#' @param infuse Parameter for opponentprocess(): Infusion duration.
#' @param verbose # Print output to console.
#'
#' @return A Bode magnitude plot.
#'
#' @examples
#' bode_plot()
#' bode_plot(gamma_a=0.5, gamma_b=c(0.5, 0.7, 0.9, 1.1), lambda_a=1, lambda_b = 1, k_apk = 0.005, k_bpk = 0.004, freq_interval = 0.0002, multiply=40, plot_frequencies=c(0.0002, 0.006))
#' bode_plot(ylim=)
#' 
#' @export
bode_plot <- function(
  # Pass on arguments to opponentprocess()
  ..., 
  
  # Set x values for biophase graphs
  freq_interval=0.0002,
  multiply=150,
  
  # If the dose frequency (or frequencies) are in plot_frequencies, then create plots for those frequencies
  plot_frequencies=c(0.002, 0.008),
  
  # Set lower y limit for Bode plot (integer). If NA, ylim is automatically set
  gg_ylim=NA,
  
  # Join graphs as subplots if TRUE, or plot separately if FALSE
  join_plots=TRUE,
  
  # If TRUE, print output on the state of the simulation
  verbose=TRUE
) {
  
  # List of dose intervals to pass to opponentprocess()
  dose_interval <- c(0,  seq(freq_interval, freq_interval*multiply, freq_interval)^-1) 
  
  # Run loop to calculate Bode magnitude plot across range of frequencies
  H_list <- list() # Create list to house graphs of hedonic outcomes vs time
  for (i in 1:length(dose_interval)) {
    
    # If first dose interval, then set up bode_data data frame
    if (i == 1) {
      
      loop_list <- opponentprocess(ii=dose_interval[2], 
                                   join_plots=join_plots,
                                   verbose=verbose,
                                   plot_frequencies=plot_frequencies,
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
                                   join_plots=join_plots,
                                   verbose=verbose,
                                   plot_frequencies=plot_frequencies,
                                   ...)
      utility_plot <- loop_list[[4]]
    } else {
      loop_list <- opponentprocess(ii=dose_interval[i],
                                   join_plots=join_plots,
                                   verbose=verbose,
                                   plot_frequencies=plot_frequencies,
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
    scale_color_manual(values=mycolors) + {
      if (!is.na(gg_ylim)) {
        coord_cartesian(ylim=c(gg_ylim, NA))
      }
    } +
    xlab(bquote('Dose frequency, f' ~ '[' * min^-1 * ']')) +
    ylab(bquote(integral(H[a*','*b](t)*dt, 0, t[sim]) ~ '[arb. units]')) +
    theme_light() +
    theme(legend.position='none')
  
  if (join_plots) {
    # Patch pharmacodynamic, temporal, and Bode plots  together with patchwork package, and print
    tryCatch(bode_patch <- utility_plot / (first(H_list) | last(H_list)) / bode_graph,
             error=function(e) {
               warning(e)
               stop("Error: graphs didn't patch together correctly. Check your inputs to bode_plot().")
             }
    )
    suppressWarnings(print(bode_patch))
  } else {
    # Plot all graphs individually
    tryCatch({
        print(utility_plot); print(first(H_list)); print(last(H_list)); print(bode_graph);
        return(invisible(list(utility_plot, first(H_list), last(H_list), bode_graph)))
      },
             error=function(e) {
               warning(e)
               stop("Error: plotting individual graphs failed. Check your inputs to bode_plot().")
             }
    )
  }
}





### TO DO: ----
# Add better comments
# Put cppcode in src folder?? For compiled codes
# Give it an option to set the color scheme
# Fix comments for 'multiply' parameter - perhaps replace with something else?
# Add roxygen comments for functions
# Come up with some examples
# Create a README.RMD file, which can then be knitted to MD to put on the front page of Github, and can be used to create examples for your article
# Write some test functions??

# THEN...
# Start here: https://ourcodingclub.github.io/tutorials/writing-r-package/
# Convert to R package project structure, or create a new project in that structure
# THEN render your page as a website!
# Once all done, use chatgpt to add comments for roxygen etc
# Also put on OSF