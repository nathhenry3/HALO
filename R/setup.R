### Setup ----
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(mrgsolve)
library(patchwork)

# Change colour palette for graphs
mycolors <- c("#08306B", "#2171B5", "#6BAED6", "#9ECAE1")

### ----

# Load in C++ model code for mrgsolve
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
    gamma_a = 0,
    lambda_b = 0,
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
    dxdt_apd = pow(apk, gamma_a) - k_apd * apd;
    dxdt_bpd = - lambda_b * pow(abs(bpk), gamma_b) - k_bpd * bpd;
    dxdt_H = apd + bpd - k_H * H;
  "

# Compile C++ code
mod <- mcode('Cppcode', cpp_code)