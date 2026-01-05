## ----------------------------------------------------------------------------
# Author: Julia Mayer
# Last modified: 23.05.2025

# This builds an MSR model in odin that can then be fitted to seroprevalence data.
# At birth, some individuals are protected from infection by maternal immunity (M).
# They become susceptible to infection as infection wanes with a rate μ (S).
# Individuals will seroconvert according to a FOI λ (R).
# The time spent in the M compartment (temporary protection) follows an Erlang 
# distribution with mean D_imm = 1/μ and shape equal to 2 and is modelled by 2 
# sub-stages M1 and M2 each being exponentially distributed with mean equal to 
# D_imm/2. The transition from M1 to M2 therefore happens at a rate 2*μ.
# (See https://sbfnk.github.io/mfiidd/play_with_seitl.html for more information 
# on modelling an exponential vs Erlang distribution.)

# We will use this model to estimate two parameters λ and the proportion of children
# born with maternal immunity p but you can adapt it to estimates other values 
# as well, for example μ.
# We will fit this model to RSV seroconversion data but it could be fitted to 
# any other type of data.
## -----------------------------------------------------------------------------
## Core equations for transitions between compartments:
update(M1) <- M1 - n_M1M2
update(M2) <- M2 + n_M1M2 - n_M2S
update(S) <- S + n_M2S - n_SR
update(R) <- R + n_SR

## Individual probabilities of transition:
n_M1M2 <- mu*2 * dt * M1 # M1 to M2
n_M2S <- mu*2 * dt * M2 # M2 to S
n_SR <- lambda * dt * S # S to R

## Initial states:
initial(M1) <- prop * M1_ini # only a proportion of children is born with maternal immunity
initial(M2) <- M2_ini
initial(S) <- (1-prop) * M1_ini
initial(R) <- R_ini

## User defined parameters - default in parentheses but it can be changed when 
# initialising the model:
M1_ini <- parameter(1 - 2 * 1e-12)
M2_ini <- parameter(1e-12)
R_ini <- parameter(1e-12)

# Transition parameters
mu <- parameter(1/59.50121)
lambda <- parameter(1e-05)

# Proportion born with maternal immunity
prop <- parameter(1)

# Comparison function - this is used for the likelihood function when fitting to the data. 
# We're telling odin that N and n_infection can be found in the data and that n_infection
# should follow a binomial distribution with parameters N (from the data) and R
# from our model.
N <- data()
n_infection <- data()
n_infection ~ Binomial(N, R)