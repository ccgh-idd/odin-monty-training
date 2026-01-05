# Initial states:
initial(S) <- N - I0
initial(I) <- I0
initial(R) <- 0
initial(incidence, zero_every = 1) <- 0 # incidence is reset to zero at each time step

# Core equations for transitions between compartments:
update(S) <- S - n_SI
update(I) <- I + n_SI - n_IR
update(R) <- R + n_IR
update(incidence) <- incidence + n_SI

# Individual probabilities of transition:
n_SI <- Binomial(S, p_SI) # S to I
n_IR <- Binomial(I, p_IR) # I to R

# Transition parameters
p_SI <- 1 - exp(-beta * I / N * dt)
p_IR <- 1 - exp(-gamma * dt)

# User defined parameters - you can add default values in parentheses and they 
# can be changed when initialising the model:
beta <- parameter()
gamma <- parameter()
I0 <- parameter()
N <- 1000

# Comparison to data - this is used for the likelihood function when fitting to the data. 
# We're telling odin that the cases can be found in the data and that they
# should follow a Poisson distribution. They are compared with the incidence output of our model.
cases <- data()
cases ~ Poisson(incidence)