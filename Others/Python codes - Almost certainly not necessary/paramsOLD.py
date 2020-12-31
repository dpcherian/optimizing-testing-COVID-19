import numpy as np
from numpy.random import randint,uniform,choice

# |  L  |  L  |  L  |  L  |  L  |  L  |  L  |  L  | 
# |  N  |  N  |  N  | *N  | *N  |  N  |  -  |  -  | 
# |  -  |  -  |  -  | *H  | *H  |  H  |  H  |  H  |
# ___Only Networks__  (*OVERLAP) _______Homes_________

n_pop = 10000    # Total population
n_loc = 2750     # Total number of locations
n_net = 250      # Total number of networks
n_overlap = 251  # Include the possibility of overlap between work and home locations
n_hospitals = 10 # Set the number of hospitals.

n_asym = 10      # Initial number of asymptomatics

# Initial populations ###########

n0 = np.array([n_pop-n_asym, n_asym, 0, 0, 0],dtype=np.int64)
           # [S,            A,      I, R, H ]

##################################


labels      = np.array([0,1,2,3,4],dtype=np.int64) # Number of states a person can be in, [S,A,I,R,H] (Not very important, just for bookkeeping)
person_attr = 4           # 4 attributes for a person: state, home, network, current location

### Rates #############

n_states = len(labels)
n_events = 6              # Different events in the model (currently) S->A, A->I, A->R, I->R, I->H, H->R
    
rate_array = np.zeros((n_states,n_states))

rate_array[0][1] = 0.45  # S -> A
rate_array[1][2] = 0.205 # A -> I
rate_array[1][3] = 0.136 # A -> R
rate_array[2][3] = 0.1   # I -> R
rate_array[2][4] = 0.25  # I -> H
rate_array[4][3] = 0.1   # H -> R


Cpars = np.array([0,  1,  1,  0,  0.1,  0.1],dtype=np.float64)  # Transmissivity of different individuals. Asymptomatic and symptomatics 
                                     # are just as infectious, hospitalised and quarantined are 10 times less infectious
       # Contact parameters


total_loc_confined_time = 10

days_bw_hq_tests = 1
days_bw_lq_tests = 1

dt = 0.1
