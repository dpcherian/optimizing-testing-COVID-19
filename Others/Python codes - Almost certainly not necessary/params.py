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
n_hospitals = 2  # Set the number of hospitals.

n_asym = 10      # Initial number of asymptomatics

# Initial populations ###########

n0 = np.array([n_pop-n_asym, n_asym, 0, 0, 0, 0, 0],dtype=np.int64)
           # [S,            A,      P,  MI,SI R, H ]

##################################


labels      = np.array([0,1,2,3,4,5,6],dtype=np.int64) # Number of states a person can be in, [S,A,P,MI,SI,R,H] (Not very important, just for bookkeeping)
person_attr = 4           # 4 attributes for a person: state, home, network, current location

### Rates #############

n_states = len(labels)
n_events = 9              # Different events in the model (currently) S->A, S->P, A->R, P-> MI, P-> SI, MI-> R, SI->R, SI->H, H->R
    
rate_array = np.zeros((n_states,n_states))

k_S   = 0.25;
gamma = 0.5; 

k_A   = 0.142
k_P   = 0.5; 
delta = 0.85;

k_MI  = 0.1; 
k_SI  = 0.5; 
sigma = 0.8; 

k_H   = 0.1;   


rate_array[0][1] = gamma*k_S;      # S -> A
rate_array[0][2] = (1-gamma)*k_S;  # S -> I
rate_array[1][5] = k_A;            # A -> R
rate_array[2][3] = delta*k_P;      # P -> MI
rate_array[2][4] = (1-delta)*k_P;  # P -> SI
rate_array[3][5] = k_MI;           # MI -> R
rate_array[4][5] = (1-sigma)*k_SI; # SI -> R
rate_array[4][6] = sigma*k_SI;     # SI -> H
rate_array[6][5] = k_H;            # H -> R


Cpars = np.array([1,  1,  1,  1, 0.1,  0.1],dtype=np.float64)  # Transmissivity of different individuals. Asymptomatic and symptomatics 
                                     # are just as infectious, hospitalised and quarantined are 10 times less infectious
       # Contact parameters


total_loc_confined_time = 14

days_bw_hq_tests = 1
days_bw_lq_tests = 1

dt = 0.1
