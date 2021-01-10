# COVID-19 on Networks

The provided code is a stochastic agent-based simulator designed to be used to model the spread of  COVID-19 (novel coronavirus, SARS-CoV-2) on a network. It has the possibility of including the effects of testing, with two different types of tests. Furthermore, testing can be started when a certain fraction of the population has recovered from the illness.

The provided code was used to produce all the results in the paper "Optimizing testing for COVID-19 in India":
 
 > Cherian, P., Krishna, S., & Menon, G. I., Optimizing testing for COVID-19 in India (2021). MedRxiv, 2020.12.31.20249106. https://doi.org/10.1101/2020.12.31.20249106

Questions or comments can be directed to dpcherian@gmail.com, or on this project's GitHub page.

## Requirements

The code is written in C++, and can be compiled using the following command: 

`g++ -O3 OptimizingTestingForCOVID19inIndia.cpp`

## Description

Individual agents are assigned a number of attributes, including their infection state, home and work locations, current location, status of confinement, last test result, date last tested, date last isolated, and so on. Each location has a similar set of attributes, such as the number of individuals in every state in that location, the quarantine status of the location, and the date it was last quarantined. We use a standard Monte Carlo algorithm where each individual transits between compartments (infection states) with an exponentially distributed time spent in each compartment. A suitably small time-step is chosen (in our case, `dt = 0.01 days`), and individuals transit out of the compartment `i` they are currently in with a probability of `lambda_i dt`. A second random number is chosen in the case of branching to decide into which compartment the individual transits. As a check on our results, we also experimented with Gillespie methods, including the tau-leap algorithm, but decided finally that the requirements of large population sizes, extensive averaging and the complexity of the dynamics favoured simpler Monte Carlo methods. However, we verified that results obtained with these alternative methods were both qualitatively and quantitatively consistent with our Monte Carlo results.

Our networks consist typically of 10,000 individuals and 2,750 locations, although a number of our simulations were carried out for 100,000 individuals and 27,500 locations to check for finite-size effects. (This can be done by modifying the global `mult` variable to 10, for example. We found that these simulations produced similar qualitative trends, and even quantitatively similar proportions of total infected individuals. In all simulations, we chose the ratio of home to work locations to be 10:1, resulting in the number of individuals per home distributed as a Poisson distribution with mean 4, and the number per work location distributed as a Poisson with mean 40. 

We move individuals between home and work locations synchronously at 12 hour intervals. (Note that the rates specified in the code are for a day i.e, 24 hours. Our update probabilities are scaled accordingly.) Once a target fraction of recovered has been crossed, at the end of every day the population is tested with a fixed number of tests. Then any pending results are declared, the numbers of the different tests available daily are reset, and the entire process repeats.

The code is separated into four basic functions:

### `createPopulation()`

The function populates the global arrays `pop`, `n_per_location` and `people_linked_to`. 

The `pop` array is a two-dimensional array (`n_pop x person_attr`), each row of which represents a single individual, with attributes: `[infection_state, home_location, work_location, current_location]`. The infection state can take values `S=0, A=1, P=2, MI=3, SI=4, R=5, H=6`. The work locations are one of those network locations between `0` and `n_net` (default: 250), while the home location is one of those network locations between `n_overlap` and `n_loc` (default: 250 to 2750). The first `n_hos` locations are assumed to be hospitals (default: 2), and those individuals assigned to this location are assumed to be designated healthcare workers (HCWs). Initially, all individuals are assumed to be at home.

The `n_per_location` array is a 2 dimensional array (`n_loc x n_states`), each row of which represents a single location, with the numbers of individuals in each state in that location: `[n_S, n_A, n_P, n_MI, n_SI, n_R, n_H]`.

The `people_linked_to` array is a 2 dimensional array (`n_loc x max_ppl`) which lists the individuals that are associated with each location (either as a home or a work location). This avoids having to loop over the entire population when looping over locations, as individuals are only either in their home or work locations at any given time. The only exception to this is when the individual is hospitalised, in which case they are moved to a random hospital, and they are added to the associated row of `people_linked_to`. They are never unlinked (in case we would like to trace exactly at which hospital an individual had been treated).

### `Targeted_Run(double Tpars[][4], int tf, bool lock_homes, bool quarantine_when_sample_taken, double begin_at, double test_frac, int iter)`

This function conducts a single run for `tf` days with the population defined in `pop`. The other variables in the array are:

- `Tpars[2][4]`, an array containing the details of the different "low quality" or "lq" or "RAT" tests, and "high quality" or "hq" or "PCR" tests. The elements of the array are:  
``{{RAT Sensitivity, RAT Specificity, RAT Test Delay, RAT Fraction in Mixture},``  
`` {PCR Sensitivity, PCR Specificity, PCR Test Delay, PCR Fraction in Mixture}}``

- The boolean flags determine:  
`lock_homes`: If the homes are quarantined or not.  
`quarantine_when_sample_taken`: If the individuals are quarantined when the sample is taken.  
The global variable `quarantine_confined` is `true` by default, it decides if individuals who are quarantined have their infectivity reduced.

- The two numbers:  
`begin_at`: is the target fraction of the population in **percentage** that has recovered before testing is begun (default: 20).  
`test_frac`: is the daily testing rate as a fraction of the population in **percentage** (default: 0, no testing). A `test_frac = 0.1` would mean that there were 10 tests available daily (for a population size of `n_pop=100 000`).

- `iter` is a number to mark different Monte Carlo runs, to avoid overwriting the output files.

### `main()`

This function defines the exit rates from each compartment, as well as the rate array. It also defines the testing parameters (specificity, sensitivity, delays, all defined in `Tpars`), and the simulation time (`tf`). Currently, the `main()` function runs one single Monte Carlo run, but this can be modified using a simple for loop. A single Monte Carlo run takes roughly 10 seconds to complete. Alternatively, the function `create_heatmap` can be called to create a single heatmap as shown in the paper, for a given number of `mc_runs`. A heatmap with `mc_runs=1` takes roughly 90 minutes to run.

### `writetofile(int output[][op_width], int tf, double Tpars[2][4], double begin_at, double test_frac, double time_taken, int details[9],int iter)`

This function creates an output file in the same working directory as the code. The output contains a log of important information (total number of tests, of movements, number of HCW who contracted the illness, etc.) as well as the total number of individuals in each state for every day. The output filename is 

`Targeted_BeginAt_<begin_at>_DTR_<test_frac>_RAT_<rat_sensitivity>_<rat_fraction_PCR_<pcr_sensitivity>_<pcr_fraction>_<random_number_string>.txt`

and its column headers are: 

`[Day,nS,nA,nP,nMI,nSI,nR,nH,PCR_conducted_today, RAT_conducted_today, Tests_remaining_today, Agents_currently_confined, Quarantines_removed_today,Locations_in_quarantine_today]`
