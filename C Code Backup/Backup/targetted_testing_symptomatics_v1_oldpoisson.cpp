#include <iostream>
#include <stdio.h>
#include <random>
#include <algorithm>
#include <math.h>
#include <time.h>

/* Random Number Generators **********************************/

int poisson(double mu)  // Generate a Poisson Random Number with mean mu
{
    std::random_device rd;
    std::mt19937 gen(rd());

    std::poisson_distribution<> d(mu);

    return d(gen);
}

double uniform()        // Generate a uniform random number between [0,1)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<> d(0.0,1.0);

    return d(gen);
}

double uniform(double a,double b)  // Generate a uniform random number in [a, b)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<> d(a,b);

    return d(gen);
}

int randint(int a,int b) // Generate a uniform random integer in [a, b)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_int_distribution<> d(a,b-1);

    return d(gen);
}


int randint(int b) // Generate a uniform random number in [0, b)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_int_distribution<> d(0,b-1);

    return d(gen);
}

/*****************************************************************/


/* Global variables *******************************************/

const int n_pop = 10000;    // Total population
const int n_loc = 2750;     // Total number of locations
const int n_net = 250;      // Total number of networks
const int n_overlap = 251;  // Start position of home locations + 1
const int n_hospitals = 10; // Set the number of hospitals.

int n_asym = 10;            // Initial number of asymptomatics

char states[][3] = {"S","A","P", "MI","SI","R","H"};
                  // 0   1   2    3    4    5   6

const int n_states = sizeof(states)/sizeof(states[0]); // Number of states
const int n_events = 8;     // Different events in the model
                            // (currently) S->A, A->I, A->R, I->R, I->H, H->R
const int person_attr = 4;  // Number of attributes in the pop array

double Cpars[] = {1,  1,   1,    1,   0.1,  0.1}; // Contact parameters for transmission
               // SA  SP  S-MI  S-SI  S-H   S-Q

int total_loc_confined_time = 10;

int days_bw_hq_tests = 1;
int days_bw_lq_tests = 1;

double rate[n_states][n_states];

int pop[n_pop][person_attr] = { };
int link_weight[n_pop] = { };
int n_per_location[n_loc][n_states] = { };
int people_linked_to[n_loc][n_pop] = {};

bool quarantine_confined = true;

/*************************************************************/

/* Other functions ***********************************************/

int sum(int array[],int len){

  int sum = 0;

  for(int i =0; i<len;i++){
    sum += array[i];
  }

  return sum;
}

int sum(double array[],int len){

  int sum = 0;

  for(int i =0; i<len;i++){
    sum += array[i];
  }

  return sum;
}

int sum(int array[][n_states], const int len, int axis){

  int sum = 0;

  for(int i =0; i<len;i++){
    sum += array[i][axis];
  }

  return sum;
}


/*****************************************************************/




void create_person(int array[], int state, int home){

  int net = randint(n_net);
  array[0] = state;
  array[1] = home;
  array[2] = net;
  array[3] = home;

}

void createPopulation(){

  for (int i=0; i<n_pop;i++){
    int home = randint(n_overlap-1, n_loc);  // Assign random homes per person
    create_person(pop[i], 0, home);

    n_per_location[pop[i][3]][0] += 1;
    link_weight[i] = 2;

  }


    // Creat list of people to be set as asymptomatics

    int temp_people[n_asym] = {};

    int random_person = 0;

    for (int i = 0; i<n_asym;i++){
      bool flag = false;
      while(flag==false){
        random_person = randint(n_pop);
        flag = true;

        for(int j=0; j<n_asym; j++){
          if(temp_people[j] == random_person){
            flag = false; break;
          }

        }
      }
      temp_people[i] = random_person;
    }

    for(int i=0; i<n_asym;i++){

      pop[ temp_people[i] ][0] = 1;

      n_per_location[ pop[temp_people[i]][3] ][0] -= 1; // Reduce number of susceptibles
                                                        // in temp_person's current locations
      n_per_location[ pop[temp_people[i]][3] ][1] += 1; // Increase number of susceptibles
                                                        // in temp_person's current locations

    }

    for(int i=0;i<n_states;i++){
    printf("%s %s %i\n", "Initial number of",states[i],sum(n_per_location, n_loc, i));
  }

    for(int i=0; i<n_loc;i++){

      int counter = 0;                              // Counter to count people in location

      for(int j=0; j<n_pop;j++){

          if(pop[j][1] == i || pop[j][2] == i){     // If the person's home or work location is i
              people_linked_to[i][counter] = j;     // Add them to this location in people_linked_to
              counter++;
          }
        people_linked_to[i][counter] = -1;          // Add -1 at the end of each list
      }
    }

}

void Targetted_Run(double dt, double rate[][n_states], double Tpars[][4], int tf, bool lock_homes, double begin_at, double test_frac){

  bool is_confined[n_pop] = {};
  bool being_tested[n_pop]= {};
  bool loc_confined[n_loc]= {};

  int test_result[n_pop] = {};
  int next_test_date[n_pop]={};

  int result_declared_date[n_pop]={};
  // for(int i =0; i<n_pop;i++){result_declared_date[i] = -1000;}
  std::fill_n(result_declared_date,n_pop,-1000);

  int loc_confined_time[n_loc]= {};
  // for(int i =0; i<n_loc;i++){loc_confined_time[i]    = -1000;}
  std::fill_n(loc_confined_time,n_pop,-1000);


  int hq_tests_conducted = 0;
  int lq_tests_conducted = 0;
  int tests_conducted    = 0;
  int results_declared   = 0;
  int locations_moved    = 0;

  double v[n_states][n_events] = {{-1, -1, 0, 0, 0, 0, 0, 0},   // Array to denote change per state per event
                                  {+1,  0,-1, 0, 0, 0, 0, 0},   // The row indexes the state variables, and the
                                  { 0, +1, 0,-1,-1, 0, 0, 0},   // column indexes the events. The number represents
                                  { 0,  0, 0,+1, 0,-1, 0, 0},   // the change in that particular state, given an event.
                                  { 0,  0, 0, 0,+1, 0,-1, 0},
                                  { 0,  0,+1, 0, 0,+1, 0,+1},
                                  { 0,  0, 0, 0, 0, 0,+1,-1}};

  double transitions[n_states][n_states] = {};
  double r[n_events] = {};   // Array to store rates per event
  int K[n_events] = {};     // Array to store numbers of event that occur in dt (Poisson distributed)

  double alpha = 1 - Cpars[5];

  int tests_available_daily = test_frac/100 * n_pop;
  int tests_today           = test_frac/100 * n_pop;

  int lq_tests_daily = Tpars[0][3]*tests_available_daily;
  int lq_tests_today = Tpars[0][3]*tests_available_daily;

  int hq_tests_daily = Tpars[1][3]*tests_available_daily;
  int hq_tests_today = Tpars[1][3]*tests_available_daily;

  // printf("%s %i\n", "LQ tests today ", lq_tests_today);
  // printf("%s %i\n", "HQ tests today ", hq_tests_today);
  // printf("%s %i\n", "All tests today", tests_today);

  double lq_sens = Tpars[0][0];
  double lq_spec = Tpars[0][1];
  double lq_delay= Tpars[0][2];

  double hq_sens = Tpars[1][0];
  double hq_spec = Tpars[1][1];
  double hq_delay= Tpars[1][2];

  double t = 0.0;
  int day  = 0;

  bool start_testing = false;

  while(t<tf){

    for(int i=0; i< n_loc; i++){
      int j=0,k=0,counter=0,p;
      bool flag = 0;
      // Find people in this location.

      int N = 0;

      for(j=0; j<n_states; j++){ N += n_per_location[i][j];} // Find total number in  location

      int ind[N];  // Indices of people in this location

      counter = 0;

      for(j=0; j<n_pop; j++){       // Loop over all people_linked_to (not pop!)

        if(people_linked_to[i][j]==-1){break;}      // If you reach the "end" of people_linked_to, break
        else if(pop[people_linked_to[i][j]][3]==i){ // Otherwise, if the person is currently here
          ind[counter] = people_linked_to[i][j];
          counter++;
        }
      }

      /* Test to check for counter==N **/
      if(N!=counter){printf("%s %i %s %i %s %i\n","Error! Number of people don't match in location",i,"N=",N,"counter=",counter );}

      if(N==0){continue;}

      double V = N - alpha*n_per_location[i][4];   // Spatial damping parameter

      int conf_by_state_in_loc[n_states]={};

      if(quarantine_confined == true){

        for(j = 0; j<N;j++){      // Go over people in location
            p = ind[j];           // index of person
            int* person = pop[p];     // slice of pop array
            if(is_confined[p] == true){conf_by_state_in_loc[person[0]]++; } // If the person is confined, increment the conf_by_state_in_loc of their state
        }
      }

      // Get the different rates (to be Poisson distributed later)

      r[0] = rate[0][1] * n_per_location[i][0] * 1/V *
                                                (Cpars[0]*(n_per_location[i][1]- conf_by_state_in_loc[1]*alpha) +   // SA
                                                 Cpars[1]*(n_per_location[i][2]- conf_by_state_in_loc[2]*alpha) +   // SP
                                                 Cpars[2]*(n_per_location[i][3]- conf_by_state_in_loc[3]*alpha) +   // S-MI
                                                 Cpars[3]*(n_per_location[i][4]- conf_by_state_in_loc[4]*alpha) +   // S-SI
                                                 Cpars[4]*n_per_location[i][6]);                                    // SH
     r[1] = rate[0][2] * n_per_location[i][0] * 1/V *
                                               (Cpars[0]*(n_per_location[i][1]- conf_by_state_in_loc[1]*alpha) +   // SA
                                                Cpars[1]*(n_per_location[i][2]- conf_by_state_in_loc[2]*alpha) +   // SP
                                                Cpars[2]*(n_per_location[i][3]- conf_by_state_in_loc[3]*alpha) +   // S-MI
                                                Cpars[3]*(n_per_location[i][4]- conf_by_state_in_loc[4]*alpha) +   // S-SI
                                                Cpars[4]*n_per_location[i][6]);                                    // SH
      r[2] = rate[1][5] * n_per_location[i][1];
      r[3] = rate[2][3] * n_per_location[i][2];
      r[4] = rate[2][4] * n_per_location[i][2];
      r[5] = rate[3][5] * n_per_location[i][3];
      r[6] = rate[4][6] * n_per_location[i][4];
      r[7] = rate[6][5] * n_per_location[i][6];

      // (Whew!)

      flag = false;
      counter = 0;  // counter to check number of runs so that poisson means aren't negative.

      while(flag==false){

        for(j=0; j<n_events; j++){ K[j] = poisson(r[j]*dt);} // Get number of people per event (poisson distributed)
        // for(int j=0; j<n_events; j++){K[j] = uniform(0,r[j]*dt) ;} // Get number of people per event (poisson distributed)


        transitions[0][1] = K[0];
        transitions[0][2] = K[1];
        transitions[1][5] = K[2];
        transitions[2][3] = K[3];
        transitions[2][4] = K[4];
        transitions[3][5] = K[5];
        transitions[4][6] = K[6];
        transitions[6][5] = K[7];

        int n_removed[n_states] ={};     // Array to store number of people ideally removed

        for(j=0; j<n_states;j++){
          for(k=0; k<n_states;k++){n_removed[j] += transitions[j][k]; }
        }

      for(j=0; j<n_states; j++){
        if((n_per_location[i][j] < n_removed[j]) && (n_removed[j]>0)){

          for(k=0; k<n_states;k++){
            transitions[j][k] = round(transitions[j][k]/n_removed[j] * n_per_location[i][j]);
          }
        }
      }

       K[0] = transitions[0][1];
       K[1] = transitions[0][2];
       K[2] = transitions[1][5];
       K[3] = transitions[2][3];
       K[4] = transitions[2][4];
       K[5] = transitions[3][5];
       K[6] = transitions[4][6];
       K[7] = transitions[6][5];

       // backup n_per_location[i]
       int backup[n_states] = {};

       for(j=0; j<n_states;j++){ backup[j] = n_per_location[i][j];}

       // Change n_per_location now
       for(j=0; j<n_states;j++){
         for(k=0; k<n_events;k++){
           n_per_location[i][j] += v[j][k]*K[k];  // n_per_location incremented by Sum_k v_{jk}*K^k
         }
       }

       if(n_per_location[i][0]<0 || n_per_location[i][1]<0 || n_per_location[i][2]<0 || n_per_location[i][3]<0 ||
          n_per_location[i][4]<0 || n_per_location[i][5]<0 || n_per_location[i][6]<0){ // If any n_per_location<0

            printf("%s\n", "Negative populations, shouldn't be happening");
            printf("%s\n", "Rolling back");
            for(j=0; j<n_states;j++){n_per_location[i][j]=backup[j];} // Roll back to backup
            printf("%s\n","Trying again...");
            counter++;

            if(counter == 9){printf("%s\n", "Fatal error: 10 retries didn't get rid of negative populations!");exit(1);}
          }
        else{flag=true;}


        // Add code here to check if all values in transitions are 0, then continue?




      }

    }





    // printf("%f\n",t);
    t+=dt;
    // End While Loop!
  }

}


int main() {

  /* Initial populations ***************************/

  int n0[n_states] = {n_pop-n_asym, n_asym, 0, 0, 0};
             // [S,            A,      I, R, H ]

  /*************************************************/


  printf("%s\n", "BufferLine");
  printf("%s\n", "BufferLine");
  printf("%s\n", "BufferLine");
  printf("%s\n", "BufferLine");
  printf("%s\n", "BufferLine");
  printf("%s\n", "BufferLine");


  // S, A, P, MI, SI, R, H
  // 0  1  2  3   4   5  6

  rate[0][1] = 0.45;  // S -> A
  rate[0][2] = 0.45;  // S -> I
  rate[1][5] = 0.136; // A -> R
  rate[2][3] = 0.205; // P -> MI
  rate[2][4] = 0.205; // P -> SI
  rate[3][5] = 0.1;   // MI -> R
  rate[4][6] = 0.25;  // SI -> H
  rate[6][5] = 0.1;   // H -> R


  double Tpars[2][4] = {{0.75, 0.98, 0, 0},
                        {1.0,  1.0,  0, 1}};


  int tf = 100; // days_bw_hq_tests

  bool lock_homes = true;
  double begin_at = 10;
  double test_frac= 0.5;

  double dt = 0.1;

  createPopulation();

  clock_t start, end; // Measuring how long the program takes to run
  double cpu_time_used;
  start = clock();
  printf("Started\n");
  Targetted_Run(dt, rate, Tpars, tf, lock_homes, begin_at, test_frac );
  printf("Done\n");
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("Total time spent Tau Leaping %lf\n",cpu_time_used);

}
