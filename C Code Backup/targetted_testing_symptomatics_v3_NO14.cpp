#include <iostream>
#include <stdio.h>
#include <random>
#include <algorithm>
#include <math.h>
#include <time.h>

/* Random Number Generators **********************************/

std::random_device rd;
std::mt19937 gen(rd());

int poisson(double mu)  // Generate a Poisson Random Number with mean mu
{
    std::poisson_distribution<> d(mu);

    return d(gen);
}

double uniform()        // Generate a uniform random number between [0,1)
{
    std::uniform_real_distribution<> d(0.0,1.0);

    return d(gen);
}

double uniform(double a,double b)  // Generate a uniform random number in [a, b)
{
    std::uniform_real_distribution<> d(a,b);

    return d(gen);
}

int randint(int a,int b) // Generate a uniform random integer in [a, b)
{
    std::uniform_int_distribution<> d(a,b-1);

    return d(gen);
}


int randint(int b) // Generate a uniform random number in [0, b)
{
    std::uniform_int_distribution<> d(0,b-1);

    return d(gen);
}

/*****************************************************************/


/* Global variables *******************************************/

const int n_pop = 10000;    // Total population
const int n_loc = 2750;     // Total number of locations
const int n_net = 250;      // Total number of networks
const int n_overlap = 251;  // Start position of home locations + 1
const int n_hospitals = 2;  // Set the number of hospitals.

int n_asym = 10;            // Initial number of asymptomatics

char states[][3] = {"S","A","P", "MI","SI","R","H"};
                  // 0   1   2    3    4    5   6

int S = 0;
int A = 1;
int P = 2;
int MI= 3;
int SI= 4;
int R = 5;
int H = 6;

const int n_states = sizeof(states)/sizeof(states[0]); // Number of states
const int n_events = 9;     // Different events in the model
                            // (currently) S->A,  S->P,  A->R, P->MI, P->SI, MI->R, SI->R, SI->H,  H->R
const int person_attr = 4;  // Number of attributes in the pop array

double Cpars[] = {1,  1,   1,    1,   0.1,  0.1}; // Contact parameters for transmission
               // SA  SP  S-MI  S-SI  S-H   S-Q

/* Initial populations ***************************/

int n[n_states];// {n_pop-n_asym, n_asym, 0, 0,  0,  0, 0 };
                     // [S,        A,    P, MI, SI, R, H ]

/*************************************************/


int total_loc_confined_time = 14;

int days_bw_hq_tests = 0;
int days_bw_lq_tests = 0;

double rate[n_states][n_states];

int pop[n_pop][person_attr] = { };
int link_weight[n_pop] = { };
int n_per_location[n_loc][n_states];
int people_linked_to[n_loc][n_pop] = {};

bool quarantine_confined = true;               // Change to false if you don't want to confine individuals who've tested positive


const int op_width = 1 + n_states + 3;         // Width of output array: <time>(1) <Number of states> <Testing details>(3)
int positives = 0;
/*****************************************************************/

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

  array[0] = state;            // State of individual
  array[1] = home;             // Home location
  array[2] = randint(n_net);   // Random work location from [0,n_net)
  array[3] = home;             // Current location (initialised as home)

}

void createPopulation(){

  n[S] = n_pop-n_asym;
  n[A] = n_asym;
  n[P] = 0;
  n[MI]= 0;
  n[SI]= 0;
  n[R] = 0;
  n[H] = 0;

  for(int i=0;i<n_loc;i++){
    for(int j=0;j<n_states;j++){n_per_location[i][j]=0;}
  }

  for (int i=0; i<n_pop;i++){
    int home = randint(n_overlap-1, n_loc);  // Assign random homes per person from [n_overlap-1,n_loc)
    create_person(pop[i], S, home);          // Assign everyone to S=0

    n_per_location[pop[i][3]][S] += 1;       // Increase number of S in current location
    link_weight[i] = 2;                      // Number of movements per person per day

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

      n_per_location[ pop[temp_people[i]][3] ][S] -= 1; // Reduce number of susceptibles
                                                        // in temp_person's current locations
      n_per_location[ pop[temp_people[i]][3] ][A] += 1; // Increase number of asymptomatics
                                                        // in temp_person's current locations

    }

    for(int i=0;i<n_states;i++){
    printf("%s %s %i\n", "Initial number of",states[i],sum(n_per_location, n_loc, i));
  }

    // Create a people_linked_to list of all people associated with a location (either as work or home).

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


void writetofile(int output[][op_width], int tf, double Tpars[2][4], double begin_at, double test_frac, double time_taken, int details[9]){
  // FOR REFERENCE: int details[] = {lock_homes, quarantine_confined, lq_tests_conducted, hq_tests_conducted, tests_conducted, results_declared, locations_moved, hcw_recovered};
  // Write output to a FILE

  time_t rawtime;
  struct tm * timeinfo;
  char dateString[80];
  time (&rawtime);
  timeinfo = localtime(&rawtime);

  strftime(dateString,sizeof(dateString),"%y%m%d%H%M%s",timeinfo);

  char str[200];

  sprintf(str,"Targetted-Testing_BeginAt_%.2f_DTR_%.2f_RAT_%d_%d_PCR_%d_%d_%s.txt",begin_at,test_frac,(int)Tpars[0][0],(int)Tpars[0][1],(int)Tpars[1][0],(int)Tpars[1][1],dateString);
  FILE *fpt=(FILE *)fopen(str,"wt"); // Open the file to print output

  fprintf(fpt,"###### TEST LOG ####################\n");
  fprintf(fpt,"# Time taken               : %.2f s\n",time_taken);

  fprintf(fpt,"# Test Parameters: \n");
  fprintf(fpt,"# %.2f %.2f %.2f %.2f\n", Tpars[0][0], Tpars[0][1], Tpars[0][2], Tpars[0][3]);
  fprintf(fpt,"# %.2f %.2f %.2f %.2f\n", Tpars[1][0], Tpars[1][1], Tpars[1][2], Tpars[1][3]);

  fprintf(fpt,"# Homes Quarantined?       : %s\n", details[0] ? "True": "False");
  fprintf(fpt,"# Confined Less Infective? : %s\n", details[1] ? "True": "False");
  fprintf(fpt,"# Testing Started When     : %.2f%% recovered\n", begin_at);
  fprintf(fpt,"# Fraction Tested Daily    : %.2f%%\n", test_frac);

  fprintf(fpt,"# LQ Tests Done in total   : %d\n", details[2]);
  fprintf(fpt,"# HQ Tests Done in total   : %d\n", details[3]);
  fprintf(fpt,"# All Tests Done in total  : %d\n", details[4]);
  fprintf(fpt,"# Results Given in total   : %d\n", details[5]);
  fprintf(fpt,"# Locations Moved in total : %d\n", details[6]);
  fprintf(fpt,"# Total recovered HCW      : %d\n", details[7]);
  fprintf(fpt,"# Total HCW                : %d\n", details[8]);

  fprintf(fpt,"# Rate Array: \n");
  for(int i=0; i<n_states;i++){fprintf(fpt,"# ");for(int j=0; j<n_states;j++){fprintf(fpt,"%5g ", rate[i][j]);}fprintf(fpt,"\n");}
  fprintf(fpt,"###### END LOG #####################\n");
  fprintf(fpt,"#\n");



  for(int i=0;i<tf+1;i++){
    for(int j=0; j<op_width;j++){
      fprintf(fpt, "%i ", output[i][j]);
    }
    fprintf(fpt, "\n");
  }
  fflush(fpt);  // Flush the buffer
}

void shuffle(int start, int end, int *array)
{
    if (end-start > 1)
    {
        for (int i = start; i < end; i++)
        {
          int j = i + rand() / (RAND_MAX / (end - i) + 1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}



void Targetted_Run(double dt, double rate[][n_states], double Tpars[][4], int tf, bool lock_homes, bool quarantine_when_sample_taken, double begin_at, double test_frac){

  clock_t start, end; // Measuring how long the function takes to run
  double cpu_time_used;
  start = clock();

  bool is_confined[n_pop] = {};
  bool being_tested[n_pop]= {};
  bool loc_confined[n_loc]= {};

  int test_result[n_pop] = {};
  int next_test_date[n_pop]={};

  int result_declared_date[n_pop]={};
  std::fill_n(result_declared_date,n_pop,-1000);  // Set all values of result_declared_date to -1000

  int loc_confined_time[n_loc]= {};
  std::fill_n(loc_confined_time,n_loc,-1000);     // Set all values of loc_confined_time to -1000

  int person_isolated_time[n_pop] = {};
  std::fill_n(person_isolated_time, n_pop,-1000);


  int hq_tests_conducted = 0;
  int lq_tests_conducted = 0;
  int tests_conducted    = 0;
  int results_declared   = 0;
  int locations_moved    = 0;
                                // SA   SP  AR  PMI PSI MIR SIR SIH  HR
  double v[n_states][n_events] = {{-1,  -1,  0,  0,  0, 0,  0,  0,  0}, // S  // Array to denote change per state per event
                                  {+1,   0, -1,  0,  0, 0,  0,  0,  0}, // A  // The row indexes the state variables, and the
                                  { 0,  +1,  0, -1, -1, 0,  0,  0,  0}, // P  // column indexes the events. The number represents
                                  { 0,   0,  0, +1,  0,-1,  0,  0,  0}, // MI // the change in that particular state, given an event.
                                  { 0,   0,  0,  0, +1, 0, -1, -1,  0}, // SI
                                  { 0,   0, +1,  0,  0,+1, +1,  0, +1}, // R
                                  { 0,   0,  0,  0,  0, 0,  0, +1, -1}};// H

  int transitions[n_states][n_states] = {};               // Array to store transitions
  double r[n_events] = {};                                // Array to store rates per event
  int K[n_events] = {};                                   // Array to store numbers of event that occur in dt (Poisson distributed)

  double alphaH = 1 - Cpars[4];                           // Quantity by which V is reduced per hospitalised individual
  double alphaQ = 1 - Cpars[5];

  int tests_available_daily = test_frac/100 * n_pop;
  int tests_remaining_today = test_frac/100 * n_pop;

  int lq_tests_daily = Tpars[0][3]*tests_available_daily;
  int lq_tests_today = Tpars[0][3]*tests_available_daily;

  int hq_tests_daily = Tpars[1][3]*tests_available_daily;
  int hq_tests_today = Tpars[1][3]*tests_available_daily;

  double lq_sens = Tpars[0][0];
  double lq_spec = Tpars[0][1];
  double lq_delay= Tpars[0][2];

  double hq_sens = Tpars[1][0];
  double hq_spec = Tpars[1][1];
  double hq_delay= Tpars[1][2];

  double t = 0.0;
  int day  = 0;

  int output[tf+1][op_width];   // Output array, to be printed to file

  // First line of output
  output[day][0] = day; for(int s=0;s<n_states;s++){output[day][s+1]=n[s];} output[day][n_states+1]=hq_tests_conducted;output[day][n_states+2]=lq_tests_conducted;output[day][n_states+3]=tests_remaining_today;

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

      int conf_by_state_in_loc[n_states]={};


      if(quarantine_confined == true){
        for(j = 0; j<N;j++){      // Go over people in location
            p = ind[j];           // index of person
            // int* person = pop[p]; // slice of pop array
            if(is_confined[p] == true){conf_by_state_in_loc[pop[p][0]]++; } // If the person is confined, increment the conf_by_state_in_loc of their state
        }
        // printf("Location %i\n", i);for(int qq=0;qq<n_states;qq++){printf("%i ",conf_by_state_in_loc[qq]);}printf("\n");
      }

      double V = N - alphaH*n_per_location[i][H];   // Spatial damping parameter (adjusted by alphaH for hospitals, and alphaQ for homes where there are people confined)


      // Get the different rates (to be Poisson distributed later)

      // FOR REFERENCE: Cpars = {0 = SA  1 = SP  2 = S-MI  3 = S-SI  4 = S-H   5 = S-Q}

      // S->A
      r[0] = rate[S][A] * n_per_location[i][S] * 1/V *
                                                (Cpars[0]*(n_per_location[i][A]- conf_by_state_in_loc[A]*alphaQ) +    // SA
                                                 Cpars[1]*(n_per_location[i][P]- conf_by_state_in_loc[P]*alphaQ) +    // SP
                                                 Cpars[2]*(n_per_location[i][MI]-conf_by_state_in_loc[MI]*alphaQ) +   // S-MI
                                                 Cpars[3]*(n_per_location[i][SI]-conf_by_state_in_loc[SI]*alphaQ) +   // S-SI
                                                 Cpars[4]*n_per_location[i][H]);                                      // SH
     // S->P
     r[1] = rate[S][P] * n_per_location[i][S] * 1/V *
                                               (Cpars[0]*(n_per_location[i][A]- conf_by_state_in_loc[A]*alphaQ) +    // SA
                                                Cpars[1]*(n_per_location[i][P]- conf_by_state_in_loc[P]*alphaQ) +    // SP
                                                Cpars[2]*(n_per_location[i][MI]-conf_by_state_in_loc[MI]*alphaQ) +   // S-MI
                                                Cpars[3]*(n_per_location[i][SI]-conf_by_state_in_loc[SI]*alphaQ) +   // S-SI
                                                Cpars[4]*n_per_location[i][H]);                                      // SH
      // A->R
      r[2] = rate[A][R] * n_per_location[i][A];
      // P->MI
      r[3] = rate[P][MI]* n_per_location[i][P];
      // P->SI
      r[4] = rate[P][SI]* n_per_location[i][P];
      // MI->R
      r[5] = rate[MI][R]* n_per_location[i][MI];
      // SI->R
      r[6] = rate[SI][R] * n_per_location[i][SI];
      // SI->H
      r[7] = rate[SI][H] * n_per_location[i][SI];
      // H->R
      r[8] = rate[H][R] * n_per_location[i][H];

      // (Whew!)

      flag = false;
      counter = 0;  // counter to check number of runs so that poisson means aren't negative.

      while(flag==false){

        for(j=0; j<n_events; j++){K[j] = poisson(r[j]*dt);} // Get number of people per event (poisson distributed)

        transitions[S][A] = K[0];
        transitions[S][P] = K[1];
        transitions[A][R] = K[2];
        transitions[P][MI]= K[3];
        transitions[P][SI]= K[4];
        transitions[MI][R]= K[5];
        transitions[SI][R]= K[6];
        transitions[SI][H]= K[7];
        transitions[H][R] = K[8];

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

       K[0] = transitions[S][A];
       K[1] = transitions[S][P];
       K[2] = transitions[A][R];
       K[3] = transitions[P][MI];
       K[4] = transitions[P][SI];
       K[5] = transitions[MI][R];
       K[6] = transitions[SI][R];
       K[7] = transitions[SI][H];
       K[8] = transitions[H][R];

       // backup n_per_location[i]
       int backup[n_states] = {};

       for(j=0; j<n_states;j++){ backup[j] = n_per_location[i][j];}

       // Change n_per_location now
       for(j=0; j<n_states;j++){
         for(k=0; k<n_events;k++){
           n_per_location[i][j] += v[j][k]*K[k];  // n_per_location incremented by Sum_k v_{jk}*K^k
         }
       }

       if(n_per_location[i][S] <0 || n_per_location[i][A]<0 || n_per_location[i][P]<0 || n_per_location[i][MI]<0 ||
          n_per_location[i][SI]<0 || n_per_location[i][R]<0 || n_per_location[i][H]<0){ // If any n_per_location<0

            printf("%s\n", "Negative populations, shouldn't be happening");
            printf("%s\n", "Rolling back");
            for(j=0; j<n_states;j++){n_per_location[i][j]=backup[j];} // Roll back to backup
            printf("%s\n","Trying again...");
            counter++;

            if(counter == 9){printf("%s\n", "Fatal error: 10 retries didn't get rid of negative populations!");exit(1);}
          }
        else{flag=true;}

        // for(int j=0;j<n_states;j++){for(int k = 0; k<n_states;k++){if(transitions[j][k]!=0){printf("%i\n",transitions[j][k] );}}}

        // Add code here to check if all values in transitions are 0, then continue?


        bool done[n_pop]={};

        for(int j=0;j<n_states;j++){
            for(int k=0;k<n_states;k++){
                int l = 0;
                for(int m=0; m<N;m++){
                  if(l==transitions[j][k]){break;}

                  if(pop[ind[m]][0]==j && done[m]==false){
                    pop[ind[m]][0] = k;
                    n[j]--; n[k]++;

                    if(k==R && j==H){ // If they recover
                      // is_confined[ind[m]] = false; // Remove confinement

                       // Slight difference from python, send ALL recovered home
                      n_per_location[ pop[ind[m]][3] ][k]--; // Decrement number of k (=R) in current locations
                      pop[ind[m]][3] = pop[ind[m]][1];       // Send them home
                      n_per_location[ pop[ind[m]][3] ][k]++; // Increment number of k (=R) in current locations

                      ////

                    }
                    else if(k==H){ // If they're hospitalised
                      int h = randint(n_hospitals);                       // Choose a hospital at random
                      while(h==pop[ind[m]][2]){h = randint(n_hospitals);} // Make sure a HCW is not hospitalised in the same hospital

                      for(int s=0;s<n_pop;s++){                           // Link them to this hospital forever. (They are never unlinked, even after they leave)
                        if(people_linked_to[h][s] == -1){people_linked_to[h][s]=ind[m];people_linked_to[h][s+1]=-1;break;}
                      }

                      is_confined[ind[m]] = false;           // Remove confinement

                      test_result[ind[m]] = 0;

                      n_per_location[ pop[ind[m]][3] ][k]--; // Decrement number of k (=H) in current location
                      pop[ind[m]][3] = h;                    // Send them to a random hospitals
                      n_per_location[ pop[ind[m]][3] ][k]++; // Increment number of k (=H) in current location

                    }

                    done[m] = true; l++;                     // Don't deal with this person again, and then increment l
                  }
                }
            }
        }
      }

    }


    if(t>=day){

      /************* TARGETTED TESTING ******************/


      if(n[R]>=begin_at/100 * n_pop){

        // Find all symptomatics (?) Right now just looking at SI and MI

        int list_of_sym[n_pop]={};
        int n_sym = 0;

        int list_of_remaining[n_pop]={};       // Array to store remaining people for random testing
        int n_remaining  = 0;                  // Number of such people

        int not_eligible_for_testing = 0;              // Number of people not eligible for testing today

        for(int i=0; i<n_pop;i++){

          /******************** Find list of targets and remaining people **********************/


          if((pop[i][0]==MI||pop[i][0]==SI) && day>=next_test_date[i] && being_tested[i]==false && is_confined[i]==false){
            list_of_sym[n_sym] = i; n_sym++;
          }
          else if(pop[i][0]!=H && being_tested[i]==false && is_confined[i]==false){//-->corrected? WRONG! Remaining should be those who are not being tested AND NOT CONFINED?
            list_of_remaining[n_remaining] = i; n_remaining++;
          }
          else{ not_eligible_for_testing++;}

          /****************************** Done finding lists **********************************/
        }

        list_of_sym[n_sym] = -1;                 // Set -1 to mark the end of this array
        list_of_remaining[n_remaining] = -1;     // Set -1 to mark end of this array

        if(n_sym>0){
          int targetted_tests_done_today =   std::min(tests_remaining_today,n_sym);  // Tests done this dt is the minimum of the tests available and
                                                                                     // the people to be tested. As the day progresses, the tests
                                                                                     // available drops lower until it's 0, and no testing happens.

          shuffle(0, n_sym, list_of_sym); // Shuffle list of people to be tested.

          // for(int lo=targetted_tests_done_today;lo<n_sym;lo++){ list_of_remaining[n_remaining] = list_of_sym[lo];n_remaining++; } // Add those who couldn't be tested to the n_remaining
                                                                                                                                  // (this is a little pointess: if there are enough tests_conducted
                                                                                                                                  // targetted individuals will always be tested, until no tests
                                                                                                                                  // remain. Meaning if these people couldn't be tested in targetted TESTING
                                                                                                                                  // no random testing is going to happen! But anyway....)

          list_of_remaining[n_remaining] = -1;  // Reset position of -1 to mark new end of this array

          for(int j=0; j<targetted_tests_done_today;j++){

            int si = list_of_sym[j];  // Individual to test

            if(day >= next_test_date[si] && being_tested[si]==false && tests_remaining_today>0){
              // If so, perform a tests
              being_tested[si] = true; tests_conducted++; tests_remaining_today--;

              if(quarantine_when_sample_taken==true){
                //MOVE THEM Home
                n_per_location[pop[si][3]][pop[si][0]]--;
                pop[si][3] = pop[si][1];
                n_per_location[pop[si][3]][pop[si][0]]++;
                is_confined[si] = quarantine_confined;
              }      // Quarantine as soon as sample is taken

              loc_confined[pop[si][1]] = lock_homes;      // Lock home depending on variable `lock_homes`.
              loc_confined_time[pop[si][1]] = t;

              // Targetted testing using HQ tests unless there are none
              int test_type = 1;

              if(hq_tests_today<=0){test_type = 0;} // If no HQ tests available, give them LQ (one of the two is guaranteed, since tests_remaining_today>0)


              if(test_type==0){
                // Do a low quality test

                lq_tests_conducted++; lq_tests_today--;

                next_test_date[si]      = day + days_bw_lq_tests;   // Next day to be considered for a test
                result_declared_date[si]= day + lq_delay;           // Days to wait before result result_declared_date

                if(pop[si][0]>S && pop[si][0]<R){                  // If the person isn't susceptible or recovered or hospitalised
                  if(uniform()<lq_sens){test_result[si]=1;}         // and the test comes back positive, set their test_result
                  else{test_result[si]=-1;}                         // otherwise it's negative
                }
                else if(pop[si][0]!=H){                             // On the other hand, if they are S or R (not testing Hospitalised)
                  if(uniform()>lq_spec){test_result[si]=1;}         // and the test comes back false positive, set their test_result
                  else{test_result[si]=-1;}                         // otherwise it's negative
                }
              }

              else if(test_type==1){
                // Do a high quality test

                hq_tests_conducted++; hq_tests_today--;

                next_test_date[si]      = day + days_bw_hq_tests;   // Next day to be considered for a test
                result_declared_date[si]= day + hq_delay;           // Days to wait before result result_declared_date

                if(pop[si][0]>S && pop[si][0]<R){                  // If the person isn't susceptible or recovered or hospitalised
                  if(uniform()<hq_sens){test_result[si]=1;}         // and the test comes back positive, set their test_result
                  else{test_result[si]=-1;}                         // otherwise it's negative
                }
                else if(pop[si][0]!=H){                             // On the other hand, if they are S or R (not testing Hospitalised)
                  if(uniform()>hq_spec){test_result[si]=1;}         // and the test comes back false positive, set their test_result
                  else{test_result[si]=-1;}                         // otherwise it's negative
                }
              }
            }
          }
        }

        /*************RANDOMLY TESTING POPULATION*************/

        // Test remaining people in population (in the array list_of_remaining) randomly

        int remaining_tests_done_today =   std::min(tests_remaining_today,n_remaining);  // Remaining tests done is the minimum of the tests available and
                                                                                         // available drops lower until it's 0, and no testing happens.
                                                                                         // the people to be tested. As the day progresses, the tests


        if(remaining_tests_done_today>0){               // If there are any tests remaining

            shuffle(0,n_remaining,list_of_remaining);  // Shuffle list of remaining people

            // Test the first "remaining_tests_done_today" people

            for(int j=0; j<remaining_tests_done_today;j++){

              int ri = list_of_remaining[j];           // Individual to test

              // ** THE PART BELOW (EXCEPT FOR LINES 722 to 726) IS ESSENTIALLY THE SAME AS FOR TARGETTED TESTING (with si -> ri)

              being_tested[ri] = true; tests_conducted++; tests_remaining_today--;

              if(quarantine_when_sample_taken==true){
                //MOVE THEM Home
                n_per_location[pop[ri][3]][pop[ri][0]]--;
                pop[ri][3] = pop[ri][1];
                n_per_location[pop[ri][3]][pop[ri][0]]++;

                is_confined[ri] = quarantine_confined;}      // Quarantine as soon as sample is taken

              loc_confined[pop[ri][1]] = lock_homes;      // Lock home depending on variable `lock_homes`.
              loc_confined_time[pop[ri][1]] = t;

              int test_type = randint(2);        // Returns either 0 and 1 with equal probability.

              if(test_type==0 && lq_tests_today<=0){test_type=1;}      // If no LQ tests, give them HQ
              else if(test_type==1 && hq_tests_today<=0){test_type=0;} // and vice versa

              if(test_type==0){
                // Do a low quality test

                lq_tests_conducted++; lq_tests_today--;

                next_test_date[ri]      = day + days_bw_lq_tests;   // Next day to be considered for a test
                result_declared_date[ri]= day + lq_delay;           // Days to wait before result result_declared_date

                if(pop[ri][0]>S && pop[ri][0]<R){                  // If the person isn't susceptible or recovered or hospitalised
                  if(uniform()<lq_sens){test_result[ri]=1;}         // and the test comes back positive, set their test_result
                  else{test_result[ri]=-1;}                         // otherwise it's negative
                }
                else if(pop[ri][0]!=H){                             // On the other hand, if they are S or R (not testing Hospitalised)
                  if(uniform()>lq_spec){test_result[ri]=1;}         // and the test comes back false positive, set their test_result
                  else{test_result[ri]=-1;}                         // otherwise it's negative
                }
              }

              else if(test_type==1){
                // Do a high quality test

                hq_tests_conducted++; hq_tests_today--;

                next_test_date[ri]      = day + days_bw_hq_tests;   // Next day to be considered for a test
                result_declared_date[ri]= day + hq_delay;           // Days to wait before result result_declared_date

                if(pop[ri][0]>S && pop[ri][0]<R){                  // If the person isn't susceptible or recovered or hospitalised
                  if(uniform()<hq_sens){test_result[ri]=1;}         // and the test comes back positive, set their test_result
                  else{test_result[ri]=-1;}                         // otherwise it's negative
                }
                else if(pop[ri][0]!=H){                             // On the other hand, if they are S or R (not testing Hospitalised)
                  if(uniform()>hq_spec){test_result[ri]=1;}         // and the test comes back false positive, set their test_result
                  else{test_result[ri]=-1;}                         // otherwise it's negative
                }
              }
              // END OF RANDOM TESTING



            }


        }


        /***********************DONE TESTING!***********************/

        for(int i=0;i<n_pop;i++){

            /*********************DECLARING RESULTS********************/

            if(day >= result_declared_date[i] && being_tested[i]==true){

              // First set them to no longer being tested
              being_tested[i] = false;

              results_declared++;

              if(loc_confined[pop[i][1]] == true){ loc_confined[pop[i][1]] = false; loc_confined_time[pop[i][1]]= -1000;} // Remove location confinement (if present)

              // Declare results_declared

              if(test_result[i]==1 && pop[i][0] != H){  // If the result is positive, and the person hasn't already moved to a hospital
                positives++;
                is_confined[i] = quarantine_confined; // Confine them if quarantine_confined==true

                person_isolated_time[i] = day;        // Added new: time of isolation.

                // Move them home
                n_per_location[pop[i][3]][pop[i][0]]--;              // Decremement number in current location
                if(pop[i][3] == pop[i][2]){ pop[i][3] = pop[i][1];}  // If they're at work, send them home.
                n_per_location[pop[i][3]][pop[i][0]]++;              // Incremement number in current location.

                loc_confined[pop[i][1]] = lock_homes;                      // Reconfine their homes
                loc_confined_time[pop[i][1]] = t;

                next_test_date[i] += 0;//14; // TO make compatible with applet.
              }
              else if(test_result[i] == -1 || pop[i][0] == H){       // If they're negative or if they've turned to hospitalised

                is_confined[i] = false;  // Remove confinement

                person_isolated_time[i] = -1000;        // Reset time of isolation.

              }

              test_result[i] = 0;        // Reset test result to 0.
            }

            /***********************DONE DECLARING!**********************/
        }



      }

      // Increment the day, write the output to an array, and reset the number of tests  //
      day++;
      output[day][0] = day; for(int s=0;s<n_states;s++){output[day][s+1]=n[s];} output[day][n_states+1]=hq_tests_conducted;output[day][n_states+2]=lq_tests_conducted;output[day][n_states+3]=tests_remaining_today;
      tests_remaining_today = tests_available_daily;
      lq_tests_today = lq_tests_daily;
      hq_tests_today = hq_tests_daily;
      // printf("%i %i\n",day, positives);

    }




    for(int i=0;i<n_pop;i++){

      /******** Remove confinement if 14 days have passed ********/

      if(is_confined[i]==true && day >= person_isolated_time[i]+14){
        is_confined[i]=false;
        person_isolated_time[i] = -1000;
      }


      /**********************MOVE PEOPLE ABOUT********************/

      if(is_confined[i] == false && loc_confined[pop[i][3]]==false){

        int locs_changed = poisson(link_weight[i]*dt); // Change this to prob. of finding even/modelled
        // printf("%i\n",locs_changed);
        locations_moved += locs_changed;
        locs_changed = locs_changed%2;

        if(locs_changed==1){
          int home_loc = pop[i][1];
          int work_loc = pop[i][2];

          // If they're at home (and aren't hospitalised), move them to work and vice versa
          if(pop[i][3]==home_loc && pop[i][0]!=H){ pop[i][3] = work_loc; n_per_location[home_loc][pop[i][0]]--; n_per_location[work_loc][pop[i][0]]++; }
          else if(pop[i][3]==work_loc && pop[i][0]!=H){ pop[i][3] = home_loc; n_per_location[work_loc][pop[i][0]]--; n_per_location[home_loc][pop[i][0]]++; }

        }
      }
    }

    // Lock or unlock homes

    for(int i=0; i<n_loc;i++){    // This can probably be added to the first location loop, if needed.

      if(loc_confined[i]==true && t-loc_confined_time[i] > total_loc_confined_time){

        for(int j = 0; j<n_pop; j++){if(pop[j][3]==i && is_confined[j]==true){printf("Problem in location %i\n", i);}}

        loc_confined[i]      = false; // Remove confinement
        loc_confined_time[i] = - 1000;
      }
    }



    t+=dt;
    // End While Loop!
  }

  int hcw_recovered = 0;
  int hcw = 0;
  for(int i=0;i<n_pop;i++){               // Find all HCW, and mark those that have recovered.
    if(pop[i][2]<n_hospitals){
      hcw++;
      if(pop[i][0]==R){hcw_recovered++;}
    }
  }


  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("Total time spent Tau Leaping %lf\n",cpu_time_used);


  int details[] = {lock_homes, quarantine_confined, lq_tests_conducted, hq_tests_conducted, tests_conducted, results_declared, locations_moved, hcw_recovered,hcw};
  writetofile(output, tf, Tpars, begin_at, test_frac,cpu_time_used,details);
}



int main() {

  // S, A, P, MI, SI, R, H
  // 0  1  2  3   4   5  6

  double k_S   = 0.3;
  double gamma = 0.5;  // Fraction going from S->A

  double k_A   = 0.143;
  double k_P   = 0.5;
  double delta = 0.85;  // Fraction going from P->MI

  double k_MI  = 0.1;
  double k_SI  = 0.5;
  double sigma = 0.8;

  double k_H = 0.1;

  rate[S][A] = gamma*k_S;        // S -> A
  rate[S][P] = (1 - gamma)*k_S;  // S -> I
  rate[A][R] = k_A;              // A -> R
  rate[P][MI]= delta*k_P;        // P -> MI
  rate[P][SI]= (1-delta)*k_P;    // P -> SI
  rate[MI][R]= k_MI;             // MI -> R
  rate[SI][R]= (1-sigma)*k_SI;   // SI -> R
  rate[SI][H]= sigma*k_SI;       // SI -> H
  rate[H][R] = k_H;              // H -> R


  double Tpars[2][4] = {{0.75, 0.98, 10, 0},
                        {1.0,  1.0,  10, 1}};


  int tf = 200;

  double test_frac= 100;
  double begin_at = 0;

  bool quarantine_when_sample_taken = false;     // Change to false if you don't want people quarantined before results arrive
  bool lock_homes = false;                       // Change to false if you don't want homes to be locked down


  double dt = 0.1;


  clock_t start, end; // Measuring how long the program takes to run
  double cpu_time_used;
  start = clock();
  printf("Started\n");

  for(int i=0;i<10;i++){
  createPopulation();
  Targetted_Run(dt,
                rate,
                Tpars,
                tf,
                lock_homes,
                quarantine_when_sample_taken,
                begin_at,
                test_frac);
  }
  printf("Done\n");
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("Total time spent Tau Leaping %lf\n",cpu_time_used);




}
