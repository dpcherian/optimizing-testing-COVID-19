#include <iostream>
#include <stdio.h>
#include <random>
#include <algorithm>
#include <math.h>
#include <time.h>

/* Random Number Generators **********************************/

std::random_device rd;
std::mt19937 gen(rd());

double uniform()        // Generate a uniform random number between [0,1)
{
    std::uniform_real_distribution<> d(0.0,1.0);
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

double old_dt;

const int mult = 2;              // Multiplier to scale up problem (1: n_pop = 1 x 10,000, 2: n_pop = 20,000, etc.)

const int n_pop = mult*10000;    // Total population
const int n_loc = mult*2750;     // Total number of locations
const int n_net = mult*250;      // Total number of networks
const int n_overlap = n_net+ 1;  // Start position of home locations + 1
const int n_hospitals = mult*2;  // Set the number of hospitals.

int n_asym = mult*10;            // Initial number of asymptomatics

const int maxppl = 10000;        // Maximum number of people who can be linked to a location (to avoid 2GB limit on arrays)

char states[][3] = {"S","A","P", "MI","SI","R","H"};
                  // 0   1   2    3    4    5   6

const int S = 0;
const int A = 1;
const int P = 2;
const int MI= 3;
const int SI= 4;
const int R = 5;
const int H = 6;

const int n_states = sizeof(states)/sizeof(states[0]); // Number of states
// const int n_events = 9;     // Different events in the model
                               // (currently) S->A,  S->P,  A->R, P->MI, P->SI, MI->R, SI->R, SI->H,  H->R
const int person_attr = 4;     // Number of attributes in the pop array

double Cpars[] = {1,  1,   1,    1,   0.1,  0.1}; // Contact parameters for transmission
               // SA  SP  S-MI  S-SI  S-H   S-Q

/* Initial populations ***************************/

int n[n_states];// {n_pop-n_asym, n_asym, 0, 0,  0,  0, 0 };
                     // [S,        A,    P, MI, SI, R, H ]

/*************************************************/

int total_loc_confined_time = 14;
int total_isolation_time = 14;
int if_positive_test_after = 14;

int days_bw_hq_tests = 0;
int days_bw_lq_tests = 0;

double rate[n_states][n_states] = {};

int pop[n_pop][person_attr] = { };
//int link_weight[n_pop] = { };
int n_per_location[n_loc][n_states];
int people_linked_to[n_loc][maxppl] = {};

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

int sum(bool array[],int len){

  int sum = 0;

  for(int i =0; i<len;i++){
    sum += array[i];
  }

  return sum;
}


double sum(double array[],int len){

  double sum = 0;

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
  
  // Reset n_per_location
  for(int i=0;i<n_loc;i++){
    for(int j=0;j<n_states;j++){n_per_location[i][j]=0;}
  }

  for (int i=0; i<n_pop;i++){
    int home = randint(n_overlap-1, n_loc);  // Assign random homes per person from [n_overlap-1,n_loc)
    create_person(pop[i], S, home);          // Assign everyone to S=0

    n_per_location[pop[i][3]][S] += 1;       // Increase number of S in current location
//    link_weight[i] = 2;                      // Number of movements per person per day

  }


    // Create list of people to be set as asymptomatics
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

    // Create a people_linked_to list of all people associated with a location (either as work or home).
   
    // Start by resetting old list
    for(int i=0;i<n_loc;i++){
      for(int j=0;j<maxppl;j++){
        if(people_linked_to[i][j]==-1){break;}
        people_linked_to[i][j]=-1;
      }
    }

    // Create new list
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


void writetofile(int output[][op_width], int tf, double Tpars[2][4], double begin_at, double test_frac, double time_taken, int details[9],int iter){
  // FOR REFERENCE: int details[] = {quarantine_confined, lock_homes, quarantine_when_sample_taken, lq_tests_conducted, hq_tests_conducted, tests_conducted, results_declared, locations_moved, hcw_recovered};
  // Write output to a FILE

  char str[210];

  sprintf(str,"./Targeted-Testing_BeginAt_%g_DTR_%g_RAT_%g_%g_PCR_%g_%g_%lf%lf-%i.txt",begin_at,test_frac,Tpars[0][0],Tpars[0][3],Tpars[1][0],Tpars[1][3],uniform(),uniform(),iter);
  FILE *fpt=(FILE *)fopen(str,"wt"); // Open the file to print output

  fprintf(fpt,"###### TEST LOG ####################\n");
  fprintf(fpt,"# Time taken               : %.2f s\n",time_taken);

  fprintf(fpt,"# Test Parameters: \n");
  fprintf(fpt,"# %.2f %.2f %.2f %.2f\n", Tpars[0][0], Tpars[0][1], Tpars[0][2], Tpars[0][3]);
  fprintf(fpt,"# %.2f %.2f %.2f %.2f\n", Tpars[1][0], Tpars[1][1], Tpars[1][2], Tpars[1][3]);

  fprintf(fpt,"# Confined Less Infective? : %s\n", details[0] ? "True": "False");
  fprintf(fpt,"# Homes Quarantined?       : %s\n", details[1] ? "True": "False");
  fprintf(fpt,"# Quarantine when sampled? : %s\n", details[2] ? "True": "False");
  fprintf(fpt,"# Testing Started When     : %.2f%% recovered\n", begin_at);
  fprintf(fpt,"# Fraction Tested Daily    : %.2f%%\n", test_frac);

  fprintf(fpt,"# LQ Tests Done in total   : %d\n", details[3]);
  fprintf(fpt,"# HQ Tests Done in total   : %d\n", details[4]);
  fprintf(fpt,"# All Tests Done in total  : %d\n", details[5]);
  fprintf(fpt,"# Results Given in total   : %d\n", details[6]);
  fprintf(fpt,"# Locations Moved in total : %d\n", details[7]);
  fprintf(fpt,"# Total recovered HCW      : %d\n", details[8]);
  fprintf(fpt,"# Total HCW                : %d\n", details[9]);

  fprintf(fpt,"# Rate Array: \n");
  for(int i=0; i<n_states;i++){fprintf(fpt,"# ");for(int j=0; j<n_states;j++){fprintf(fpt,"%5g ", rate[i][j]);}fprintf(fpt,"\n");}
  fprintf(fpt,"###### END LOG #####################\n");
  fprintf(fpt,"#\n");



  for(int i=0;i<tf;i++){
    for(int j=0; j<op_width;j++){
      fprintf(fpt, "%i ", output[i][j]);
    }
    fprintf(fpt, "\n");
  }
  fflush(fpt);  // Flush the buffer (Slows down by roughly 4%)
}


void Targetted_Run(double dt, double Tpars[][4], int tf, bool lock_homes, bool quarantine_when_sample_taken, double begin_at, double test_frac,int iter){

  clock_t start, end; // Measuring how long the function takes to run
  double cpu_time_used;
  start = clock();

  bool is_confined[n_pop]  = {};
  bool being_tested[n_pop] = {};
  bool loc_confined[n_loc] = {};

  //* Added to keep track of retesting RAT negative symptomatics **//
  int last_test_type[n_pop]={};
  std::fill_n(last_test_type,n_pop,-1000);  // Set all values of last_test_type to -1000
  int last_test_result[n_pop]={};
  std::fill_n(last_test_result,n_pop,-1000);  // Set all values of last_test_result to -1000
  /*****************************************************************/

  int test_result[n_pop]   = {};
  int next_test_date[n_pop]= {};

  int result_declared_date[n_pop]={};
  std::fill_n(result_declared_date,n_pop,-1000);  // Set all values of result_declared_date to -1000

  int loc_confined_time[n_loc] = {};
  std::fill_n(loc_confined_time,n_loc,-1000);     // Set all values of loc_confined_time to -1000

  int person_isolated_time[n_pop] = {};
  std::fill_n(person_isolated_time, n_pop,-1000);


  int hq_tests_conducted = 0;
  int lq_tests_conducted = 0;
  int tests_conducted    = 0;
  int results_declared   = 0;
  int locations_moved    = 0;

  double r[n_states][n_states] = {};                      // Array to store rates per event

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
  bool midday_move_completed = false;

  int output[tf+1][op_width];   // Output array, to be printed to file

  // First line of output
  output[day][0] = day; for(int s=0;s<n_states;s++){output[day][s+1]=n[s];} output[day][n_states+1]=hq_tests_conducted;output[day][n_states+2]=lq_tests_conducted;output[day][n_states+3]=tests_remaining_today;

  while(t<tf){

    // Moving people around deterministically (HOME TO WORK)
    if(midday_move_completed == false && t - day >  0.5){
      midday_move_completed = true;

      for(int i=0;i<n_pop;i++){

        if(pop[i][0]!=H && is_confined[i] == false && loc_confined[pop[i][3]]==false){
          locations_moved++;
          int home_loc = pop[i][1];
          int work_loc = pop[i][2];
          if(pop[i][3]==home_loc){ pop[i][3] = work_loc; n_per_location[home_loc][pop[i][0]]--; n_per_location[work_loc][pop[i][0]]++; }
        }
      }
    }

    for(int i=0; i< n_loc; i++){
      printf("Location %i",i);
      int N = 0;  for(int j=0; j<n_states; j++){ N += n_per_location[i][j];} // Total number currently in this location according to n_per_location

      if(N==0){continue;}

      int backup[n_states]={};
      printf("N %i",N);
      int ind[N]={};  // Indices of people in this location

      int counter = 0;
      for(int j=0; j<n_pop; j++){       // Loop over all people_linked_to (not pop!)

        if(people_linked_to[i][j]==-1){break;}      // If you reach the "end" of people_linked_to, break
        else if(pop[people_linked_to[i][j]][3]==i){ // Otherwise, if the person is currently here
          ind[counter] = people_linked_to[i][j];
          backup[pop[ind[counter]][0]]++;
          counter++;
        }
      }

      if(counter!=N){printf("Error in location %i, N = %i, counter = %i \n",i,N,counter );}

      // bool temp_flag=false;
      // for(int j=0; j<n_states;j++){if(n_per_location[i][j]-backup[j]!=0)temp_flag=true;}  // TEST: To remove
      // if(temp_flag==true){printf("Error in location %i\n",i); printf("N_per_loc\n");for(int k =0; k<n_states; k++){printf("%i ",n_per_location[i][k]);}printf("\n");printf("Backup\n");for(int k =0; k<n_states; k++){printf("%i ",backup[k]);}printf("\n");}

      int conf_by_state_in_loc[n_states]={};
      int tot_conf = 0;
      if(quarantine_confined == true){
        for(int j = 0; j<N;j++){      // Go over people in location
            int p = ind[j];           // index of person
            if(is_confined[p] == true){conf_by_state_in_loc[pop[p][0]]++; tot_conf++; } // If the person is confined, increment the conf_by_state_in_loc of their state
        }
      }

      if(loc_confined[i]==true && conf_by_state_in_loc[S]+conf_by_state_in_loc[A]+conf_by_state_in_loc[P]+conf_by_state_in_loc[MI]+conf_by_state_in_loc[SI]+conf_by_state_in_loc[R]+conf_by_state_in_loc[H] == 0){loc_confined[i]=false;loc_confined_time[i]=-1000;} // Unlock the house if there are no confined people


      shuffle(0,N,ind); // Shuffle list of people currently in locations

      // Rates that don't change
      // A->R
      r[A][R] = rate[A][R];
      // P->MI
      r[P][MI] = rate[P][MI];
      // P->SI
      r[P][SI] = rate[P][SI];
      // MI->R
      r[MI][R] = rate[MI][R];
      // SI->R
      r[SI][R] = rate[SI][R];
      // SI->H
      r[SI][H] = rate[SI][H];
      // H->R
      r[H][R] = rate[H][R] ;



      for(int j=0;j<N;j++){                    // Loop over the people currently in the location

        int newN = n_per_location[i][S]+n_per_location[i][A]+n_per_location[i][P]+n_per_location[i][MI]+n_per_location[i][SI]+n_per_location[i][R]+n_per_location[i][H];


        double V = newN - alphaH*n_per_location[i][H];   // Spatial damping parameter (adjusted by alphaH for hospitals, and alphaQ for homes where there are people confined)


        // S->A
       r[S][A] = rate[S][A]  * 1/V * (Cpars[0]*(n_per_location[i][A]- conf_by_state_in_loc[A]*alphaQ) +    // SA
                                      Cpars[1]*(n_per_location[i][P]- conf_by_state_in_loc[P]*alphaQ) +    // SP
                                      Cpars[2]*(n_per_location[i][MI]-conf_by_state_in_loc[MI]*alphaQ) +   // S-MI
                                      Cpars[3]*(n_per_location[i][SI]-conf_by_state_in_loc[SI]*alphaQ) +   // S-SI
                                      Cpars[4]*n_per_location[i][H]);                                      // SH

      // S->P
      r[S][P] = rate[S][P]   * 1/V * (Cpars[0]*(n_per_location[i][A]- conf_by_state_in_loc[A]*alphaQ) +    // SA
                                      Cpars[1]*(n_per_location[i][P]- conf_by_state_in_loc[P]*alphaQ) +    // SP
                                      Cpars[2]*(n_per_location[i][MI]-conf_by_state_in_loc[MI]*alphaQ) +   // S-MI
                                      Cpars[3]*(n_per_location[i][SI]-conf_by_state_in_loc[SI]*alphaQ) +   // S-SI
                                      Cpars[4]*n_per_location[i][H]);                                      // SH


        double exit_rate[n_states]={};                                    // Array to store total exit rates from a state
        for(int jj=0;jj<n_states;jj++){
          for(int k=0;k<n_states;k++){exit_rate[jj] += r[jj][k];}
        }



        int from = pop[ind[j]][0];
        if(uniform()<exit_rate[from]*dt){            // If the person is selected to move
          double p = uniform();
          double temp=0;

          r[S][A] = (rate[S][A]/V)*(n_per_location[i][A]+n_per_location[i][P]+n_per_location[i][MI]+n_per_location[i][SI] - 0.9*tot_conf + 0.1*n_per_location[i][H]);
          r[S][P] = (rate[S][P]/V)*(n_per_location[i][A]+n_per_location[i][P]+n_per_location[i][MI]+n_per_location[i][SI] - 0.9*tot_conf + 0.1*n_per_location[i][H]);
        
          // S->A
          r[S][A] = rate[S][A]  * 1/V * (Cpars[0]*(n_per_location[i][A]- conf_by_state_in_loc[A]*alphaQ) +    // SA
                                         Cpars[1]*(n_per_location[i][P]- conf_by_state_in_loc[P]*alphaQ) +    // SP
                                         Cpars[2]*(n_per_location[i][MI]-conf_by_state_in_loc[MI]*alphaQ) +   // S-MI
                                         Cpars[3]*(n_per_location[i][SI]-conf_by_state_in_loc[SI]*alphaQ) +   // S-SI
                                         Cpars[4]*n_per_location[i][H]);                                      // SH

          // S->P
          r[S][P] = rate[S][P]  * 1/V * (Cpars[0]*(n_per_location[i][A]- conf_by_state_in_loc[A]*alphaQ) +    // SA
                                         Cpars[1]*(n_per_location[i][P]- conf_by_state_in_loc[P]*alphaQ) +    // SP
                                         Cpars[2]*(n_per_location[i][MI]-conf_by_state_in_loc[MI]*alphaQ) +   // S-MI
                                         Cpars[3]*(n_per_location[i][SI]-conf_by_state_in_loc[SI]*alphaQ) +   // S-SI
                                         Cpars[4]*n_per_location[i][H]);                                      // SH



          for(int to=0;to<n_states;to++){    // Loop over possible "to" states
            temp += (r[from][to]/exit_rate[from]);
            if(p<temp){     // If such a transition must occur,
              pop[ind[j]][0] = to;           // Send this person to the "to" state.
              n_per_location[i][from]--; n_per_location[i][to]++;
              n[from]--; n[to]++;            // Change the values of n[from] and n[to]

              if(to==H){
                // Move them to the hospital

                  int h = randint(n_hospitals);                       // Choose a hospital at random
                  while(h==pop[ind[j]][2]){h = randint(n_hospitals);} // Make sure a HCW is not hospitalised in the same hospital

                  for(int s=0;s<n_pop;s++){                           // Link them to this hospital forever. (They are never unlinked, even after they leave)
                    if(people_linked_to[h][s] == -1){people_linked_to[h][s]=ind[j];people_linked_to[h][s+1]=-1;break;}
                  }

                  if(is_confined[ind[j]]==true){
                    is_confined[ind[j]] = false;           // Remove confinement
                    person_isolated_time[ind[j]] = -1000;  // Added new: reset time of isolation.
                    conf_by_state_in_loc[from]--;          // Reduce number of confined in this location if the person was confined
                  }

                  n_per_location[ pop[ind[j]][3] ][H]--; // Decrement number of "from" in current location
                  pop[ind[j]][3] = h;                    // Send them to a random hospitals
                  n_per_location[ pop[ind[j]][3] ][H]++; // Increment number of "to" in current location

              }
              else if(from==H && to == R){
                // Remove confinement and move them home
//                  is_confined[ind[j]] = false; // Remove confinement

                 // Send recovered who were hospitalised home
                n_per_location[ pop[ind[j]][3] ][R]--; // Decrement number of k (=R) in current locations
                pop[ind[j]][3] = pop[ind[j]][1];       // Send them home
                n_per_location[ pop[ind[j]][3] ][R]++; // Increment number of k (=R) in current locations
              }
              break; // Exit the "to" loop, move to next person
            } // End if condition (if p<temp)
          } // End loop over "to" states
        } // End if condition (infection state == from)
      } // End loop over people in location



      /************ END CHANGE STATE OF POP ***********/


      // if(n_per_location[i][S] <0 || n_per_location[i][A]<0 || n_per_location[i][P]<0 || n_per_location[i][MI]<0 ||
      //    n_per_location[i][SI]<0 || n_per_location[i][R]<0 || n_per_location[i][H]<0){ printf("HORRIBLE, HORRIBLE ERROR! NEGATIVE POPULATIONS!\n");}

    }
    // End of location loop.

    /**************************** THE END OF DAY(S) ****************************/

    if(t>=day+1){


            /************* TARGETTED TESTING ******************/

            if(n[R]>=begin_at/100 * n_pop){

              // Find all symptomatics

              int list_of_sym[n_pop]={};
              int n_sym = 0;

              int sym_rat_neg[n_pop]={};            // Symptomatics who tested negative on an RAT test that was declared greater than or equal to 7 days ago
              int n_srn = 0;

              int list_of_remaining[n_pop]={};       // Array to store remaining people for random testing
              int n_remaining  = 0;                  // Number of such people

              int not_eligible_for_testing = 0;      // Number of people not eligible for testing today

              for(int i=0; i<n_pop;i++){

                /******************** Find list of targets and remaining people **********************/

                if((pop[i][0]==MI||pop[i][0]==SI) && day>=next_test_date[i] && being_tested[i]==false && is_confined[i]==false){
                  list_of_sym[n_sym] = i; n_sym++;

                  // NEW: (CHECK!!!) Retesting symptomatics (possible problem: they may not have been SI or MI when they were tested a week ago!)
                  if(day>=result_declared_date[i]+7 && last_test_type[i]==0 && last_test_result[i]==-1){sym_rat_neg[n_srn] = i; n_srn++;} // Make separate list of symptomatics who
                                                                                                                                          // last tested negative on a RAT, more than 7 days ago
                }
                else if(pop[i][0]!=H && being_tested[i]==false && is_confined[i]==false){
                  list_of_remaining[n_remaining] = i; n_remaining++;
                }
                else{ not_eligible_for_testing++;}

                /****************************** Done finding lists **********************************/
              }

              list_of_sym[n_sym] = -1;                 // Set -1 to mark the end of this array
              list_of_remaining[n_remaining] = -1;     // Set -1 to mark end of this array

              sym_rat_neg[n_srn] = -1;                 // Set -1 to mark end of this Array


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

                // list_of_remaining[n_remaining] = -1;  // Reset position of -1 to mark new end of this array

                  /*** REARRANGE list_of_sym SO THAT sym_rat_neg COME FIRST! *****/
                  int counter = 0;

                  for(int j=0;j<n_srn;j++){                 // For every sym_rat_neg
                    for(int k=0;k<n_sym;k++){               // Go over the remaining symptomatics
                      if(list_of_sym[k] == sym_rat_neg[j]){ // When you find this symptomatic in that list
                        int temp = list_of_sym[counter];       // Swap the lowest person on the list eligible
                        list_of_sym[counter] = list_of_sym[k];
                        list_of_sym[k] = temp;

                        counter++; break;                      // Increase the counter and break this loop
                      }
                    }
                  }
                  // Commented out lines below just to check whether it behaves as expected (it does!)
                  // printf("sym_rat_negs :"); for(int j=0;j<n_srn;j++){printf("%i ",sym_rat_neg[j]);}printf("\n");
                  // printf("Individual details :\n"); for(int j=0;j<n_srn;j++){printf("%i: Today %i Result declared %i  last_test_type %i last_test_result %i\n",sym_rat_neg[j], day, result_declared_date[sym_rat_neg[j]],last_test_type[sym_rat_neg[j]], last_test_result[sym_rat_neg[j]]);}printf("\n");
                  // printf("Pos rearranged:"); for(int j=0;j<n_sym;j++){printf("%i ",list_of_sym[j]);}printf("\n");

                for(int j=0; j<targetted_tests_done_today;j++){

                  int si = list_of_sym[j];  // Individual to test

                  if(day >= next_test_date[si] && being_tested[si]==false && tests_remaining_today>0){
                    // If so, perform a tests
                    being_tested[si] = true; tests_conducted++; tests_remaining_today--;

                    if(quarantine_when_sample_taken==true){   // Quarantine as soon as sample is taken
                      // Move them Home
                      n_per_location[pop[si][3]][pop[si][0]]--;
                      pop[si][3] = pop[si][1];
                      n_per_location[pop[si][3]][pop[si][0]]++;

                      is_confined[si] = quarantine_confined;
                      person_isolated_time[si] = day;             // Added new: time of isolation.

                      loc_confined[pop[si][1]] = lock_homes;      // Lock home depending on variable `lock_homes`.
                      loc_confined_time[pop[si][1]] = day;
                    }

                    // Targetted testing using HQ tests unless there are none
                    int test_type = 1;

                    if(hq_tests_today<=0){test_type = 0;} // If no HQ tests available, give them LQ (one of the two is guaranteed, since tests_remaining_today>0)

                    last_test_type[si] = test_type;       // NEW: Keeping track of last test type

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

                    // ** THE PART BELOW EXCEPT FOR SELECTING TESTS IS ESSENTIALLY THE SAME AS FOR TARGETTED TESTING (with si -> ri)

                    being_tested[ri] = true; tests_conducted++; tests_remaining_today--;

                    if(quarantine_when_sample_taken==true){     // Quarantine as soon as sample is taken
                      // Move them Home
                      n_per_location[pop[ri][3]][pop[ri][0]]--;
                      pop[ri][3] = pop[ri][1];
                      n_per_location[pop[ri][3]][pop[ri][0]]++;

                      is_confined[ri] = quarantine_confined;
                      person_isolated_time[ri] = day;             // Added new: time of isolation.

                      loc_confined[pop[ri][1]] = lock_homes;      // Lock home depending on variable `lock_homes`.
                      loc_confined_time[pop[ri][1]] = day;

                    }

                    int test_type = randint(2);        // Returns either 0 and 1 with equal probability.

                    last_test_type[ri] = test_type;       // NEW: Keeping track of last test type

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

                  if(day == result_declared_date[i] && being_tested[i]==true){

                    // First set them to no longer being tested
                    being_tested[i] = false;
                    results_declared++;

                    // Declare results_declared

                    if(test_result[i]==1 && pop[i][0] != H){  // If the result is positive, and the person hasn't already moved to a hospital
                      positives++;
                      is_confined[i] = quarantine_confined; // Confine them if quarantine_confined==true
                      person_isolated_time[i] = day;        // Added new: time of isolation.

                      last_test_result[i] = 1;              // NEW: Added to keep track of retesting RAT negative symptomatics

                      // Move them home
                      n_per_location[pop[i][3]][pop[i][0]]--;              // Decremement number in current location
                      if(pop[i][3] == pop[i][2]){ pop[i][3] = pop[i][1];}  // If they're at work, send them home.
                      n_per_location[pop[i][3]][pop[i][0]]++;              // Incremement number in current location.

                      loc_confined[pop[i][1]] = lock_homes;                      // Reconfine their homes
                      loc_confined_time[pop[i][1]] = day;

                      next_test_date[i] = day + if_positive_test_after; // TO make compatible with applet. (CHECK!) NEW: changed from += 14 to today + 14 days.
                    }
                    else if(test_result[i] == -1 && pop[i][0] != H){       // If they're negative and have not been hospitalised (CHANGED!)
                      is_confined[i] = false;  // Remove confinement
                      person_isolated_time[i] = -1000;        // Reset time of isolation.

                      last_test_result[i] = -1;              // NEW: Added to keep track of retesting RAT negative symptomatics

                      // if(last_test_type[i]==0){next_test_date[i]=day+7;} // NEW: (CHECK!!) If they were tested negative using an RAT test (THIS IS PROBABLY WRONG, and not needed)
                    }

                    if(pop[i][0]==H){is_confined[i]=false; person_isolated_time[i]=-1000;last_test_result[i]=0;} // If they've been hospitalised, do the same as for negative results,
                                                                                                                 // only set their last test result to 0.(CHECK!!!!)

                    test_result[i] = 0;        // Reset test result to 0.
                  }

                  /***********************DONE DECLARING!**********************/

                  /******** Remove confinement if 14 days have passed ********/

                  if(is_confined[i]==true && day >= person_isolated_time[i]+total_isolation_time){ // NEW: changed 14 -> total_isolation_time
                    is_confined[i]=false;
                    person_isolated_time[i] = -1000;
                  }

              }

            }

            // int unlockedtoday = 0;
            // Lock or unlock homes
            for(int i=0; i<n_loc;i++){    // This can probably be added to the first location loop, if needed.

              if(loc_confined[i]==true && day-loc_confined_time[i] >= total_loc_confined_time){
                // for(int j = 0; j<n_pop; j++){if(pop[j][3]==i && is_confined[j]==true){printf("Problem in location %i\n", i);}}
                loc_confined[i]      = false; // Remove confinement
                loc_confined_time[i] = - 1000;
                // unlockedtoday++;
              }

              // Work this in a while
              // if(loc_confined[i]==true && conf_by_state_in_loc[S]+conf_by_state_in_loc[A]+conf_by_state_in_loc[P]+conf_by_state_in_loc[MI]+conf_by_state_in_loc[SI]+conf_by_state_in_loc[R]+conf_by_state_in_loc[H] == 0){loc_confined[i]=false;loc_confined_time[i]=-1000;} // Unlock the house if there are no confined people
            }


      // Moving people around deterministically (WORK TO HOME)

      for(int i=0;i<n_pop;i++){
        if(pop[i][0]!=H && is_confined[i] == false && loc_confined[pop[i][3]]==false){

          locations_moved++;

          int home_loc = pop[i][1];
          int work_loc = pop[i][2];
          if(pop[i][3]==work_loc){ pop[i][3] = home_loc; n_per_location[work_loc][pop[i][0]]--; n_per_location[home_loc][pop[i][0]]++; }
        }
      }

      midday_move_completed = false;

      // Increment the day, write the output to an array, and reset the number of tests  //
      day++;

      output[day][0] = day; for(int s=0;s<n_states;s++){output[day][s+1]=n[s];} output[day][n_states+1]=hq_tests_conducted;output[day][n_states+2]=lq_tests_conducted;output[day][n_states+3]=tests_remaining_today;
      tests_remaining_today = tests_available_daily;
      lq_tests_today = lq_tests_daily;
      hq_tests_today = hq_tests_daily;

      // int locked = 0;                                                         // TEST: to remove
      // for(int hq=0;hq<n_loc;hq++){if(loc_confined_time[hq]>=0){locked++;}}    // TEST: to remove

      // printf("%i Confined: %i  Confined homes: %i Unlocked today: %i\n", day, sum(is_confined,n_pop), sum(loc_confined,n_loc) ,unlockedtoday);  // TEST: to remove
    }

    /********************* END OF END OF DAYS ********************************/


    t+=dt;
  }  // End of While Loop above.

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

  int details[] = {quarantine_confined, lock_homes, quarantine_when_sample_taken, lq_tests_conducted, hq_tests_conducted, tests_conducted, results_declared, locations_moved, hcw_recovered,hcw};
  writetofile(output, tf, Tpars, begin_at, test_frac,cpu_time_used,details,iter);



}// End of TargettedTesting function.



int main() {

  // S, A, P, MI, SI, R, H
  // 0  1  2  3   4   5  6

  double k_S   = 0.25;
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

  // exit_rate[S] = k_S;
  // exit_rate[A] = k_A;
  // exit_rate[P] = k_P;
  // exit_rate[MI]= k_MI;
  // exit_rate[SI]= k_SI;
  // exit_rate[R] = 0;
  // exit_rate[H] = k_H;


  old_dt = 0.01;
  double dt = old_dt;

  double Tpars[2][4] = {{0.6, 0.98, 0, 0},
                        {1.0,  1.0, 0, 1}};

  int tf = 200;
  bool lock_homes=false;
  bool quarantine_when_sample_taken=false;
  double begin_at = 20;
  double test_frac=0.0;


 // double p[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
 double p[] = {0.0};
 int p_len  = sizeof(p)/sizeof(p[0]);

 // double s[] = {0.5, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95};
 double s[] = {0.50};
 int s_len  = sizeof(s)/sizeof(s[0]);

// double dtr[]={0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
 double dtr[]={0.5};

 int dtr_len = sizeof(dtr)/sizeof(dtr[0]);

  int mc_runs = 1;



  for(int i=0;i<mc_runs;i++){
   for(int j=0;j<s_len;j++){

     Tpars[0][0] = s[j];

     for(int k=0; k<p_len;k++){

       Tpars[0][3] = p[k];
       Tpars[1][3] = 1-p[k];

       for(int l =0; l<dtr_len;l++){

         test_frac = dtr[l];
            createPopulation();
            Targetted_Run(dt,
                        Tpars,
                        tf,
                        lock_homes,
                        quarantine_when_sample_taken,
                        begin_at,
                        test_frac,
                        i);

       }
     }
   }
  }

}
