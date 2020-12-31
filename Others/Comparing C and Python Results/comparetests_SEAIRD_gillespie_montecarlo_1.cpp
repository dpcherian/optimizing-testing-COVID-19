// Overall logic: location network + infection model + testing/tracing/isolation

// Location network:
// This consists of NL nodes and NE edges
// nodes correspond to locations, such as homes, work locations, schools etc
// and people move between locations that are connected by the edges
// Each person has one home location and one non-home location between which they move at some rate
// Thus each edge (or link - I use the terms interchangeably) is specied by two numbers: a person ID and a rate

// Infection model:
// Each person can be in state S,E,A,I,R or D
// Each transition is specified by a rate

// Testing/tracing/isolation:
// At the moment a person is tested at a certain rate if they are not already confined or not already dead
// Testing to declaring result is also modelled as a Poisson process with a specified rate (probably a bad 
// approximation, but easier to implement in Gillespie)
// At the moment the test has a specified sensitivity (can give false negative) and 
// given specificity (cannot give false positive) 
// If a person is declared positive, they are immediately moved to their home and confined there
// However, infection can continue, i.e. if another person at the same home location is susceptible, 
// they could get infected

// TODO:
// 1. build different kinds of (nonrandom) networks
// 2. add age to each person, all infection model parameters will also become age dependent
// 3. implement some simple contact tracing, or at least confine everyone in a particular home if one person 
// from that home tests positive
// 4. modify code to allow using two different tests, e.g. one with high sensitivity but long delay and the
// other with lower sensitivity but less delay to declaring the result

// some necessary headers - ignore
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

// random number generator -------------------------------------------------
// (copied from Numerical Recipes in C)
// not critical - any decent random number generator can be used

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
int j;
long k;
static long iy=0;
static long iv[NTAB];
float temp;

if (*idum<=0 || !iy)
	{
	if (-(*idum)<1) *idum=1;
	else *idum=-(*idum);
	for (j=NTAB+7;j>=0;j--)
		{
		k=(*idum)/IQ;
		*idum=IA*(*idum-k*IQ)-IR*k;
		if (*idum<0) *idum+=IM;
		if (j<NTAB) iv[j]=*idum;
		}
	iy=iv[0];
	}

k=(*idum)/IQ;
*idum=IA*(*idum-k*IQ)-IR*k;
if (*idum<0) *idum+=IM;
j=iy/NDIV;
iy=iv[j];
iv[j]=*idum;
if ((temp=AM*iy)>RNMX) return RNMX;
else return temp;
}

// -------------------------------------------------------------------------
long idum1=-13; // set random number generator seed
// -------------------------------------------------------------------------

const int Numstates=6,NL=2750,NP=10000,NE=NP; // number of infection states, number of locations, number of people, total number of links(edges in network) - right now the simulation won't work if NE is not the same as NP, code needs to be generalized
int links[NE][3]; // array of links/edges [i][0]&[i][1] are the two locations the ith edge connects, [i][2] is the person ID for that link
double linkweight[NE]; // the rate at which each link is used by the corresponding person
double infectedtime[NP]; // time at which a person first got infected (not using this anywhere right now, but may be useful later)
double rate[Numstates][Numstates]; // rate from one state to another
double backgroundinfectionrate=0; // (small) chance for S->E even without being near an I person - this is to model a small rate at which people meet someone outside their usual network and get infected
double randomtestrate=1000,testsensitivity=0.99,testspecificity=1,finishtestrate=24.0; // rate of testing people, test sensitivity, rate at which result is declared
int infectionstate[NP]; // state of each person; 0=S;1=E;2=A;3=I;4=R;5=D
int testingstate[NP]; // 0=not being tested; -1 or 1 = being tested (result not out yet) but will turn out -ve or +ve; 2=test result came +ve and is quarantined
int homeloc[NP],currentloc[NP]; // home location for each person (doesn't change once set at the beginning of the simulation), current location for each person (changes as people move along links/edges)
int numinstateinloc[Numstates][NL];// number of people who are a given state in each location 
int numinstate[Numstates]; // number of people in a given state (across all locations)

void initialize_location_network(void) // sets up location network
{
// set up links[i][0..2] and linkweight[i] and homeloc[i]
// arbitrary choice for now, needs to be generalized to build other kinds of (non-random) networks
for (int i=0;i<NP;i++) {homeloc[i]=(int)(2500*ran1(&idum1));links[i][0]=homeloc[i];links[i][1]=250+(int)(2500*ran1(&idum1));links[i][2]=i;linkweight[i]=1;}
}

void initialize_infection(void)
{
int i,j;

// initialize rates
rate[0][1]=0.2;
rate[1][2]=0.2;
rate[2][3]=0.5;
rate[2][4]=0.1;
rate[3][4]=0.05;
rate[3][5]=0.05;

// initialize all to be susceptible except person ten random person, who are marked infectious (I=2)
for (i=0;i<NL;i++) for (j=0;j<Numstates;j++) numinstateinloc[j][i]=0;
for (i=0;i<NP;i++) {testingstate[i]=infectionstate[i]=0;infectedtime[i]=-1000;currentloc[i]=homeloc[i];numinstateinloc[0][currentloc[i]]++;}
for (j=0;j<Numstates;j++) numinstate[j]=0;
numinstate[0]=NP;
for (i=0;i<10;i++)
	{
	do j=(int)(NP*ran1(&idum1)); while (infectionstate[j]>0);
	infectionstate[j]=2;
	numinstate[0]--;numinstate[2]++;
	numinstateinloc[0][currentloc[j]]--;numinstateinloc[2][currentloc[j]]++;
	}

}

double infectionmodel_gillespierun(double tf) // runs the Gillespie simulation
{
int i,j;
int newloc;
int numStoE,numEtoA,numAtoI,numAtoR,numItoR,numItoD,nummove;

double t,deltat,tout;
double psum,a0,a,a1,a2,a3,a4,a5,a6,a7,a8;
char str[200];

sprintf(str,"gillespie_highdensity_test_%d_%d_%dhr.txt",(int)(randomtestrate*100),(int)(testsensitivity*100),(int)(24.0/finishtestrate));
FILE *fpx,*fpt=(FILE *)fopen(str,"wt"); // for output

initialize_infection();
t=0;tout=0; // t tracks current time, tout tracks the time at which to print data to output file
numStoE=numEtoA=numAtoI=numAtoR=numItoR=numItoD=nummove=0;

while (numinstate[1]+numinstate[2]+numinstate[3]>0) // run the infection until it completely dissipates, alternatively could use while (t>tf) to run till a specified final time
//while (t<1)
	{
	if (t>=tout) {fprintf(fpt,"%lf ",t);for (j=0;j<Numstates;j++) fprintf(fpt,"%d ",numinstate[j]);fprintf(fpt,"%d %d %d %d %d %d %d\n",numStoE,numEtoA,numAtoI,numAtoR,numItoR,numItoD,nummove);tout+=1;numStoE=numEtoA=numAtoI=numAtoR=numItoR=numItoD=nummove=0;} // print data to output file

	a0=0; // a0 will keep track of the sum of all rates for the Gillespie
	for (i=0;i<NP;i++) if (infectionstate[i]==0) a0+=rate[0][1]*(numinstateinloc[2][currentloc[i]]+numinstateinloc[3][currentloc[i]])+backgroundinfectionrate;// S->E
	a1=a0;
	for (i=0;i<NP;i++) if (infectionstate[i]==1) a0+=rate[1][2];//E->A
	a2=a0;
	for (i=0;i<NP;i++) if (infectionstate[i]==2) a0+=rate[2][3];//A->I
	a3=a0;
	for (i=0;i<NP;i++) if (infectionstate[i]==2) a0+=rate[2][4];//A->R
	a4=a0;
	for (i=0;i<NP;i++) if (infectionstate[i]==3) a0+=rate[3][4];//I->R
	a5=a0;
	for (i=0;i<NP;i++) if (infectionstate[i]==3) a0+=rate[3][5];//I->D
	a6=a0;
	
	for (i=0;i<NP;i++) if (testingstate[i]==0 && infectionstate[i]!=5) a0+=randomtestrate;//test a person if they are not already under testing or confinement and they are not dead
	a7=a0;
	for (i=0;i<NP;i++) if (testingstate[i]==1 || testingstate[i]==-1) a0+=finishtestrate; // declare test result
	a8=a0;
	for (i=0;i<NE;i++) if (testingstate[links[i][2]]!=2 && infectionstate[links[i][2]]!=5) a0+=linkweight[i]; // move people along links if they are not confined and they are not dead
	deltat=-log(ran1(&idum1))/a0; // pick time step from exponential distribution with mean 1/a0

	// now choose which transition happens at time t+deltat
	a=ran1(&idum1)*a0;
	if (a<a1) // S->E
		{
		psum=0;
		for (i=0;i<NP;i++) if (infectionstate[i]==0) {psum+=rate[0][1]*(numinstateinloc[2][currentloc[i]]+numinstateinloc[3][currentloc[i]])+backgroundinfectionrate;if (psum>a) break;}
		// at this point i is the individual who needs to be exposed
		if (infectionstate[i]!=0) {printf("error picking S or I or R to become infectious %lf %d %lf %lf %lf %lf",t,i,a0,a1,a2,a);exit(0);}
		infectionstate[i]=1;infectedtime[i]=t+deltat;
		numinstateinloc[0][currentloc[i]]--;numinstateinloc[1][currentloc[i]]++;
		numinstate[0]--;numinstate[1]++;
		numStoE++;
		}
	else if (a<a2) // E->A
		{
		psum=a1;	
		for (i=0;i<NP;i++) if (infectionstate[i]==1) {psum+=rate[1][2];if (psum>a) break;}
		if (infectionstate[i]!=1) {printf("error picking S or I or R to become infectious %lf %d %lf %lf %lf %lf",t,i,a0,a1,a2,a);exit(0);}
		infectionstate[i]=2;
		numinstateinloc[1][currentloc[i]]--;numinstateinloc[2][currentloc[i]]++;
		numinstate[1]--;numinstate[2]++;
		numEtoA++;
		}
	else if (a<a3) // A->I
		{
		psum=a2;
		for (i=0;i<NP;i++) if (infectionstate[i]==2) {psum+=rate[2][3];if (psum>a) break;}
		if (infectionstate[i]!=2) {printf("error picking S or E or R to recover\n");exit(0);}
		infectionstate[i]=3;
		numinstateinloc[2][currentloc[i]]--;numinstateinloc[3][currentloc[i]]++;	
		numinstate[2]--;numinstate[3]++;	
		numAtoI++;
		}
	else if (a<a4) // A->R
		{
		psum=a3;
		for (i=0;i<NP;i++) if (infectionstate[i]==2) {psum+=rate[2][4];if (psum>a) break;}
		if (infectionstate[i]!=2) {printf("error picking S or E or R to recover\n");exit(0);}
		testingstate[i]=0;
		infectionstate[i]=4;
		numinstateinloc[2][currentloc[i]]--;numinstateinloc[4][currentloc[i]]++;	
		numinstate[2]--;numinstate[4]++;	
		numAtoR++;
		}
	else if (a<a5) // I->R
		{
		psum=a4;
		for (i=0;i<NP;i++) if (infectionstate[i]==3) {psum+=rate[3][4];if (psum>a) break;}
		if (infectionstate[i]!=3) {printf("error picking S or E or R to recover\n");exit(0);}
		testingstate[i]=0;
		infectionstate[i]=4;
		numinstateinloc[3][currentloc[i]]--;numinstateinloc[4][currentloc[i]]++;	
		numinstate[3]--;numinstate[4]++;	
		numItoR++;
		}
	else if (a<a6) // I->D
		{
		psum=a5;
		for (i=0;i<NP;i++) if (infectionstate[i]==3) {psum+=rate[3][5];if (psum>a) break;}
		if (infectionstate[i]!=3) {printf("error picking S or E or R to recover\n");exit(0);}
		infectionstate[i]=5;
		numinstateinloc[3][currentloc[i]]--;numinstateinloc[5][currentloc[i]]++;	
		numinstate[3]--;numinstate[5]++;	
		numItoD++;
		}
		
	else if (a<a7) // (random) test
		{	
		psum=a6;
		for (i=0;i<NP;i++) if (testingstate[i]==0 && infectionstate[i]!=5) {psum+=randomtestrate;if (psum>a) break;}
		if (((infectionstate[i]==2 || infectionstate[i]==3) && ran1(&idum1)<testsensitivity) || ((infectionstate[i]<=1 || infectionstate[i]==4) && ran1(&idum1)>testspecificity)) testingstate[i]=1;
		else testingstate[i]=-1;
		}
	else if (a<a8) // declare result and isolate if positive
		{
		psum=a7;
		for (i=0;i<NP;i++) if (testingstate[i]==1 || testingstate[i]==-1) {psum+=finishtestrate;if (psum>a) break;}
		testingstate[i]++;
		if (testingstate[i]==2 && currentloc[i]!=homeloc[i])
			{
			numinstateinloc[infectionstate[i]][currentloc[i]]--;numinstate[infectionstate[i]]--;
			currentloc[i]=homeloc[i];
			numinstateinloc[infectionstate[i]][currentloc[i]]++;numinstate[infectionstate[i]]++;	
			}
		}
	else // move people along links
		{
		psum=a8;
		for (i=0;i<NE;i++) if (testingstate[links[i][2]]!=2 && infectionstate[links[i][2]]!=5) {psum+=linkweight[i];if (psum>a) break;}
		if (links[i][0]==currentloc[links[i][2]]) newloc=links[i][1];
		else newloc=links[i][0];
		// p.s. this is why the simulation will not run properly if NE (number of edges) is not equal to NP (number of people) - the code assumes a person is definitely at one of the nodes of any link that is labelled by that person's ID. If NE>NP this may not be true.
		numinstateinloc[infectionstate[links[i][2]]][currentloc[links[i][2]]]--;
		numinstateinloc[infectionstate[links[i][2]]][newloc]++;
		currentloc[links[i][2]]=newloc;
		nummove++;
		}

	t+=deltat;
	
	} // end of while loop

// output data for final time
fprintf(fpt,"%lf ",t);for (j=0;j<Numstates;j++) fprintf(fpt,"%d ",numinstate[j]);
fprintf(fpt,"%d %d %d %d %d %d %d\n",numStoE,numEtoA,numAtoI,numAtoR,numItoR,numItoD,nummove);
fclose(fpt);
return t;
}

void infectionmodel_montecarlorun(double tf)
{
int i,j;
int newloc;
int nummove,numtests,numtestsfinished;

double t,deltat,tout;
double r;

char str[200];

sprintf(str,"montecarlo_highdensity_test_%d_%d_%dhr.txt",(int)(randomtestrate*100),(int)(testsensitivity*100),(int)(24.0/finishtestrate));
FILE *fpx,*fpt=(FILE *)fopen(str,"wt"); // for output

initialize_infection();
deltat=0.01;
tout=0;nummove=numtests=numtestsfinished=0;

for (t=0;t<tf;t+=deltat)
	{
	if (t>=tout-1e-6) {fprintf(fpt,"%lf ",t);for (j=0;j<Numstates;j++) fprintf(fpt,"%d ",numinstate[j]);fprintf(fpt,"%d %d %d\n",nummove,numtests,numtestsfinished);tout+=1;nummove=numtests=numtestsfinished=0;} // print data to output file

	// infection transitions
	for (i=0;i<NP;i++)
		switch(infectionstate[i])
			{
			case 0:// Susc to Exposed
					if (ran1(&idum1)<rate[0][1]*(numinstateinloc[2][currentloc[i]]+numinstateinloc[3][currentloc[i]])*deltat) 
						{
						infectionstate[i]=1;infectedtime[i]=t;
						numinstateinloc[0][currentloc[i]]--;numinstateinloc[1][currentloc[i]]++;
						numinstate[0]--;numinstate[1]++;
						}
					break;
			case 1:// Exposed to Asymptomatic
					if (ran1(&idum1)<rate[1][2]*deltat) 
						{
						infectionstate[i]=2;
						numinstateinloc[1][currentloc[i]]--;numinstateinloc[2][currentloc[i]]++;
						numinstate[1]--;numinstate[2]++;
						}
					break;
			case 2:// Asymptomatic to Symptomatic or Recovered
					r=ran1(&idum1);
					if (ran1(&idum1)<rate[2][3]*deltat) 
						{
						infectionstate[i]=3;
						numinstateinloc[2][currentloc[i]]--;numinstateinloc[3][currentloc[i]]++;
						numinstate[2]--;numinstate[3]++;
						}
					else if (ran1(&idum1)<rate[2][4]*deltat)
						{
						infectionstate[i]=4;
						numinstateinloc[2][currentloc[i]]--;numinstateinloc[4][currentloc[i]]++;
						numinstate[2]--;numinstate[4]++;
						testingstate[i]=0;
						}	
					break; 
			case 3:// Symptomatic to Recovered or Dead
					r=ran1(&idum1);
					if (ran1(&idum1)<rate[3][4]*deltat) 
						{
						infectionstate[i]=4;
						numinstateinloc[3][currentloc[i]]--;numinstateinloc[4][currentloc[i]]++;
						numinstate[3]--;numinstate[4]++;
						testingstate[i]=0;
						}
					else if (ran1(&idum1)<rate[3][5]*deltat) 
						{
						infectionstate[i]=5;
						numinstateinloc[3][currentloc[i]]--;numinstateinloc[5][currentloc[i]]++;
						numinstate[3]--;numinstate[5]++;
						}	
					break;
			case 4:break; // Recovered
			case 5:break; // Dead
			}

	// testing
	// (i) declare results
	for (i=0;i<NP;i++) if (testingstate[i]==1 || testingstate[i]==-1) if (ran1(&idum1)<finishtestrate*deltat) 
		{
		testingstate[i]++;
		if (testingstate[i]==2 && currentloc[i]!=homeloc[i])
			{
			numinstateinloc[infectionstate[i]][currentloc[i]]--;numinstate[infectionstate[i]]--;
			currentloc[i]=homeloc[i];
			numinstateinloc[infectionstate[i]][currentloc[i]]++;numinstate[infectionstate[i]]++;	
			}
		numtestsfinished++;
		}
	// (ii) do tests
	for (i=0;i<NP;i++)
	  if (testingstate[i]==0 && infectionstate[i]!=5)
		if (ran1(&idum1)<randomtestrate*deltat)
			{
			numtests++;
			if (((infectionstate[i]==2 || infectionstate[i]==3) && ran1(&idum1)<testsensitivity) || ((infectionstate[i]<=1 || infectionstate[i]==4) && ran1(&idum1)>testspecificity)) testingstate[i]=1;
			else testingstate[i]=-1;
			}
		
	// move
	for (i=0;i<NE;i++) if (ran1(&idum1)<linkweight[i]*deltat && infectionstate[links[i][2]]!=5 && testingstate[links[i][2]]!=2)
		{
		nummove++;
		// move person along link i if they are not quarantined and not dead
		if (links[i][0]==currentloc[links[i][2]]) newloc=links[i][1];
		else newloc=links[i][0];
		// p.s. this is why the simulation will not run properly if NE (number of edges) is not equal to NP (number of people) - the code assumes a person is definitely at one of the nodes of any link that is labelled by that person's ID. If NE>NP this may not be true.
		numinstateinloc[infectionstate[links[i][2]]][currentloc[links[i][2]]]--;
		numinstateinloc[infectionstate[links[i][2]]][newloc]++;
		currentloc[links[i][2]]=newloc;
		}	
	
	} // end of for t loop

// output data for final time
fprintf(fpt,"%lf ",t);for (j=0;j<Numstates;j++) fprintf(fpt,"%d ",numinstate[j]);fprintf(fpt,"%d %d %d\n",nummove,numtests,numtestsfinished);
fclose(fpt);
}

int main(void)
{
int i;
double tf=170;

clock_t start, end, temp1, temp2; // this is just to measure how long the program takes to run
double cpu_time_used;
start = clock();
printf("started\n");

initialize_location_network();
//tf=infectionmodel_gillespierun(tf);
temp1 = clock();
cpu_time_used = ((double) (temp1 - start)) / CLOCKS_PER_SEC;
printf("time spent on gillespie %lf\n",cpu_time_used);

infectionmodel_montecarlorun(tf);

end = clock();
cpu_time_used = ((double) (end - temp1)) / CLOCKS_PER_SEC;
printf("time spent on monte carlo %lf\n",cpu_time_used);
cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
printf("total time spent %lf\n",cpu_time_used);
}
