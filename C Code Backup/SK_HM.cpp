// v5 vs v4: adding a list that keeps track of people in each location, so you can simply loop over them
// this is to make the logic easier to follow and modify (e.g. when we need to add age)

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// random number generator -------------------------------------------------

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

long idum=-1; // random number seed

// -----------------------------------------------------------------------

const int Maxtf=400;
int tf=200,steptime=1,currenttime;
double deltat=0.01; // time step

// SevenCompartment model parameters and stuff
int peopleperhome=4;
double WHratio=0.1;
const int Numstates=7,NP=10000,NE=NP,MaxNumloc=NP;
int Numhospitals=(int)floor(NP/5000),Numhomes=(int)floor(NP/peopleperhome),Numloc=(int)floor(Numhomes*(1+WHratio)); // number of infection states, number of locations, number of people, total number of links(edges in network) - right now the simulation won't work if NE is not the same as NP, code needs to be generalized
const int S=0,A=1,P=2,MI=3,SI=4,H=5,R=6; // 7 states

int links[NE][3]; // array of links/edges [i][0]&[i][1] are the two locations the ith edge connects, [i][2] is the person ID for that link
double linkweight[NE]; // the rate at which each link is used by the corresponding person
double infectedtime[NP]; // time at which a person first got infected (not using this anywhere right now, but may be useful later)
double exitrate[Numstates]; // rate of exit from each state
double Irecrate;
double pSP,pPSI,pSIH; // probability that S goes to P, P goes to SI and SI goes to H
double hdamping,qdamping; // multiplying factor in [0 1] that reduces transmissability in hospitals and quarantines
double totaltestrate=0.05;
double testsensitivity[2]={1.0, 0.6},testspecificity[2]={1.0, 0.98};
int testdelay[2]={0, 0}; // test sensitivity, rate at which result is declared
double RATfrac=0;
int numwaitingforresult_symp,numwaitingforresult_other,numtestedpos_symp,numtestedpos_other,testprob_symp,testprob_other;
int totalPCRleft,totalRATleft;
int testtime[NP],testtype[NP];
int totaltests;
double teststartfrac=0.2;

int qonstarttest; // whether to quarantine upon start of test or only when result comes out
int qhouse;
int infectionstate[NP],numinfectedby[NP],numinfectedbystate[Numstates]; // state of each person; 0=S;1=E;2=A;3=I;4=R;5=D
int testingstate[NP]; // 0=not being tested; -1 or 1 = being tested (result not out yet) but will turn out -ve or +ve; 2=test result came +ve and is quarantined
int homeloc[NP],currentloc[NP]; // home location for each person (doesn't change once set at the beginning of the simulation), current location for each person (changes as people move along links/edges)
int numquarantinedinloc[MaxNumloc],numinstateinloc[Numstates][MaxNumloc],numAtoSIquarantinedinloc[MaxNumloc];// number of people who are a given state in each location
int numinloc[MaxNumloc];
int firstinloc[MaxNumloc],previnloc[NP],nextinloc[NP];
int numinstate[Numstates]; // number of people in a given state (across all locations)
int maxinfected,maxinfectedtime;
int nummove,numtests,numtestsfinished;
int firstrun=1;
int teststarttime,teststartmarked;
int symplist[NP],otherlist[NP];
int numsymp,numother,numothertests,numsymptests,leftovertests;

int checkhealthofnumbers(void)
{
int i,j,k,l;
int n[Numstates],nn,nnn;
int nq[Numloc],nqASI[Numloc],nsl[Numstates][Numloc];
int errflag=0;

for (i=0;i<Numstates;i++) n[i]=0;
for (l=0;l<Numloc;l++)
  {
  for (i=0;i<Numstates;i++) n[i]+=numinstateinloc[i][l];
  }
nn=nnn=0;
for (i=0;i<Numstates;i++) {nn+=n[i];nnn+=numinstate[i];if (n[i]!=numinstate[i]) {errflag=1;printf("Error in numinstate %d on day %d (%d,%d)\n",i,currenttime,n[i],numinstate[i]);}}
if (nn!=NP) {errflag=1;printf("Error in numinstateinloc total %d on day %d\n",nn,currenttime);}
if (nnn!=NP) {errflag=1;printf("Error in numinstate total %d on day %d\n",nn,currenttime);}

// loop over all people and calc numquarantinedinloc and numAtoSIquarantinedinloc
/*for (l=0;l<Numloc;l++) nq[l]=nqASI[l]=0;
for (i=0;i<NP;i++) if ((qonstarttest==0 && testingstate[i]==2) || (qonstarttest==1 && testingstate[i]!=0))
  {
  nq[currentloc[i]]++;
  if (infectionstate[i]>=A && infectionstate[i]<=SI) nqASI[currentloc[i]]++;
  }
for (l=0;l<Numloc;l++)
  {
  if (nq[l]!=numquarantinedinloc[l]) printf("Error in numquarantinedinloc (%d,%d) on day %d\n",nq[l],numquarantinedinloc[l],currenttime);
  if (nqASI[l]!=numAtoSIquarantinedinloc[l]) printf("Error in numAtoSIquarantinedinloc (%d,%d) on day %d\n",nqASI[l],numAtoSIquarantinedinloc[l],currenttime);
  }

// loop over people in each loc using linked list
for (l=0;l<Numloc;l++) {nq[l]=nqASI[l]=0;for (i=0;i<Numstates;i++) nsl[i][l]=0;}
for (l=0;l<Numloc;l++)
  {
  j=firstinloc[l];
  while (j>=0)
    {
      for (i=0;i<Numstates;i++) if (infectionstate[j]==i) nsl[i][l]++;
      if ((qonstarttest==0 && testingstate[j]==2) || (qonstarttest==1 && testingstate[j]!=0))
      {
      nq[l]++;
      if (infectionstate[j]>=A && infectionstate[j]<=SI) nqASI[l]++;
      }
    j=nextinloc[j];
    }
  }
for (l=0;l<Numloc;l++)
  {
  for (i=0;i<Numstates;i++) if (nsl[i][l]!=numinstateinloc[i][l]) printf("Error in linked list state %d (%d,%d) in loc %d on day %d\n",i,nsl[i][l],numinstateinloc[i][l],l,currenttime);
  if (nq[l]!=numquarantinedinloc[l]) printf("Error in linked list numquarantinedinloc (%d,%d) in loc %d on day %d\n",nq[l],numquarantinedinloc[l],l,currenttime);
  if (nqASI[l]!=numAtoSIquarantinedinloc[l]) printf("Error in linked list numAtoSIquarantinedinloc (%d,%d) in loc %d on day %d\n",nqASI[l],numAtoSIquarantinedinloc[l],l,currenttime);
  }
*/
return errflag;
}

void checkhealthofpeopleinloclist(int a)
{
int i,j,n;
for (i=0;i<Numloc;i++)
  {
  n=0;
  j=firstinloc[i];
  while (j>=0) {n++;j=nextinloc[j];}
  if (n!=numinloc[i]) {printf("loc list error caused by %d on day %d at loc %d - %d %d\n",a,currenttime,i,n,numinloc[i]);exit(0);}
  }
}

void addtoloc(int loc,int person)
{
numinloc[loc]++;
if (numinloc[loc]==1) {firstinloc[loc]=person;previnloc[person]=-1;nextinloc[person]=-1;}
else {previnloc[person]=-1;nextinloc[person]=firstinloc[loc];previnloc[firstinloc[loc]]=person;firstinloc[loc]=person;}
//checkhealthofpeopleinloclist(0);
}

void removefromloc(int loc,int person)
{
int i,j;
i=nextinloc[person];j=previnloc[person];
if (j>=0) nextinloc[j]=i;
if (i>=0) {previnloc[i]=j;if (j==-1) firstinloc[loc]=i;}
if (i==-1 && j==-1) firstinloc[loc]=-1;
numinloc[loc]--;
//checkhealthofpeopleinloclist(1);
}


void initialize_location_network_ER(void) // sets up location network
{
int i;
// set up links[i][0..2] and linkweight[i] and homeloc[i]
for (i=0;i<NP;i++)
  {
  homeloc[i]=(int)(Numhomes*ran1(&idum));//homeloc[i]=i%Numhomes;
  links[i][0]=homeloc[i];links[i][1]=Numhomes+(int)((Numloc-Numhomes)*ran1(&idum));links[i][2]=i;linkweight[i]=2.0;
  }
}


void initialize_infection(void)
{
int i,j;

currenttime=teststarttime=-1000;teststartmarked=0;

// initialize all to be susceptible except ten random persons, who are marked infectious/asymptomatic/presymptomatic (I=2)
for (i=0;i<Numloc;i++) {numquarantinedinloc[i]=numAtoSIquarantinedinloc[i]=numinloc[i]=0;firstinloc[i]=-1;for (j=0;j<Numstates;j++) numinstateinloc[j][i]=0;}
for (i=0;i<NP;i++) {numinfectedby[i]=testingstate[i]=infectionstate[i]=testtype[i]=0;testtime[i]=infectedtime[i]=-1000.0;currentloc[i]=homeloc[i];numinstateinloc[0][currentloc[i]]++;addtoloc(currentloc[i],i);}
checkhealthofpeopleinloclist(-1);
//printf("all ok so far\n");
for (j=0;j<Numstates;j++) numinstate[j]=numinfectedbystate[j]=0;
numinstate[S]=NP;
for (i=0;i<10;i++)
  {
  do j=(int)(NP*ran1(&idum)); while (infectionstate[j]>0);
  infectionstate[j]=A;
  numinstate[S]--;numinstate[A]++;
  numinstateinloc[S][currentloc[j]]--;numinstateinloc[A][currentloc[j]]++;
  }
checkhealthofnumbers();

currenttime=0;nummove=numtests=numtestsfinished=totaltests=0;
maxinfected=10;maxinfectedtime=currenttime;
numwaitingforresult_symp=numwaitingforresult_other=0;
numtestedpos_symp=numtestedpos_other=0;
}


void infectionmodel_montecarlo_steponeday(double deltat,int tf,FILE *fpout,int iter)
{
int i,j,tempi;
double timewithinday;
int newloc;

double r,psum;
double ninf,ntot,oldninf,oldntot;
int ns[Numstates];
int endtime=currenttime+steptime,middaymove=0;
double midtime=currenttime+0.5*steptime;
int numq,numqh;
int numHCW,numHCWinf,numnonHCW,numnonHCWinf;

fprintf(fpout,"%d ",currenttime);
for (i=0;i<Numstates;i++)
  {
  fprintf(fpout,"%d ",numinstate[i]);
  }
fprintf(fpout,"\n");

if (currenttime>=tf)
  {
  numHCW=numHCWinf=numnonHCW=numnonHCWinf=0;
  for (i=0;i<NP;i++)
      {
      if (links[i][1]>=Numloc-Numhospitals-1) {numHCW++;if (infectionstate[i]!=S) numHCWinf++;}
      else {numnonHCW++;if (infectionstate[i]!=S) numnonHCWinf++;}
      }

  return;
  }

// testing
if (numinstate[R]>=floor(teststartfrac*NP)) // no need for testing if it's too early
  {
  if (teststarttime<0) teststarttime=currenttime; // just for displaying when testing started

  // calc total tests - assuming all are PCR for now
  totaltests=(int)floor(totaltestrate*NP+0.5);
  totalRATleft=(int)floor(totaltests*RATfrac+0.5);
  totalPCRleft=totaltests-totalRATleft;

  // find how many symps can be targetted
  numsymp=0;
  for (i=0;i<NP;i++) if (testingstate[i]==0 && (infectionstate[i]==MI || infectionstate[i]==SI)) {symplist[numsymp]=i;numsymp++;}
  // all symps that can be tested are now in symplist
  if (numsymp<totaltests) // if more tests than symps, then make list of others and pick out totaltests-numsymp of these randomly
    {
    numsymptests=numsymp;
    numother=0;
    for (i=0;i<NP;i++) if (testingstate[i]==0 && (infectionstate[i]==S || infectionstate[i]==A || infectionstate[i]==P || infectionstate[i]==R)) {otherlist[numother]=i;numother++;}
    // now all non-symps who can be tested are in otherlist
    if (numother>totaltests-numsymp)
      {
      for (i=0;i<totaltests-numsymp;i++)
        {
        j=(int)(ran1(&idum)*(numother-i));
        tempi=otherlist[numother-1-i];
        otherlist[numother-1-i]=otherlist[j];
        otherlist[j]=tempi;
        }
      leftovertests=0;numothertests=totaltests-numsymp;
      }
    else {leftovertests=totaltests-numsymp-numother;numothertests=numother;}
    }
  else // otherwise pick out totaltests from the symps (no testing of others will occur)
    {
    for (i=0;i<totaltests;i++)
      {
      j=(int)(ran1(&idum)*(numsymp-i));
      tempi=symplist[numsymp-1-i];
      symplist[numsymp-1-i]=symplist[j];
      symplist[j]=tempi;
      }
    leftovertests=0;numothertests=0;
    numsymptests=totaltests;
    }

  // (i) do tests
  // test numsymptests symps
    for (j=0;j<numsymptests;j++)
      {
      i=symplist[numsymp-1-j];
      if (testingstate[i]!=0) {printf("testing error!\n");exit(0);}
      testtime[i]=currenttime;
      if (totalPCRleft>0) {testtype[i]=0;totalPCRleft--;}
      else {testtype[i]=1;totalRATleft--;}
      numtests++;
      numwaitingforresult_symp++;
      if (((infectionstate[i]==A || infectionstate[i]==P || infectionstate[i]==MI || infectionstate[i]==SI) && ran1(&idum)<testsensitivity[testtype[i]]) || ((infectionstate[i]==S || infectionstate[i]==R) && ran1(&idum)>testspecificity[testtype[i]])) testingstate[i]=1;
      else testingstate[i]=-1;
      if (qonstarttest==1) // quarantine upon starting test if qonstarttest=1
        {
        numinstateinloc[infectionstate[i]][currentloc[i]]--;
        removefromloc(currentloc[i],i);
        currentloc[i]=homeloc[i];numquarantinedinloc[currentloc[i]]++;
        if (infectionstate[i]>=A && infectionstate[i]<=SI) numAtoSIquarantinedinloc[currentloc[i]]++;
        numinstateinloc[infectionstate[i]][currentloc[i]]++;
        addtoloc(currentloc[i],i);
        }
      }

    // test numothertests others
    for (j=0;j<numothertests;j++)
      {
      i=otherlist[numother-1-j];
      if (testingstate[i]!=0) {printf("testing error!\n");exit(0);}
      testtime[i]=currenttime;
      if (ran1(&idum)<totalPCRleft/(totalPCRleft+totalRATleft)) {testtype[i]=0;totalPCRleft--;}
      else {testtype[i]=1;totalRATleft--;}
      numtests++;
      numwaitingforresult_other++;
      if (((infectionstate[i]==A || infectionstate[i]==P || infectionstate[i]==MI || infectionstate[i]==SI) && ran1(&idum)<testsensitivity[testtype[i]]) || ((infectionstate[i]==S || infectionstate[i]==R) && ran1(&idum)>testspecificity[testtype[i]])) testingstate[i]=1;
      else testingstate[i]=-1;
      if (qonstarttest==1) // quarantine upon starting test if qonstarttest=1
        {
        numinstateinloc[infectionstate[i]][currentloc[i]]--;
        removefromloc(currentloc[i],i);
        currentloc[i]=homeloc[i];numquarantinedinloc[currentloc[i]]++;
        if (infectionstate[i]>=A && infectionstate[i]<=SI) numAtoSIquarantinedinloc[currentloc[i]]++;
        numinstateinloc[infectionstate[i]][currentloc[i]]++;
        addtoloc(currentloc[i],i);
        }
      }

  // (ii) declare results
  for (i=0;i<NP;i++) if ((testingstate[i]==1 || testingstate[i]==-1) && currenttime>=testtime[i]+testdelay[testtype[i]])
    {
    testingstate[i]++;//if (infectionstate[i]==R) testingstate[i]=0;
    numtestsfinished++;
    if (infectionstate[i]==MI || infectionstate[i]==SI) {numwaitingforresult_symp--;numtestedpos_symp+=testingstate[i]/2;}
    else {numwaitingforresult_other--;numtestedpos_other+=testingstate[i]/2;}

    testtime[i]=currenttime;
    if (qonstarttest==0 && testingstate[i]==2)
      {
      numinstateinloc[infectionstate[i]][currentloc[i]]--;
      removefromloc(currentloc[i],i);
      currentloc[i]=homeloc[i];numquarantinedinloc[currentloc[i]]++;
      if (infectionstate[i]>=A && infectionstate[i]<=SI) numAtoSIquarantinedinloc[currentloc[i]]++;
      numinstateinloc[infectionstate[i]][currentloc[i]]++;
      addtoloc(currentloc[i],i);
      }
    if (qonstarttest==1 && testingstate[i]==0) {numquarantinedinloc[currentloc[i]]--;if (infectionstate[i]>=A && infectionstate[i]<=SI) numAtoSIquarantinedinloc[currentloc[i]]--;}
    }

  // check to lift quarantine
  for (i=0;i<NP;i++) if (testingstate[i]==2 && currenttime>=testtime[i]+14)// && infectionstate[i]==R)
    {
    testingstate[i]=0;numquarantinedinloc[currentloc[i]]--;if (infectionstate[i]>=A && infectionstate[i]<=SI) numAtoSIquarantinedinloc[currentloc[i]]--;
    if (infectionstate[i]==MI || infectionstate[i]==SI) numtestedpos_symp--;
    else numtestedpos_other--;
    }

    numq=0.0;numqh=0.0;for (i=0;i<Numloc;i++) {numq+=numquarantinedinloc[i];if (numquarantinedinloc[i]>0) numqh++;}


  } // end if (numinstate[R]>teststartfrac*NP)

// move
  for (i=0;i<NE;i++) if (infectionstate[links[i][2]]!=H)
    {
    if (qhouse==1 && numquarantinedinloc[currentloc[links[i][2]]]>0) continue;
    if (qhouse==0 && ((qonstarttest==0 && testingstate[links[i][2]]==2) || (qonstarttest==1 && testingstate[links[i][2]]!=0))) continue;

    // move person along link i if they are not quarantined and not dead
    if (links[i][0]==currentloc[links[i][2]]) newloc=links[i][1];
    else if (links[i][1]==currentloc[links[i][2]]) newloc=links[i][0];
    else continue; // if person is on neither end of this link (say they are in hospital or somewhere else) then don't do anything

    newloc=links[i][1];
    nummove++;
    numinstateinloc[infectionstate[links[i][2]]][currentloc[links[i][2]]]--;
    numinstateinloc[infectionstate[links[i][2]]][newloc]++;
    removefromloc(currentloc[links[i][2]],links[i][2]);
    currentloc[links[i][2]]=newloc;
    addtoloc(currentloc[links[i][2]],links[i][2]);
    }
middaymove=0;


// update infection and move people
timewithinday=currenttime;

while (timewithinday<endtime)
  {
  timewithinday+=deltat;

  // infection transitions
  for (i=0;i<NP;i++)
    {
    switch(infectionstate[i])
      {
      case S:// Susc to Asymptomatic/Presymptomatic
          // Three ways to calc infection rate:
          // 1) Fastest: Use numinstateinloc and numquarantinedinloc and numAtoSIquarantinedinloc
          ninf=numinstateinloc[A][currentloc[i]]+numinstateinloc[P][currentloc[i]]+numinstateinloc[MI][currentloc[i]]+numinstateinloc[SI][currentloc[i]]+hdamping*numinstateinloc[H][currentloc[i]]-(1.0-qdamping)*numAtoSIquarantinedinloc[currentloc[i]];
          ntot=0;for (j=0;j<Numstates;j++) ntot+=numinstateinloc[j][currentloc[i]];
          ntot-=(1-hdamping)*numinstateinloc[H][currentloc[i]];
          oldninf=ninf;oldntot=ntot;
          //ntot=ninf+numinstateinloc[S][currentloc[i]]+numinstateinloc[R][currentloc[i]]; // alternate denominator (now not doing this)
          // 2) 10x slower: Use the linked list which keeps track of who is in each location
          /*for (j=0;j<Numstates;j++) ns[j]=0;
          ninf=0;ntot=numinloc[currentloc[i]];
          j=firstinloc[currentloc[i]];
          while (j>=0)
            {
            ns[infectionstate[j]]++;
            if (infectionstate[j]>=A && infectionstate[j]<=SI) ninf+=((qonstarttest==0 && testingstate[j]==2) || (qonstarttest==1 && testingstate[j]!=0))?qdamping:1.0;
            else if (infectionstate[j]==H) {ninf+=hdamping;ntot-=(1-hdamping);}
            j=nextinloc[j];
            }
          if (abs(oldninf-ninf)>1e-6 || abs(oldntot-ntot)>1e-6)
            {
              printf("%d %d %lf %lf %lf %lf %d\n",currenttime,i,oldninf,ninf,oldntot,ntot,numinloc[currentloc[i]]);
              for (j=0;j<Numstates;j++) printf("%d ",numinstateinloc[j][currentloc[i]]);printf("\n");
              for (j=0;j<Numstates;j++) printf("%d ",ns[j]);printf("\n");
              exit(0);
            }*/
          //ntot=numinloc[currentloc[i]]-(1-hdamping)*numinstateinloc[H][currentloc[i]];// alternate way of calculating denom - uses numinstateinloc, so it's a mixed method
          // 3) 100x slower: loop over all people every time
          /*ninf=ntot=0;
          for (j=0;j<NP;j++) if (currentloc[j]==currentloc[i])
            {
            ntot=ntot+1;
            if (infectionstate[j]>=A && infectionstate[j]<=SI) ninf+=((qonstarttest==0 && testingstate[j]==2) || (qonstarttest==1 && testingstate[j]!=0))?qdamping:1.0;
            else if (infectionstate[j]==H) {ninf+=hdamping;ntot-=(1-hdamping);}
            }*/
          r=ran1(&idum);
          if (r<ninf*exitrate[S]*deltat/ntot)
            {
            if (ninf*exitrate[S]*deltat/ntot>1) printf("prob exceeded %d %lf\n",currenttime,(ninf*exitrate[S]*deltat/ntot));
            if (ran1(&idum)<pSP) infectionstate[i]=P;
            else infectionstate[i]=A;
            infectedtime[i]=currenttime;
            numinstateinloc[S][currentloc[i]]--;numinstateinloc[infectionstate[i]][currentloc[i]]++;
            numinstate[S]--;numinstate[infectionstate[i]]++;

            /*psum=0;for (j=0;j<NP;j++) if (currentloc[j]==currentloc[i] && j!=i)
              {
              if (infectionstate[j]==A || infectionstate[j]==P || infectionstate[j]==MI || infectionstate[j]==SI)
                {
                psum+=((qonstarttest==0 && testingstate[i]==2) || (qonstarttest==1 && testingstate[i]!=0))?qdamping:1.0;
                }
              else if (infectionstate[j]==H) psum+=hdamping;
              if (psum>r*ninf) break;
              }
            numinfectedby[j]++;
            numinfectedbystate[infectionstate[j]]++;*/
            if ((qonstarttest==0 && testingstate[i]==2) || (qonstarttest==1 && testingstate[i]!=0)) numAtoSIquarantinedinloc[currentloc[i]]++;
            }
          break;
      case A:// Asymptomatic to Recovered
          if (ran1(&idum)<exitrate[A]*deltat)
            {
            infectionstate[i]=R;
            numinstateinloc[A][currentloc[i]]--;numinstateinloc[R][currentloc[i]]++;
            numinstate[A]--;numinstate[R]++;
            //if (testingstate[i]==2) {numquarantinedinloc[currentloc[i]]--;testingstate[i]=0;}
            if ((qonstarttest==0 && testingstate[i]==2) || (qonstarttest==1 && testingstate[i]!=0)) numAtoSIquarantinedinloc[currentloc[i]]--;
            }
          break;
      case P:// Presymptomatic to MI or SI
          if (ran1(&idum)<exitrate[P]*deltat)
            {
            if (ran1(&idum)<pPSI) infectionstate[i]=SI;
            else infectionstate[i]=MI;
            numinstateinloc[P][currentloc[i]]--;numinstateinloc[infectionstate[i]][currentloc[i]]++;
            numinstate[P]--;numinstate[infectionstate[i]]++;
            }
          break;
      case MI:// Mildly Symptomatic to Recovered
          if (ran1(&idum)<exitrate[MI]*deltat)
            {
            infectionstate[i]=R;
            numinstateinloc[MI][currentloc[i]]--;numinstateinloc[R][currentloc[i]]++;
            numinstate[MI]--;numinstate[R]++;
            //if (testingstate[i]==2) {numquarantinedinloc[currentloc[i]]--;testingstate[i]=0;}
            if ((qonstarttest==0 && testingstate[i]==2) || (qonstarttest==1 && testingstate[i]!=0)) numAtoSIquarantinedinloc[currentloc[i]]--;
            }
          break;
      case SI: // SI to recovered or hospitalised
          if (ran1(&idum)<exitrate[SI]*deltat)
            {
            if (ran1(&idum)<pSIH) infectionstate[i]=H;
            else infectionstate[i]=R;
            numinstateinloc[SI][currentloc[i]]--;
            numinstate[SI]--;
            if (infectionstate[i]==R)
              {
              numinstate[R]++;numinstateinloc[R][currentloc[i]]++;
              //if (testingstate[i]==2) numquarantinedinloc[currentloc[i]]--;
              if ((qonstarttest==0 && testingstate[i]==2) || (qonstarttest==1 && testingstate[i]!=0)) numAtoSIquarantinedinloc[currentloc[i]]--;
              }
            else
              {// move to hospital and update numinstateinloc there
              if ((qonstarttest==0 && testingstate[i]==2) || (qonstarttest==1 && testingstate[i]!=0)) {testingstate[i]=0;numquarantinedinloc[currentloc[i]]--;numAtoSIquarantinedinloc[currentloc[i]]--;}
              removefromloc(currentloc[i],i);
              currentloc[i]=(int)(ran1(&idum)*Numhospitals)+Numloc-Numhospitals;
              addtoloc(currentloc[i],i);
              numinstate[H]++;numinstateinloc[H][currentloc[i]]++;
              //if (checkhealthofnumbers()==1) {printf("iter %d currenttime %d transition of %d sent to %d home is %d work is %d; %lf\n",iter,currenttime,i,currentloc[i],links[i][0],links[i][1],r);exit(0);}
              }
            }
          break;
      case H:// Hospitalized to Recovered
          if (ran1(&idum)<exitrate[H]*deltat)
            {
            infectionstate[i]=R;
            numinstateinloc[H][currentloc[i]]--;
            numinstate[H]--;numinstate[R]++;
            // move back home and update numinstateinloc there
            removefromloc(currentloc[i],i);
            currentloc[i]=homeloc[i];
            addtoloc(currentloc[i],i);
            numinstateinloc[R][currentloc[i]]++;
            }
          break;
      }
    }


  // move
  if (timewithinday>midtime && middaymove==0)
    {
    for (i=0;i<NE;i++) if (infectionstate[links[i][2]]!=H)
      {
      if (qhouse==1 && numquarantinedinloc[currentloc[links[i][2]]]>0) continue;
      if (qhouse==0 && ((qonstarttest==0 && testingstate[links[i][2]]==2) || (qonstarttest==1 && testingstate[i]!=0))) continue;

      // move person along link i if they are not quarantined and not dead
      if (links[i][0]==currentloc[links[i][2]]) newloc=links[i][1];
      else if (links[i][1]==currentloc[links[i][2]]) newloc=links[i][0];
      else continue; // if person is on neither end of this link (say they are in hospital or somewhere else) then don't do anything

      newloc=links[i][0];
      nummove++;
      numinstateinloc[infectionstate[links[i][2]]][currentloc[links[i][2]]]--;
      numinstateinloc[infectionstate[links[i][2]]][newloc]++;
      removefromloc(currentloc[links[i][2]],links[i][2]);
      currentloc[links[i][2]]=newloc;
      addtoloc(currentloc[links[i][2]],links[i][2]);
      }
    middaymove=1;
    }

  /*for (i=0;i<NE;i++) if (ran1(&idum)<linkweight[i]*deltat && infectionstate[links[i][2]]!=H && ((qhouse==0 && testingstate[links[i][2]]!=2) || (qhouse==1 && numquarantinedinloc[currentloc[links[i][2]]]==0)))
   */
  } // end at endtime

checkhealthofpeopleinloclist(100);
currenttime++;
//printf("%d\n",currenttime);
checkhealthofnumbers();
}

void infectionmodel_montecarlo_onerun(int tf,int iter)
{
int i;
FILE *fpout;
char str[200];


// init testing protocol
qonstarttest=0;
qhouse=0;

// initialize rates
hdamping=0.1;
qdamping=0.1;

pSP=0.5;pPSI=0.15;pSIH=0.8;
Irecrate=0.1;
exitrate[S]=0.25;
exitrate[A]=1.0/7.0;
exitrate[P]=0.5;
exitrate[MI]=exitrate[H]=Irecrate;
exitrate[SI]=Irecrate/(1-pSIH);

initialize_location_network_ER();
initialize_infection();
//"Targeted-Testing_BeginAt_%g_DTR_%g_RAT_%g_%g_PCR_%g_%g_%s-%i.txt"
sprintf(str,"../tmp/SK_BeginAt_%g_DTR_%g_RAT_%g_%g_PCR_%g_%g_%lf%lf-%i.txt",100*teststartfrac,100*totaltestrate,testsensitivity[1],RATfrac,testsensitivity[0],1-RATfrac,ran1(&idum),ran1(&idum),iter);
//sprintf(str,"beta_0.5_DTR_100_qd_0.5.txt");
fpout=(FILE *)fopen(str,"wt");
currenttime=0;
while(currenttime<tf) infectionmodel_montecarlo_steponeday(deltat,tf,fpout,iter);
fclose(fpout);

}

int main(void)
{
srand(time(0));
idum=-rand();

// double p[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
// int p_len  = sizeof(p)/sizeof(p[0]);
//
// double s[] = {0.5, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95};
// int s_len  = sizeof(s)/sizeof(s[0]);
//
// double dtr[]={0.001, 0.005};
// int dtr_len = sizeof(dtr)/sizeof(dtr[0]);

int mc_runs = 50;


//double p[] = {0,2,4,6,8,10};
//int p_len  = sizeof(p)/sizeof(p[0]);
//
//double s[] = {1, 5, 10, 20, 50};
//int s_len  = sizeof(s)/sizeof(s[0]);
//
//double dtr[]={0.1, 0.5, 1, 10, 20, 40, 50};
//int dtr_len = sizeof(dtr)/sizeof(dtr[0]);
//
//int mc_runs =20;
//



for(int i=0;i<mc_runs;i++){
  // for(int j=0;j<s_len;j++){
  //
  //   // Tpars[0][0] = s[j];
  //   testsensitivity[2] = s[j];
  //
  //
  //   for(int k=0; k<p_len;k++){
  //
  //     // Tpars[0][3] = p[k];
  //     // Tpars[1][3] = 1-p[k];
  //
  //     RATfrac = p[k];
  //
  //     for(int l =0; l<dtr_len;l++){
  //       totaltestrate = dtr[l];
  //       //  idum=-i; // set random number seed for each run
  //         //printf("iter %d\n",i);
    infectionmodel_montecarlo_onerun(200,i);
//     }
// }
// }
}

}
