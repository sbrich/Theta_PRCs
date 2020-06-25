#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define ton 0.
#define toff 20000
#define dt .01
#define tend 20000.
#define repnum 1

#define numcells 1

//#define cm 90.0
//#define vr -60.6
//#define vt -43.1
//#define vpeak 2.5
//#define aizh 0.1
//#define bizh -0.1
//#define cizh -67.0
//#define dizh 0.1
//#define klow 1.7
//#define khigh 14.0

#define Esyn -75.0
#define alphainv 0.27
#define betainv 3.0

void matrixmult(double** A, int rowsA, int colsA, double** B, int rowsB, int colsB, double** answer)
{
    int i, j, k;
    for (i=0; i<rowsA; i++)
    {
        for (j=0; j<colsB; j++)
        {
            answer[i][j] = 0;
            for (k=0; k<colsA; k++)
            {
                answer[i][j]=answer[i][j]+(A[i][k])*(B[k][j]);
            }
        }
    }
}


void ScalarMultMatrix(double** A, int rowsA, int colsA, double scalar, double** answer)
{
    int i, j;
    for (i=0; i<rowsA; i++)
    {
        for (j=0; j<colsA; j++)
        {
            answer[i][j]=A[i][j]*scalar;
        }
    }
}


double** MatrixAdd(double** A, int rowsA, int colsA, double** B, int rowsB, int colsB)
{
    double** answer;
    {
        answer=(double**) malloc(rowsA*sizeof(double*));
        int a;
        for (a=0;a<rowsA; a++)
        {
            answer[a]=(double*) malloc(colsB*sizeof(double));
        }
        
        int i, j;
        for (i=0; i<rowsA; i++)
        {
            for (j=0; j<colsB; j++)
            {
                answer[i][j]=A[i][j]+B[i][j];
            }
        }
    }
    return answer;
}

double dvdt (double V, double u, double Iapp, double vt, double klow, double khigh, double vr, double cm)
{
    double dvdt;
    double k;
    if (V<=vt)
    {
        k=klow;
    }
    else if (V>vt)
    {
        k=khigh;
    }
    dvdt=(k*(V-vr)*(V-vt)-u+Iapp)/cm;
    return dvdt;
}

double dudt (double V, double u, double Iapp, double aizh, double bizh, double vr)
{
    double dudt;
    dudt=aizh*(bizh*(V-vr)-u);
    return dudt;
}

double randn (double mu, double sigma)
{
    double U1, U2, W, mult;
    static double X1, X2;
    static int call = 0;
    
    if (call == 1)
    {
        call = !call;
        return (mu + sigma * (double) X2);
    }
    
    do
    {
        U1 = -1 + ((double) rand () / RAND_MAX) * 2;
        U2 = -1 + ((double) rand () / RAND_MAX) * 2;
        W = pow (U1, 2) + pow (U2, 2);
    }
    while (W >= 1 || W == 0);
    
    mult = sqrt ((-2 * log (W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;
    
    call = !call;
    
    return (mu + sigma * (double) X1);
}

double ds0dt (double s0, double sinf, double taus)
{
    double ds0dt;
    ds0dt=-(s0-sinf)/taus;
    return ds0dt;
}

double ds1dt (double s1, double beta)
{
    double ds1dt;
    ds1dt=-s1*beta;
    return ds1dt;
}







int main (int argc, char *argv[])
{
    //INPUT PARAMETERS
    
    //USE gsyn=-500, perturblength=1, Iapp 20 or 30, precision=100
    double gsyn;
    double Iapp;
    int precision;
    int number;
    int state;
    double perturblength;
    
    gsyn=atof(argv[1]);
    Iapp=atof(argv[2]);
    precision=atoi(argv[3]);
//    number=atof(argv[4]);
    //    state=atoi(argv[5]);
    perturblength=atof(argv[4]);
    
    double cm;
    double vr;
    double vt;
    double vpeak;
    double aizh;
    double bizh;
    double cizh;
    double dizh;
    double klow;
    double khigh;
    
    int a;
    int i;
    int j;
    
    double** params;
    params=(double**) malloc((25)*sizeof(double*));
    for (a=0; a<25; a++)
    {
        params[a]=(double*) malloc(4*sizeof(double));
    }
    
    FILE *fp = fopen("MMH_model_params.csv", "r");
    char buf[1024];
    char *eptr;
    int row_count = 0;
    int field_count = 0;
    
    while (fgets(buf, 1024, fp)) {
        field_count = 0;
        //        row_count++;
        
        char *field = strtok(buf, ",");
        while (field){
            params[row_count][field_count]=strtod(field, &eptr);
            field_count++;
            field = strtok(NULL, ",");
        }
        
        row_count++;
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    //    printf("%f, %f, %f, %d\n", p, gsyn, probii, number);
    int steps;
    steps=tend/dt;
    
    double** now;
    now=(double**) malloc((3)*sizeof(double*));
    for (a=0; a<3; a++)
    {
        now[a]=(double*) malloc(numcells*sizeof(double));
    }
    
    double*** answer;
    answer=(double***) malloc((steps+1)*sizeof(double**));
    for (a=0; a<steps+1; a++)
    {
        answer[a]=(double**) malloc((3)*sizeof(double*));
        int b;
        for (b=0; b<3; b++)
        {
            answer[a][b]=(double*) malloc(numcells*sizeof(double));
        }
    }
    
    double** k1;
    k1=(double**) malloc((2)*sizeof(double*));
    for (a=0; a<2; a++)
    {
        k1[a]=(double*) malloc(numcells*sizeof(double));
    }
    
    double** k2;
    k2=(double**) malloc((2)*sizeof(double*));
    for (a=0; a<2; a++)
    {
        k2[a]=(double*) malloc(numcells*sizeof(double));
    }
    
    double** k3;
    k3=(double**) malloc((2)*sizeof(double*));
    for (a=0; a<2; a++)
    {
        k3[a]=(double*) malloc(numcells*sizeof(double));
    }
    
    double** k4;
    k4=(double**) malloc((2)*sizeof(double*));
    for (a=0; a<2; a++)
    {
        k4[a]=(double*) malloc(numcells*sizeof(double));
    }
    
    double* times;
    times=(double*) malloc((steps+1)*sizeof(double));
    
    double** Isyn;
    Isyn=(double**) malloc((numcells)*sizeof(double*));
    for (a=0; a<numcells; a++)
    {
        Isyn[a]=(double*) malloc(1*sizeof(double));
    }
    
    double** Ssyn;
    Ssyn=(double**) malloc((numcells)*sizeof(double*));
    for (a=0; a<numcells; a++)
    {
        Ssyn[a]=(double*) malloc(1*sizeof(double));
    }
    
    
    double** appcurr;
    appcurr=(double**) malloc((numcells)*sizeof(double*));
    for (a=0; a<numcells; a++)
    {
        appcurr[a]=(double*) malloc(1*sizeof(double));
    }
    
    double** spiketimes;
    spiketimes=(double**) malloc((steps+1)*sizeof(double*));
    int s;
    for (s=0;s<(steps+1); s++)
    {
        spiketimes[s]=(double*) malloc(numcells*sizeof(double));
    }
    
    int* spikecount;
    spikecount=(int*) malloc(numcells*sizeof(int));
    
    double* Inow;
    Inow=(double*) malloc(numcells*sizeof(double));
    
    double* Itotal;
    Itotal=(double*) malloc(numcells*sizeof(double));
    
    double* I;
    I=(double*) malloc(numcells*sizeof(double));
    
    int rep;
    
    struct timeval time;
    gettimeofday(&time,NULL);
    
    // microsecond has 1 000 000
    // Assuming you did not need quite that accuracy
    // Also do not assume the system clock has that accuracy.
    srand((time.tv_sec * 1000) + (time.tv_usec / 1000)+ number);
    
    // The trouble here is that the seed will repeat every
    // 24 days or so.
    
    // If you use 100 (rather than 1000) the seed repeats every 248 days.
    
    // Do not make the MISTAKE of using just the tv_usec
    // This will mean your seed repeats every second.
    
    //    printf("Opening output file\n");
    
    for (int paramnum=0; paramnum<25;paramnum++){
        
        char f1[100], f2[100], f3[100], f4[100], f5[100];
        
        int Iapptemp=Iapp;
//        sprintf(f1, "InhibitoryNetwork_%d_TrackVariables_MMH.csv", paramnum);
//        sprintf(f2, "InhibitoryNetwork_%d_SpikeTimes_MMH.csv", paramnum);
//        sprintf(f3, "InhibitoryNetwork_%d_InputCurrents_MMH.csv", paramnum);
        sprintf(f4, "InhibitoryNetwork_PRC_MMH_%d_%d.csv", paramnum, Iapptemp);
        sprintf(f5, "InhibitoryNetwork_Parameters_MMH_%d_%d.csv", paramnum, Iapptemp);
        
//        FILE *Output;
//        Output = fopen(f1, "wt");
//
//        FILE *Output2;
//        Output2 = fopen(f2, "wt");
//
//        FILE *Output3;
//        Output3 = fopen(f3, "wt");
        
        FILE *Output4;
        Output4 = fopen(f4, "wt");
        
        FILE *Output5;
        Output5 = fopen(f5, "wt");
        
//        printf("%s\n", f1);
//        printf("%s\n", f2);
//        printf("%s\n", f3);
        printf("%s\n", f4);
        printf("%s\n", f5);
        
        cm=115;
        vr=-61.8;
        vt=-57.0;
        vpeak=22.6;
        aizh=params[paramnum][0];
        bizh=params[paramnum][1];
        cizh=-65.8;
        dizh=params[paramnum][2];
        klow=params[paramnum][3];
        khigh=3.3;

        
        double** connectivity;
        connectivity=(double**) malloc((numcells)*sizeof(double*));
        int w;
        for (w=0;w<numcells; w++)
        {
            connectivity[w]=(double*) malloc(numcells*sizeof(double));
        }
        
        int b;
        int k;
        
        double TotalAverageCurrent;
        double TotalAverageVoltage;
        
        double Iappavg;
        
        
        
        
        
        // FIND ME
        double alpha;
        double beta;
        alpha=1/alphainv;
        beta=1/betainv;
        double sinf;
        double taus;
        sinf=alpha/(alpha+beta);
        taus=1/(alpha+beta);
        double temptime1;
        int temptime2;
        
        // INITIAL CONDITIONS
        for(j=0;j<numcells;j++)
        {
            now[0][j]=-60.0;
            now[1][j]=0;
            now[2][j]=0;
        }
        
        // HETEROGENEITY
        for(j=0; j<numcells; j++)
        {
            I[j]=randn(Iapp, 0);
        }
        
        for(j=0;j<numcells;j++)
        {
            Ssyn[j][0]=0;
        }
        
        double period;
        int perturbtime;
        int perturbtimeend;
        int p;
        double phaseshift;
        
        for (p=0; p<=precision; p++)
        {
            for (j=0; j<numcells; j++)
            {
                spikecount[j]=0;
            }
            
            for (s=0;s<(steps+1); s++)
            {
                for (j=0;j<numcells;j++)
                {
                    spiketimes[s][j]=0;
                }
            }
            perturbtime=-1;
            perturbtimeend=-1;
            
            // INITIAL CONDITIONS
            for(j=0;j<numcells;j++)
            {
                now[0][j]=-60.0;
                now[1][j]=0;
                now[2][j]=0;
            }
            
            //    printf("Entering time loop\n");
            for(i=0; i<=steps; i++)
            {
                //        printf("Move now into answer\n");
                for(j=0; j<numcells; j++)
                {
                    answer[i][0][j]=now[0][j];
                    //        printf("a\n");
                    answer[i][1][j]=now[1][j];
                    //        printf("b\n");
                    answer[i][2][j]=now[2][j];
                }
                times[i]=i*dt;
                
                //        printf("Set Inow\n");
                Inow[0]=Iapp;
                
                if (i>=perturbtime && i<=perturbtimeend && spikecount[0]==10)
                {
                    Inow[0]=Inow[0]+gsyn;
                }
                
                // FORWARD EULER
                for (j=0; j<numcells; j++)
                {
                    now[0][j]=now[0][j]+dvdt(now[0][j],now[1][j],Inow[j], vt, klow, khigh, vr, cm)*dt;
                    now[1][j]=now[1][j]+dudt(now[0][j],now[1][j],Inow[j], aizh, bizh, vr)*dt;
                    if (now[0][j]>=vpeak)
                    {
                        now[0][j]=cizh;
                        now[1][j]=now[1][j]+dizh;
                        spiketimes[spikecount[j]][j]=(i)*dt;
                        spikecount[j]=spikecount[j]+1;
                        if (spikecount[0]==10)
                        {
                            period=spiketimes[9][0]-spiketimes[8][0];
                            perturbtime=i+p*((period/dt)/precision);
                            perturbtimeend=perturbtime+perturblength/dt;
                            printf("perturbtime=%3d\n", perturbtime);
                        }
                        if (spikecount[0]==11)
                        {
                            phaseshift=(period-(spiketimes[10][0]-spiketimes[9][0]))/period;
                            fprintf(Output4, "%3d, %3f, %3f\n", p, phaseshift, period);
                            printf("phaseshift=%3f\n", phaseshift);
                        }
                    }
                    
                    //            if (spikecount[j]>0)
                    //            {
                    //                // SYNAPSE TURN ON TIME
                    //                //                        if (spiketimes[spikecount[j]-1][j]>100)
                    //                if (spiketimes[spikecount[j]-1][j]>0)
                    //                {
                    //                    temptime1=times[i]-spiketimes[spikecount[j]-1][j];
                    //                    if (temptime1<1)
                    //                    {
                    //                        now[2][j]=now[2][j]+ds0dt(now[2][j],sinf,taus)*dt;
                    //                    }
                    //                    else
                    //                    {
                    //                        now[2][j]=now[2][j]+ds1dt(now[2][j], beta)*dt;
                    //                    }
                    //                }
                    //            }
                }
                
                
                //            if ((i % 2000) == 0)
                //            {
                //                printf("%3f\n", i*dt);
                //            }
                
                //        fprintf(Output, "%3f,", times[i]);
                //        //                fprintf(Output, "%3f,", TotalAverageCurrent);
                //        //                fprintf(Output, "%3f,", TotalAverageVoltage);
                //        for (j=0; j<10; j++)
                //        {
                //            fprintf(Output, "%3f, %3f, %3f,", now[0][j], appcurr[j][0], Itotal[j]);
                //        }
                //        fprintf(Output, "\n");
                
                if (i==steps)
                {
                    fprintf(Output5, "%3f, %3f, %3d, %3f, %3d", gsyn, Iapp, precision, perturblength, paramnum);
                    fprintf(Output5, "\n");
                    
//                    for(j=0; j<numcells; j++)
//                    {
//                        fprintf(Output3, "%3f\n", I[j]);
//                    }
//
//                    //    printf("Print spiketimes data \n");
//                    for(j=0; j<numcells; j++)
//                    {
//                        if (spikecount[j]>0)
//                        {
//                            int k;
//                            for (k=0; k<spikecount[j]; k++)
//                            {
//                                fprintf(Output2, "%3f,", spiketimes[k][j]);
//                            }
//                        }
//                        if (spikecount[j]==0)
//                        {
//                            fprintf(Output2, "%3f,", 0.0);
//                            //            fprintf(Output2, "\n");
//                        }
//                        fprintf(Output2, "\n");
//                    }
//                    fprintf(Output2, "\n");
//                    //    printf("End printing spiketimes data\n");
                    
                    printf("DONE %d\n", p);
                }
            }
        }
        
        
        
        //    printf("Closing output file\n");
//        fclose(Output);
//        fclose(Output2);
//        fclose(Output3);
        fclose(Output4);
        fclose(Output5);
    }
    
    return 0;
}
