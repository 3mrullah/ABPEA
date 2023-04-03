/*
Coded by:   Emrullah Sonuç
            http://3mrullah.com

E-mail:     esonuc@karabuk.edu.tr
            emrullah.sonuc@nottingham.ac.uk

Paper:      Sonuç, E., & Özcan, E. (2023). 
			An adaptive parallel evolutionary algorithm for solving the uncapacitated facility location problem. 
			Expert Systems with Applications, 119956.

DOI:        https://doi.org/10.1016/j.eswa.2023.119956


*/
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>
#include<float.h>
#include <omp.h>

#define SHAKE 25
#define phiMAX 0.9
#define phiMIN 0.5
#define W 5
//#define PMIN_ 0.1
//#define RT_ 1 // Reward Type 1: Instant, 2:Avg, 3:Ext
#define ALPHA 0.9
#define MAX_ITER 80000

#define THREAD_NUM  512           // Number Of Threads 

double PMIN = 0.0;
int RT = 0;

static inline double closed_interval_rand(double x0, double x1)
{
    return x0 + (x1 - x0) * rand() / ((double) RAND_MAX);	
}

void CalcDisSim(int loc, bool xi[loc], bool xk[loc], bool ind_x[][loc], int tid, int follower, int max_iter, int iter)
{
    // N modifiye edilecek, max_iter ve iter kontrol edilecek  
    int i, j, jj, k;
    int M11=0, M01=0, M10=0;
    double sim, phi, A, x;
    int n1=0, n0=0;
    double z = DBL_MAX, zt;
    int m11v = 0, m10v = 0;

    int rand1[loc];
    int rand0[loc];
    memset( rand1, -1, loc*sizeof(int) );
    memset( rand0, -1, loc*sizeof(int) );
	
	
    for(i=0; i<loc; i++)
    {
        if (xi[i] == 1 && xk[i] == 1)
            M11++;
        else if (xi[i] == 0 && xk[i] == 1)
            M01++;
        else if (xi[i] == 1 && xk[i] == 0)
            M10++;

        if (xi[i] == 1)
            n1++;
        else
            n0++;

        ind_x[tid][i] = xi[i];
    }
	
	//printf("M11 = %d, M01 = %d, M10 = %d, n1 = %d, n0 = %d\n", M11, M01, M10, n1, n0);
	
    /*printf("\n");
	for(i=0; i<loc; i++)
        printf("%d ", ind_x[tid][i]);
	printf(" (%d) ", tid);*/
	
	//getchar();
	
    sim = M01 + M10 + 2*M11;
	
	
    if(sim != 0 && M11!=0 && (M01!=0 || M10!=0))
    {

    sim = 2*M11 / sim;

    phi = phiMAX - ( ((phiMAX - phiMIN) / max_iter) * iter );
    A = phi * (1-sim);

    for(i=0; i<n1; i++)
        for(j=0; j<n0; j++)
        {
            x = i * 1.0 / (n1+j)* 1.0;
            zt = fabs(1 - x - A);

            if (zt < z)
            {
                z = zt;
                m11v = i;
                m10v = j;
            }
			//printf("a");
        }


    // m11! & m10!
    j = 0;    
    jj = 0;    
	for(i=0; i<loc; i++)
	{
		if (ind_x[tid][i] == 1)
		{
			rand1[j] = i;
			j++;
			//printf("b1");
		}
		else
		{
			rand0[jj] = i;
			jj++;
		}
		//printf("b2");
	}
	//printf("\nj = %d , jj = %d", j, jj);
   
	//printf("\nm11v = %d , nm10v = %d", m11v, m10v);
	//printf("\n");
	/*for(i=0; i<loc; i++)
		printf(" %d ",ind_x[tid][i]);*/
	
	// m11!
    k = n1 - m11v;
    while (k>0)
    {
        i=ceil(j*closed_interval_rand(0,1))-1;
		if(i==-1) // cok threadli de i sacmaliyor -1 donuyor o hatayi duzeltmek icin
			i=0;
        if (ind_x[tid][rand1[i]] == 1)
        {
            ind_x[tid][rand1[i]] = 0;
            k = k -1;
			
			rand1[i] = rand1[j-1];
			j = j -1;
        }
		
    }	
	
	/*printf("\nm11v = %d , nm10v = %d", m11v, m10v);
	printf("\nR");
	for(i=0; i<jj; i++)
		printf(" %d ",rand0[i]);*/
	

	//printf("\nm10v = %d , jj = %d", m10v, jj);
    k = m10v;	
    while (k>0)
    {
        i=ceil(jj*closed_interval_rand(0,1))-1;
		if(i==-1) // cok threadli de i sacmaliyor -1 donuyor o hatayi duzeltmek icin
			i=0;
		//printf("jj %d k %d",jj, k);
		//printf("\nr[%d]%d", i, rand0[i]);
        if (ind_x[tid][rand0[i]] == 0)
        {
            ind_x[tid][rand0[i]] = 1;
            k = k -1;
			
			rand0[i] = rand0[jj-1];
			jj = jj -1;
			
        }
		
		//printf("\nk = %d (%d -> %d)",k,i, jj);
		//printf("\n\n");
		//for(i=0; i<loc; i++)
		//	printf(" %d ",ind_x[tid][i]);
		//printf("e");
    }
	
    // m10!

     // return new individual
	 
	/*printf("\n");
	for(i=0; i<loc; i++)
        printf("%d ", ind_x[tid][i]);
	printf(" (_%d) ", tid);*/
	//getchar();
	 
	}
}


void shuffle(int *array, size_t n)
{
    size_t i,j;
    int t;
	
    if (n > 1)
    {
        for (i = 0; i < n - 1; i++)
        {
          j = i + rand() / (RAND_MAX / (n - i) + 1);
          t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}

double UFLP(int loc, int cus, double customer[][loc], double location[loc], bool *per)
{
    double cost = 0;
    double min_cost = 0;
    int i,j;

    for(i=0; i<loc; i++)
        if (per[i] == 1)
            cost = cost + location[i];
	
    for(i=0; i<cus; i++)
    {
        min_cost = 9999999999;
        for(j=0; j<loc; j++)
        {
            if (per[j] == 1 && customer[i][j]<min_cost)
                min_cost = customer[i][j];
        }
        cost = cost + min_cost;
    }
	
    return cost;
}



int main(int argc, char**argv)
{
    clock_t begin, end;
    double time_spent;
    begin = clock();
    srand(time(NULL));
    int cus, loc;
    int i,j,iter,follower,best_agent, chck;

    char infile[80],indirectory[]="ORLIB-uncap\\\\";
    if (argc <= 2){
	  printf("Usage:\n\t%s input file\n", argv[0]);
	  exit(1);
    }
    strcpy(infile, argv[1]);
    strcat(indirectory,infile);
    FILE *fp = fopen(indirectory, "r");
    fscanf(fp, "%d", &loc);
    fscanf(fp, "%d", &cus);
    double location[loc];
    double customer[cus][loc];
    double temp;
    int tt;
    for(i=0;i<loc;i++)
    {
        fscanf(fp, "%s", &temp);
        fscanf(fp, "%lf", &location[i]);
    }
	
    for(i=0;i<cus;i++)
    {
        fscanf(fp, "%d",&tt);
        for(j=0;j<loc;j++)
        {
            fscanf(fp, "%lf", &customer[i][j]);
        }
    }
	PMIN = atof(argv[2]);
	RT = atoi(argv[3]);
	
	omp_set_num_threads(THREAD_NUM);
	
    double obj_agents[THREAD_NUM] = {0};
    double obj_agents_new[THREAD_NUM] = {0};
    double obj_agents_best[THREAD_NUM] = {0};
    double reward[THREAD_NUM][W] = {0};
    int opselected[THREAD_NUM] = {0};
	
    bool ind_agents[THREAD_NUM][loc];
    bool ind_agents_new[THREAD_NUM][loc];
    bool ind_agents_best[THREAD_NUM][loc];
	
	bool ind_dissim_xi[loc];
	bool ind_dissim_xk[loc];
	
	
	int tid;
	int best_agent_;
	int control=0;
	int shrink=0;
	
    memset( ind_agents, 0, THREAD_NUM*loc*sizeof(bool) );
    memset( ind_agents_new, 0, THREAD_NUM*loc*sizeof(bool) );
    memset( ind_agents_best, 0, THREAD_NUM*loc*sizeof(bool) );
	
    memset( ind_dissim_xi, 0, loc*sizeof(bool) );
    memset( ind_dissim_xk, 0, loc*sizeof(bool) );

    double global_best = 9999999999;
    double global_best_ = 9999999999;
	
    closed_interval_rand(0,1);

	

		
	// INITIALIZATION BEGINS
	#pragma omp parallel shared(ind_agents, ind_agents_best, loc, cus, customer, location, best_agent, global_best, reward) private(i, j, tid)
    {
        tid = omp_get_thread_num();
		//printf("%d ", tid);
		srand(time(NULL) * tid * tid);
		
		if(THREAD_NUM>loc)
		{	
		
			if(tid<loc)
			{
				ind_agents[tid][tid%loc] = 1;
				ind_agents_new[tid][tid%loc] = 1;
				ind_agents_best[tid][tid%loc] = 1;		
				
			}
			else
			{	
				for (j=0;j< loc; j++)
				{				
					if ( closed_interval_rand(0,1) < 0.5 )
					{
						ind_agents[tid][j] = 0;
						ind_agents_new[tid][j] = 0;
						ind_agents_best[tid][j] = 0;
					}
					else
					{
						ind_agents[tid][j] = 1;
						ind_agents_new[tid][j] = 1;
						ind_agents_best[tid][j] = 1;
					}	
				}
				
			}
		}
		else
		{
			for (j=0;j< loc; j++)
				{				
					if ( closed_interval_rand(0,1) < 0.5 )
					{
						ind_agents[tid][j] = 0;
						ind_agents_new[tid][j] = 0;
						ind_agents_best[tid][j] = 0;
					}
					else
					{
						ind_agents[tid][j] = 1;
						ind_agents_new[tid][j] = 1;
						ind_agents_best[tid][j] = 1;
					}	
				}
		}
		
 		double q = 99999999;
		int t = 0;
		for (j=0;j< loc; j++)	
			if (location[j]<q)
			{
				q = location[j];
				t = j;
			}		
		
			ind_agents[tid][t] = 1;
			ind_agents_new[tid][t] = 1;
			ind_agents_best[tid][t] = 1; 
				
		
		obj_agents_best[tid] = UFLP(loc, cus, customer, location, ind_agents[tid]);
		obj_agents[tid] = obj_agents_best[tid];
		obj_agents_new[tid] = obj_agents_best[tid];
		for(i=0; i<W; i++)
			reward[tid][i] = 0;
		
		#pragma omp barrier
		if (obj_agents_best[tid] < global_best)
		{
			global_best = obj_agents_best[tid];
			best_agent = tid;
		}
    }
	//printf("\nBest: \t %.3f (%d)\n", global_best, best_agent);	
	// INITIALIZATION ENDS
	
	double c[3] = {0, 0, 0};
	double pr[3] = {0.33, 0.33, 0.33};
	double q = 0;
	int sr[3] = {0, 0, 0};
	int opcount[3] = {0, 0, 0};
	// op1_c0 - bincsa
	// op2_c1 - dicesim
	// op3_c2 - elitemut
	for (iter=0; iter<MAX_ITER/THREAD_NUM; iter++)
    {
				
		#pragma omp parallel shared(ind_agents, ind_agents_new, ind_agents_best, obj_agents, obj_agents_new, obj_agents_best, loc, cus, customer, location, global_best, global_best_, best_agent_, best_agent, reward, opselected, opcount, c, pr, sr, shrink, control) private(i, j, follower, tid, q)
		{
			tid = omp_get_thread_num();
			follower=ceil(THREAD_NUM*closed_interval_rand(0,1))-1;
			follower = (follower < 0) ? 0 : follower;
			follower = (follower >= THREAD_NUM) ? THREAD_NUM-1 : follower;
								
			double rnd = closed_interval_rand(0,1);
			rnd = closed_interval_rand(0,1);
			
			opselected[tid] = 0;	
			j = 0;			
			if(pr[0]>pr[1] && pr[0] > pr[2])
				j = 0;
			else if(pr[1]>pr[0] && pr[1]>pr[2])
				j = 1;
			else if(pr[2]>pr[0] && pr[2]>pr[1])
				j = 2;	

			
			if(rnd < 1 - PMIN) 
			{
				opselected[tid] = j;
			}
			else
			{	
				opselected[tid]=ceil(3*closed_interval_rand(0,1))-1;
				
				while(j==opselected[tid])
					opselected[tid]=ceil(3*closed_interval_rand(0,1))-1;
				
				opselected[tid] = (opselected[tid] < 0) ? 0 : opselected[tid];
				
				
			}
			
			
			if(opselected[tid] == 0) 
			{
				/****************** BINCSA OP1 ******************/
				#pragma omp atomic
				opcount[0]++;				
				opselected[tid] = 0;
				
				
				for(j=0; j<loc; j++)
				{
					
					ind_agents_new[tid][j] = ( ind_agents_best[tid][j] ^ ( ( (rand() ^ tid ^ tid) & 1) & ( (ind_agents_best[follower][j] ^ ind_agents_best[tid][j]) ) ) );					
					
				}
				
				obj_agents_new[tid] = UFLP(loc, cus, customer, location, ind_agents_new[tid]);
				
				// Reward
				reward[tid][iter%W] = (loc / obj_agents_best[best_agent]) * (obj_agents[tid]-obj_agents_new[tid]);
				reward[tid][iter%W] = (reward[tid][iter%W] < 0) ? 0 : reward[tid][iter%W];		
				
				q = 0.0;
				if(RT == 1)
				{
					q=reward[tid][iter%W];
				}
				else if(RT == 2)
				{
					// AVG
					for(i=0; i<W; i++)
					{
						q+=reward[tid][i];
					}
					q = q / W;
				}
				else if(RT == 3)
				{
					// EXTREME 				
					for(i=0; i<W; i++)
					{	
						if(reward[tid][i]>q)
							q=reward[tid][i];
					}
				}
				
				
				
				#pragma omp critical
				c[0] = (1-ALPHA) * c[0] + ALPHA * q;
				
				if(reward[tid][iter%W] !=0)
				{
					obj_agents[tid] = obj_agents_new[tid];
						for(j=0; j<loc; j++)
							ind_agents[tid][j] = ind_agents_new[tid][j];
					
					#pragma omp atomic
					sr[0]++;
				}
					

				
			}
			else if(opselected[tid] == 1)  
			{
				/****************** IBINABC OP2 ******************/
				#pragma omp atomic
				opcount[1]++;
				opselected[tid] = 1;
				int dr = loc * (0.8 - ( ((0.8 - 0.2) / (MAX_ITER/THREAD_NUM)) * iter ));
				
								
				int darr[dr];
				for(j=0; j<dr; j++) 
					darr[j] = j;
				shuffle(darr,dr);
								
				for(j=0; j<dr; j++)
				{
					
					if(obj_agents_best[follower]<obj_agents_best[tid])
						ind_agents_new[tid][darr[j]] = ind_agents_best[follower][darr[j]];
					else
						ind_agents_new[tid][darr[j]] = ind_agents_best[tid][darr[j]];	// BEST_AGENT or TID	
					
				}
				
				obj_agents_new[tid] = UFLP(loc, cus, customer, location, ind_agents_new[tid]);				
				
				// Reward
				reward[tid][iter%W] = (loc / obj_agents_best[best_agent]) * (obj_agents[tid]-obj_agents_new[tid]);
				//reward[tid] = (obj_agents[tid] / obj_agents_best[best_agent]) * (obj_agents[tid]-obj_agents_new[tid]);
				reward[tid][iter%W] = (reward[tid][iter%W] < 0) ? 0 : reward[tid][iter%W];	
				q = 0.0;
				if(RT == 1)
				{
					q=reward[tid][iter%W];
				}
				else if(RT == 2)
				{
					// AVG
					for(i=0; i<W; i++)
					{
						q+=reward[tid][i];
					}
					q = q / W;
				}
				else if(RT == 3)
				{
					// EXTREME 				
					for(i=0; i<W; i++)
					{	
						if(reward[tid][i]>q)
							q=reward[tid][i];
					}
				}
				
				#pragma omp critical
				c[0] = (1-ALPHA) * c[0] + ALPHA * q;
				
				if(reward[tid][iter%W] !=0)
				{
					obj_agents[tid] = obj_agents_new[tid];
						for(j=0; j<loc; j++)
							ind_agents[tid][j] = ind_agents_new[tid][j];
					
					#pragma omp atomic
					sr[1]++;
				}
			}
			else if(opselected[tid] == 2)
			{				
				/****************** ELITEMUT OP3 ******************/
								
				#pragma omp atomic
				opcount[2]++;
				opselected[tid] = 2;
				
				
				for(j=0; j<loc; j++)
				{
					
					ind_dissim_xi[j] = ind_agents[tid][j];				
					ind_dissim_xk[j] = ind_agents[follower][j];					
				} 
				
				// Dice
				CalcDisSim(loc, ind_dissim_xi, ind_dissim_xk, ind_agents_new, tid, follower, MAX_ITER/THREAD_NUM, iter);
				
				obj_agents_new[tid] = UFLP(loc, cus, customer, location, ind_agents_new[tid]);
				
				// Reward
				reward[tid][iter%W] = (loc / obj_agents_best[best_agent]) * (obj_agents[tid]-obj_agents_new[tid]);
				reward[tid][iter%W] = (reward[tid][iter%W] < 0) ? 0 : reward[tid][iter%W];	
				q = 0.0;
				if(RT == 1)
				{
					q=reward[tid][iter%W];
				}
				else if(RT == 2)
				{
					// AVG
					for(i=0; i<W; i++)
					{
						q+=reward[tid][i];
					}
					q = q / W;
				}
				else if(RT == 3)
				{
					// EXTREME 				
					for(i=0; i<W; i++)
					{	
						if(reward[tid][i]>q)
							q=reward[tid][i];
					}
				}
				
				#pragma omp critical
				c[0] = (1-ALPHA) * c[0] + ALPHA * q;
				
				if(reward[tid][iter%W] !=0)
				{
					obj_agents[tid] = obj_agents_new[tid];
						for(j=0; j<loc; j++)
							ind_agents[tid][j] = ind_agents_new[tid][j];
						
					#pragma omp atomic
					sr[2]++;
				}
			}
			else
				printf("\n error detected");
					
			#pragma omp barrier		
			#pragma omp critical			
			if(obj_agents[tid] < obj_agents_best[tid])
            {
                obj_agents_best[tid] = obj_agents[tid];

                for(j=0; j<loc; j++)
                    ind_agents_best[tid][j] = ind_agents[tid][j];
				
				if (obj_agents_best[tid] < global_best)
				{
					global_best = obj_agents_best[tid];
					best_agent = tid;
					
					control = 1;
				}
            }
			
			//#pragma omp barrier	
			if (tid==0)
			{
				if (control == 1)
				{
					shrink = 0;
					control = 0;
				}
				else
				{	
					#pragma omp atomic
					shrink++;
				}
			}
				
			
			#pragma omp single
			{				
				pr[0] = 0.1 + (1 - 3 * 0.1) * (c[0] / (c[0]+c[1]+c[2]));
				pr[1] = 0.1 + (1 - 3 * 0.1) * (c[1] / (c[0]+c[1]+c[2]));
				pr[2] = 0.1 + (1 - 3 * 0.1) * (c[2] / (c[0]+c[1]+c[2]));
			}			
			
			#pragma omp barrier	
			if (shrink > SHAKE) 
			{
				if (tid==0)
				{
					if (global_best <= global_best_)
					{
						global_best_ = global_best;
						best_agent_ = best_agent;
						//printf("\tBest: %s \t %.3f (%d) | %.3f\n", argv[1], global_best_, best_agent_, ((double)(clock() - begin) / CLOCKS_PER_SEC));
					}
				}
								
				#pragma omp barrier
				
				
					for (j=0; j<loc; j++)
					{				
						if ( closed_interval_rand(0,1) < 0.5 )
						{
							ind_agents[tid][j] = 0;
							ind_agents_new[tid][j] = 0;
							ind_agents_best[tid][j] = 0;
						}
						else
						{
							ind_agents[tid][j] = 1;
							ind_agents_new[tid][j] = 1;
							ind_agents_best[tid][j] = 1;
						}	
					}
				
				
				obj_agents_best[tid] = UFLP(loc, cus, customer, location, ind_agents[tid]);
				obj_agents[tid] = obj_agents_best[tid];
				obj_agents_new[tid] = obj_agents_best[tid];		

				global_best = 9999999999;
				#pragma omp barrier
				if (obj_agents_best[tid] < global_best)
				{
					global_best = obj_agents_best[tid];
					best_agent = tid;
				}		
				if (tid==0)
					shrink = 0;
			}					
		}
	}
	if (global_best < global_best_)
	{
		global_best_ = global_best;
		best_agent_ = best_agent;
	}
	
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	

    //printf("PMIN: %.2f  \t W___: %d \t RT__: %d", PMIN, W, RT);
    printf("\t%s \t %.3f \t (%.3f)\n", argv[1], global_best_, time_spent);
    
    return 0;
}
