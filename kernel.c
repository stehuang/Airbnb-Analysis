
#include <R.h>
#include <Rmath.h>

// problem1


double pareto_density(double x, double c, double p){
	double fx;
	fx = (p-1) * pow(c, p-1) * pow(x+c, -p);
	return(fx);
}


//double pareto_simp(double start,double stop,long n, double c, double p)
//{
//    double mult,x,t,inc;
//    long i;
//    inc = (stop - start) / (double)n; // interval length
//    x = start;
//    t = pareto_density(x, c, p); //f(x)
//    mult = 4.0;
//    for(i=1; i<n; i++){
//		x += inc;
//		t += mult * pareto_density(x, c, p);
//		if(mult == 4.0) mult = 2.0;
//		else mult = 4.0;
//    }
//    t += pareto_density(stop, c, p);
//    return(t * inc / 3.0);
//}


void paretoint(int *xmax, double *c, double *p, double *y){
	int n = 1000000;
    double start = 0;
    double stop = *xmax;
    
    double mult,x,t,inc;
        long i;
        inc = (stop - start) / (double)n; // interval length
        x = start;
        t = pareto_density(x, *c, *p); //f(x)
        mult = 4.0;
        for(i=1; i<n; i++){
            x += inc;
            t += mult * pareto_density(x, *c, *p);
            if(mult == 4.0) mult = 2.0;
            else mult = 4.0;
        }
        t += pareto_density(stop, *c, *p);
        *y = (t * inc / 3.0);

}




// problem 2


void kernel_reg(int *n, double *x, double *y, double *b, int* m, 
	double *g2, double *res2)
{
	int i,j;
	double c;
	double est_sum_y;
	double est_sum;


	// loop through each grid index
    for(i=0; i<*m; i++){
    	est_sum = 0.0;
    	est_sum_y = 0.0;

    	// calculate density; loop through each data index
    	for(j=0; j<*n; j++){
    		c = dnorm(x[j]-g2[i], 0, *b, 0)/ *n; //compute kernel at x on grid
    		est_sum += c; // denominator
    		est_sum_y += y[j] * c; //multiply by y and add; numerator
    	}
    	
    	res2[i] = est_sum_y / est_sum;
    }

}


// // for some kernel function
// void kernreg2 (double *x, double *y, int *n, double *b,
// 	       double *g, int *m, double *est)
// {
//     int i,j;
//     double a1,a2,c;
//     for(i = 0; i < *m; i++){
// 		a1 = 0.0;
// 		a2 = 0.0;
// 		for(j=0; j < *n; j++){
// 		    if(fabs((x[j] - g[i])/ *b) <= 1.0){
// 			c = 0.75 * (1.0-pow((x[j]-g[i])/ *b,2.0))/ *b;
// 			a1 += y[j] * c;	
// 			a2 += c;
// 		    }
// 		}
// 		if(a2 > 0.0) est[i] = a1/a2; 
// 		else est[i] = 0.0;
// 		Rprintf("I'm on gridpoint %f and index %d right now.\n", g[i], i);
// 	}
// }


// // for gaussian kernel
// void kernreg2 (double *x, int *n, double *b,
// 	       double *g, int *m, double *est)
// {

//     int i, j;
//     double est_sum;

//     // loop through each grid index
//     for(i=0; i<*m; i++){
//     	est_sum = 0.0;

//     	// calculate density; loop through each data index
//     	for(j=0; j<*n; j++){
//     		est_sum += dnorm(x[j]-g[i], 0, *b, 0)/ *n; //using dnorm in C
//     	}
//     	est[i] = est_sum;
//     }

// }






























