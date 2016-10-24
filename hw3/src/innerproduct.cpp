/* Name: Tanmoy Sanyal
   Perm: 7550049
 */

#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>

#include <algorithm>
#include <numeric>
#include <string>
#include <iostream>
#include <cmath>
#include <iterator>
#include <functional>

#include "example_util_gettime.h"

#include <stdio.h>

#define COARSENESS input_coarseness
#define ITERS 10

int input_coarseness = 50;

double rec_cilkified(double * a, double * b, int n)
{
	double r1 = 0.0;
	double r2 = 0.0;
	if (n > COARSENESS) {
		r1 = cilk_spawn rec_cilkified(a, b, n/2);
		r2 = rec_cilkified(a + n/2, b + n/2, n - n/2);
		cilk_sync;
		return (r1+r2);
	}

	for (int i=0; i < n; i++) {
		r1 += a[i] * b[i];
	}

	return r1;
}


double loop_cilkified(double * a, double * b, int n)
{
	double result = 0.0; double inner_result = 0.0;
    	cilk_for (int i = 0; i < (n/COARSENESS); i++) {
    		inner_result = 0.0;
    		for (int j = 0; j < COARSENESS; j++)
    			inner_result += a[i*COARSENESS + j] * b[i*COARSENESS + j];
    		result += inner_result;
    }
	return result;
}


double hyperobject_cilkified(double * a, double * b, int n)
{
	cilk::reducer_opadd<double> result;
	cilk_for (int i = 0; i < n; i++)
		result += a[i] * b[i];
    
    return result.get_value();
}


int close(double x, double y, int n)
{
        double relative_diff = (x>y? (x-y)/x : (y-x)/x);
        return (relative_diff < sqrt((double) n) * exp2(-42))? 1 : 0;
}


// A simple test harness 
int inn_prod_driver(int n)
{
	double * a = new double[n];
	double * b = new double[n];
	for (int i = 0; i < n; ++i)
	{
        	a[i] = i;
		b[i] = i;
	}
   	//std::random_shuffle(a, a + n);
	//std::random_shuffle(b, b + n);

	double seqresult = std::inner_product(a, a+n, b, (double)0 );	

	long t1 = example_get_time();
	for(int i=0; i< ITERS; ++i)
	{
		seqresult = std::inner_product(a, a+n, b, (double)0);	
	}
	long t2 = example_get_time();

	double seqtime = (t2-t1)/(ITERS*1000.f);
	std::cout << "Sequential time: " << seqtime << " seconds" << std::endl;	
	
	/***********************************************************/
	/********  START TESTING RECURSIVE CILKFIED VERSION  *******/
	/***********************************************************/

	double parresult = rec_cilkified(a, b, n);   
	t1 = example_get_time();
	for(int i=0; i< ITERS; ++i)
	{
		parresult = rec_cilkified(a, b, n);   
	}
 	t2 = example_get_time();

	double partime = (t2-t1)/(ITERS*1000.f);
	std::cout << "Recursive cilkified time:" << partime << " seconds" << std::endl;
	std::cout << "Speedup is: " << seqtime/partime << std::endl;
	std::cout << "Sequential result is: "<<seqresult<<std::endl;
	std::cout << "Recursive cilkified result is: "<<parresult<<std::endl;
	std::cout << "Result is " << (close(seqresult,parresult,n)  ? "correct":"incorrect") << std::endl; 
	
	/****************************************************************/
	/********  START TESTING NESTED LOOPED CILKIFIED VERSION  *******/
	/****************************************************************/
	parresult = loop_cilkified(a, b, n);   
	
	t1 = example_get_time();
	for(int i=0; i< ITERS; ++i)
	{
		parresult = loop_cilkified(a, b, n);   
	}
 	t2 = example_get_time();


	partime = (t2-t1)/(ITERS*1000.f);
	std::cout << "Nested loop cilkified time: " << partime << " seconds" << std::endl;
	std::cout << "Speedup is: " << seqtime/partime << std::endl;
	std::cout << "Sequential result is: "<<seqresult<<std::endl;
	std::cout << "Loop cilkified result is: "<<parresult<<std::endl;
	std::cout << "Result is " << (close(seqresult,parresult,n)  ? "correct":"incorrect") << std::endl; 
	
	/**************************************************************/
	/********  START TESTING HYPEROBJECT CILKIFIED VERSION  *******/
	/**************************************************************/

	parresult = hyperobject_cilkified(a, b, n);   
	
	t1 = example_get_time();
	for(int i=0; i< ITERS; ++i)
	{
		parresult = hyperobject_cilkified(a, b, n);   
	}
 	t2 = example_get_time();

	partime = (t2-t1)/(ITERS*1000.f);
	std::cout << "Hyperobject cilkified time:" << partime << " seconds" << std::endl;
	std::cout << "Speedup is: " << seqtime/partime << std::endl;
	std::cout << "Sequential result is: "<<seqresult<<std::endl;
	std::cout << "Hyperobject result is: "<<parresult<<std::endl;
	std::cout << "Result is " << (close(seqresult,parresult,n)  ? "correct":"incorrect") << std::endl; 
    	
        
	delete [] a;
	delete [] b;
    	return 0;
}

int main(int argc, char* argv[])
{
    int n = 1 * 100 * 100;
    if (argc > 1) {
        n = std::atoi(argv[1]);
    }
   
    // user input coarseness
    if (argc > 2) {
	input_coarseness = std::atoi(argv[2]);
    }

    return inn_prod_driver(n);
}
