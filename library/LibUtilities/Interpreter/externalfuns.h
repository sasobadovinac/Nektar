#ifndef _EXTERNALFUNS_H
#define _EXTERNALFUNS_H
#include <iostream>
#include <cmath>
#include "blasius.h"

//------------------------------------------------------------------------------
//  For blasius profile
//------------------------------------------------------------------------------
inline double blau(double y){
        double t, a, b;
	int N,i; 
        int Nbla;
        Nbla = sizeof(eta)/sizeof(double);
        N = 1;
        // find the interval which contain Y
        if(y>=eta[Nbla-1]) return u[Nbla-1];
	if(y<=double(0)) return double(0);
        for(i=0;i<Nbla-1;i++)
	{
		if(y>=eta[i] && y<=eta[i+1])
		{
			N=i;
			i=Nbla+1;	
		}
	}
	// interpolation
	t =  (y-eta[N])/(eta[N+1]-eta[N]);
	a =  du[N]*(eta[N+1]-eta[N])-(u[N+1]-u[N]);
	b = -du[N+1]*(eta[N+1]-eta[N])+(u[N+1]-u[N]);
        
        return ((1-t)*u[N]+t*u[N+1]+t*(1-t)*(a*(1-t)+b*t));
}
inline double blav(double y){
        double t, a, b;
	int N,i; 
        int Nbla;
        Nbla = sizeof(eta)/sizeof(double);
        N = 1;
        // find the interval which contain Y
        if(y>=eta[Nbla-1]) return v[Nbla-1];
        if(y<double(0)) return double(0);
        for(i=0;i<Nbla-1;i++)
	{
		if(y>=eta[i] && y<=eta[i+1])
		{
			N=i;
			i=Nbla+1;		
		}
	}
	// interpolation
	t =  (y-eta[N])/(eta[N+1]-eta[N]);
	a =  dv[N]*(eta[N+1]-eta[N])-(v[N+1]-v[N]);
	b = -dv[N+1]*(eta[N+1]-eta[N])+(v[N+1]-v[N]);
        return ((1-t)*v[N]+t*v[N+1]+t*(1-t)*(a*(1-t)+b*t));
}
inline double blat(double y){
        double t, a, b;
	int N,i; 
        int Nbla;
        Nbla = sizeof(eta)/sizeof(double);
        N = 1;
        // find the interval which contain Y
        if(y>=eta[Nbla-1]) return tr[Nbla-1];
        if(y<double(0)) return tr[0];
        for(i=0;i<Nbla-1;i++)
	{
		if(y>=eta[i] && y<=eta[i+1])
		{
			N=i;
			i=Nbla+1;		
		}
	}
	// interpolation
	t =  (y-eta[N])/(eta[N+1]-eta[N]);
	a =  dtr[N]*(eta[N+1]-eta[N])-(tr[N+1]-tr[N]);
	b = -dtr[N+1]*(eta[N+1]-eta[N])+(tr[N+1]-tr[N]);
        return ((1-t)*tr[N]+t*tr[N+1]+t*(1-t)*(a*(1-t)+b*t));
}
inline double hump(double x){
	double re;
	if (x<=0){	
		re = double(0);
	}else{
         	if (x>=1){
		 	re = double(1);
		}else{
			re = 1.0/(1.0+exp(1.0/(x-1.0)+1.0/x));
		}
        }
	return re;
}
// Modified by Hui Xu 23 May 2013
// Added by Hui Xu 2013 Jan
#endif
