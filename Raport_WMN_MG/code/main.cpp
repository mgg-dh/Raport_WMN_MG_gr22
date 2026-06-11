#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "winbgi2.h"
#include "rk4.h"

#define G 9.81
#define L 1.0
#define M 0.1
#define n 2

#define K1 5.0   
#define K3 5.0   

void vrhs(double x0, double* y0, double* y1) {
    
    y1[0] = y0[1];

    double x = L * sin(y0[0]);

   
    y1[1] = -(G / L) * sin(y0[0]) - (cos(y0[0]) / (M * L)) * (K1 * x + K3 * x * x * x);
}


double energia(double alpha, double alpha_prim) {
    double x = L * sin(alpha);

    double E_k = 0.5 * M * L * L * alpha_prim * alpha_prim;     
    double E_pg = M * G * L * (1.0 - cos(alpha));               
    double E_ps = 0.5 * K1 * x * x + 0.25 * K3 * x * x * x * x; 

    return E_k + E_pg + E_ps;
}

int main() {
    
    graphics(500, 500);
    scale(0, -3, 10, 3);

    double y[2];
   
    y[0] = 0.5;
    y[1] = 0.0;

    double y1[2];
    double dt = 0.01;

    
    FILE* fp = fopen("wyniki_rk4.csv", "w");
    if (fp == NULL) {
        printf("Blad krytyczny: Nie mozna utworzyc pliku wynikowego.\n");
        return 1;
    }

    
    fprintf(fp, "czas,alpha,omega,energia\n");

    for (double t = 0; t <= 10; t += dt) {
        vrk4(t, y, dt, n, vrhs, y1);

        fprintf(fp, "%lf,%lf,%lf,%lf\n", t, y1[0], y1[1], energia(y1[0], y1[1]));

       
        setcolor(0.1);
        point(t, y1[0]); 

        setcolor(0.4);
        point(t, y1[1]); 

        setcolor(0.7);
        point(t, energia(y1[0], y1[1])); 

       
        printf("t = %2.2lf s | E_mech = %lf J\n", t, energia(y1[0], y1[1]));

        
        y[0] = y1[0];
        y[1] = y1[1];
    }

    fclose(fp);
    printf("Symulacja zakonczona sukcesem. Dane zapisano do pliku 'wyniki_rk4.csv'.\n");

    wait();
    return 0;
}