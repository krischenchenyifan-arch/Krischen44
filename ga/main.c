#include <stdio.h>
#include <stdlib.h>     // This holds the definition of RAND_MAX
#include <math.h>
#include <time.h>       // Needed to use the time to set the seed

// Return a normally distributed random fraction
// using the Box Muller transform.
// Prof. M. Smith, NCKU
float randn() {
    float TWO_PI = 2.0*M_PI;
    float u1 = rand() / (float)RAND_MAX;  // Generate a random fraction from 0 to 1.
    float u2 = rand() / (float)RAND_MAX;
    float z1 = sqrt(-2.0*log(u1))*cos(TWO_PI*u2); // This is on the Wiki Page!
    return z1; 
}

float fitness(float x) {
        x = fabs(x * sin(x) - 0.2 );
        return x;
}

float calcFittest(float *F, int NO_P){
    float bestF = 10000;
    int bestID = -1;
    for (int p =0; p<NO_P; p++){
        if (F[p] < bestF){
        bestF = F[p];
        bestID = p;
        }
    }
    return bestID;
}

int main() {
    // Using C/C++, we need to set the seed properly
    srand((unsigned int)time(NULL));
    int G = 0;
    int NO_G = 100;
    int NO_P = 40;
    int Fittest;   
    float x[NO_P], F[NO_P];

    for (int i=0; i<NO_P; i++){
        x[i] = rand() / (float)RAND_MAX;
        printf("x = %g\n", x[i]);
        F[i] = fitness(x[i]);
        printf("fitness = %g\n", F[i]); 
    } 
    Fittest = calcFittest(F, NO_P);
    printf("Fittest = %d\n", Fittest);
    while (G < NO_G){
        int p = 0;
        int ParantA;
        int ParantB;
        float x_baby, F_baby;
        float New_F[NO_P], New_x[NO_P];

        for (p; p<NO_P; p++){
            if( p == Fittest){
                ParantA = p;
                ParantB = Fittest;
                while(ParantB == Fittest){
                    ParantB = (int)(rand() / (float)RAND_MAX)*NO_P;
                }

            }else {
                ParantA = p;
                ParantB = Fittest;
            }
            x_baby = 0.5*(x[ParantA]+x[ParantB])+ 0.5*randn();
            F_baby = fitness(x_baby);
            if(F_baby < F[p]){
                New_x[p] = x_baby;
                New_F[p] = F_baby;
            }else{
                New_x[p] = x[p];
                New_F[p] = F[p];
            }
        }
        for (int p=0; p<NO_P; p++){
        x[p] = New_x[p];
        F[p] = New_F[p];
        }
        G = G+1;
        Fittest = calcFittest( F, NO_P);
    }
    printf("x_value = %g, Fitness = %g\n", x[Fittest], F[Fittest]);
    //printf("Rand = %f\n", randn());
}