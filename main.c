#include <stdio.h>
#include <math.h>

typedef struct _jfet {
    double BETA;
    double V_TO;
    double LAMBDA;
    double R_D;
    double R_S;
    double V_TOtc;
    double BETA_tce;
} jfet;

typedef struct _calculation {
    double V_DD;
    double V_DS;
    double T;
    double T_ref;
    double V_GSStep;
    double V_DSStep;
    double solverTolerance;
    double solverMaxIteration;
} calculation;

typedef struct _characteristic {
    double I_D;
    double V_GS;
} characteristic;

typedef struct _values {
    jfet *_jfet;
    calculation *_calculation;
    characteristic *_characteristic;
} values;

double jfetOutputCharacteristics(double V_GS, double V_DS, double LAMBDA, double BETA, double V_TO) {
    if (V_GS < V_TO)
        return 0.0;
    return (V_DS < (V_GS - V_TO)) ? BETA * ((V_GS - V_TO) * V_DS - 0.5 * pow(V_DS,2)) : BETA * pow(V_GS - V_TO, 2) * (1 + LAMBDA * V_DS);
}

double solveI_D(calculation* cal, jfet* fet, double V_GSActual) {
    double
        I_DActual = 0.0,
        I_DPrevius,
        V_TOCorrected,
        BETACorrected;
    int iteration = 0;

    // correction
    V_TOCorrected = fet->V_TO + fet->V_TOtc * (cal->T - cal->T_ref);
    BETACorrected = fet->BETA * (1 + fet->BETA_tce * (cal->T - cal->T_ref));

    // Newton method
    do {
        I_DPrevius = I_DActual;
        double tempV_DS = cal->V_DS - I_DActual * (fet->R_D + fet->R_S);
        I_DActual = jfetOutputCharacteristics(V_GSActual, tempV_DS, fet->LAMBDA, BETACorrected, V_TOCorrected);

        iteration++;
    } while (fabs(I_DActual - I_DPrevius) > cal->solverTolerance && iteration < cal->solverMaxIteration);

    return I_DActual;
}

int main() {

  // 2N3819 JFET spice parameters (only necessary)
    jfet _2N3819;

    _2N3819.BETA = 1.304e-3;
    _2N3819.V_TO = -3.0;
    _2N3819.LAMBDA = 2.25e-3;
    _2N3819.R_D = 1;
    _2N3819.R_S = 1;
    _2N3819.BETA_tce = -0.5e-2;
    _2N3819.V_TOtc = -2.5e-3;

    // calculation setup
    calculation calc;

    calc.V_DD = 10;
    calc.T = 26.85;
    calc.T_ref = 26.85;
    calc.V_GSStep = 1.0;
    calc.V_DSStep = 0.1;
    calc.solverTolerance = 1e-6;
    calc.solverMaxIteration = 100;

    // characteristic
    characteristic output;

    // collecting structure
    values value;
    value._calculation = &calc;
    value._jfet = &_2N3819;
    value._characteristic = &output;


    printf ("There are in V and mA.\n\n");

    output.V_GS = -3.0;
    printf("V_DS\t");
    for (output.V_GS ; output.V_GS <= 0.0 ; output.V_GS += calc.V_GSStep) {
        printf("I_D(U_GS=%.1fV)\t", output.V_GS);
    }
    printf("\n");


    calc.V_DS = 0.0;
    for (calc.V_DS ; calc.V_DS <= calc.V_DD ; calc.V_DS += calc.V_DSStep) {
        printf("%.2f\t", calc.V_DS);
        output.V_GS = -3.0;
        for (output.V_GS ; output.V_GS <= 0.0 ; output.V_GS += calc.V_GSStep) {
                double I_D;
                //I_D = jfetOutputCharacteristics(output.V_GS, calc.V_DS, _2N3819.LAMBDA, _2N3819.BETA, _2N3819.V_TO);
                I_D = solveI_D(&calc, &_2N3819, output.V_GS);
                printf("%.4f\t\t", I_D * 1000); // I_D mA-ban
        }
        printf("\n");
    }
    //calculate_characteristic();
    return 0;
}
