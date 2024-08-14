# JFET outpus characteristics

This C project calculates the output characteristics of a JFET based on a SPICE model. In this case the model is the 2N3819 N-channel JFET. Optimization by Rd and Rs requires a nonlinear approximate solution. I use the newton method.

## SPICE model

.model 2N3819 NJF(Beta=1.304m Betatce=-.5 Rd=1 Rs=1 Lambda=2.25m Vto=-3 Vtotc=-2.5m Is=33.57f Isr=322.4f N=1 Nr=2 Xti=3 Alpha=311.7u Vk=243.6 Cgd=1.6p M=.3622 Pb=1 Fc=.5 Cgs=2.414p Kf=9.882E-18 Af=1 mfg=Vishay)

## Equations with optimization (in ohmic and saturation region)

$$\textbf{If } (V_{DS} < (V_{GS} - V_{TO})) \textbf{ , then: } I_{D} = \beta' * V_{DS} * (2 * (V_{GS} - V_{TO}') - V_{DS}) * (1 + \lambda * V_{DS})$$


$$\textbf{else: } I_{D} = \beta' * (V_{GS} - V_{TO}')^2 * (1 + \lambda * V_{DS})$$


$$V_{DS} = V_{DD} - I_{D} * (R_{D} + R_{S})$$


$$\beta' = \beta * (1 + \beta_{tce} * (T - T_{ref}))$$


$$V_{TO}' = V_{TO} * V_{TOtc} * (T - T_{ref})$$


## Example images

This image was created using Excel from raw output data.

![alt text](http://www.vargalaszlo.com/images/out/JFET_output_characteristics-01.jpg)

This image is part of the output of C program.

![alt text](http://www.vargalaszlo.com/images/out/JFET_output_characteristics-02.jpg)

## Source code

```C
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

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
    double V_GS;
    double T;
    double T_ref;
    double V_GSStep;
    double V_GSStepCount;
    double V_DSStep;
    double V_DSStepCount;
    double solverTolerance;
    double solverMaxIteration;
} calculation;

typedef struct _characteristic {
    double** I_D;
    double* V_GS;
    double* V_DS;
} characteristic;

typedef struct _values {
    jfet* _jfet;
    calculation* _calculation;
    characteristic* _characteristic;
} values;

double jfetOutputCharacteristics(double V_GS, double V_DS, double LAMBDA, double BETA, double V_TO) {
    return (V_DS < (V_GS - V_TO)) ? BETA * V_DS *(2 * (V_GS - V_TO) - V_DS)*(1 + LAMBDA * V_DS) : BETA * pow(V_GS - V_TO, 2) * (1 + LAMBDA * V_DS);
}

double solverForI_D(calculation* cal, jfet* fet, bool type) {
    /* type:    0 - for transfer characteristic
                1 - for output characteristics */
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
        double tempV_DS = ((type) ? cal->V_DS : cal->V_DD) - I_DActual * (fet->R_D + fet->R_S);
        I_DActual = jfetOutputCharacteristics(cal->V_GS, tempV_DS, fet->LAMBDA, BETACorrected, V_TOCorrected);

        iteration++;
    } while (fabs(I_DActual - I_DPrevius) > cal->solverTolerance && iteration < cal->solverMaxIteration);
    return I_DActual;
}

bool jfetOutputCharacteristicsMake(values* value, double V_GSLow, double V_GSUp, double V_DSLow, double V_DSUp) {
    value->_calculation->V_GSStepCount = abs((int)((V_GSUp - V_GSLow)/value->_calculation->V_GSStep)) + 1;
    value->_calculation->V_DSStepCount = abs((int)((V_DSUp - V_DSLow)/value->_calculation->V_DSStep)) + 1;
    if (value->_calculation->V_GSStepCount < 1 || value->_calculation->V_DSStepCount < 1)
        return false;
    value->_characteristic->V_DS = (double *)calloc(value->_calculation->V_DSStepCount, sizeof(double));
    value->_characteristic->V_GS = (double *)calloc(value->_calculation->V_GSStepCount, sizeof(double));
    if (!value->_characteristic->V_DS || !value->_characteristic->V_GS)
        return false;
    // fill the V_GS vector
    for (int i = 0 ; i < value->_calculation->V_GSStepCount ; i++) {
        *(value->_characteristic->V_GS + i) = V_GSLow + (i * value->_calculation->V_GSStep);
    }
    // fill the V_DS vector
    for (int i = 0 ; i < value->_calculation->V_DSStepCount ; i++) {
        *(value->_characteristic->V_DS + i) = V_DSLow + (i * value->_calculation->V_DSStep);
    }
    // make and fill the I_D matrix
    value->_characteristic->I_D = (double**)calloc(value->_calculation->V_DSStepCount, sizeof(double*));
    if (!value->_characteristic->I_D)
        return false;
    for (int i = 0; i < value->_calculation->V_DSStepCount; ++i) {
        value->_characteristic->I_D[i] = (double*)calloc(value->_calculation->V_GSStepCount, sizeof(double));
        if (!value->_characteristic->I_D[i])
            return false;
    }
    int i = 0, j;
    for (double d = V_DSLow; d <= V_DSUp ; d += value->_calculation->V_DSStep) {
        j = 0;
        for (double g = V_GSLow ; g <= V_GSUp ; g += value->_calculation->V_GSStep) {
            value->_calculation->V_DS = d;
            value->_calculation->V_GS = g;
            value->_characteristic->I_D[i][j] = solverForI_D(value->_calculation, value->_jfet, 1) * 1000;
            j++;
        }
        i++;
    }
return true;
}

void jfetOutputCharacteristicsShow(values* value) {
    printf ("There are in V and mA.\n\nV_DS\t");
    for (int i = 0 ; i < value->_calculation->V_GSStepCount ; i++) {
        printf("I_D(U_GS=%.1fV)\t", value->_characteristic->V_GS[i]);
    }
    printf("\n");
    for (int i = 0 ; i < value->_calculation->V_DSStepCount ; i++) {
        printf("%.2f\t", value->_characteristic->V_DS[i]);
        for (int j = 0 ; j < value->_calculation->V_GSStepCount ; j++) {
                printf("%.4f\t\t", value->_characteristic->I_D[i][j]);
        }
    printf("\n");
    }
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
    value._jfet = &_2N3819;
    value._calculation = &calc;
    value._characteristic = &output;

    /* calculating
    *  jfetOutputCharacteristicsMake(values*, V_GSLow, V_GSUp, VDSLow, VDSUp)
    *
    *  Example:
    *    -3.0    <=  V_GS    <=  0.0
    *    0.0     <=  V_DS    <=  10.0
    */
    jfetOutputCharacteristicsMake(&value, -3.0, 0.0, 0.0, 10.0);

    // showing
    jfetOutputCharacteristicsShow(&value);

    return 0;
}


```

