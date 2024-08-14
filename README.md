# JFET transfer characteristic.

This C project calculates the output characteristics of a JFET based on a SPICE model. In this case the model is the 2N3819 N-channel JFET. Optimization by Rd and Rs requires a nonlinear approximate solution. I use the newton method.

## SPICE model

.model 2N3819 NJF(Beta=1.304m Betatce=-.5 Rd=1 Rs=1 Lambda=2.25m Vto=-3 Vtotc=-2.5m Is=33.57f Isr=322.4f N=1 Nr=2 Xti=3 Alpha=311.7u Vk=243.6 Cgd=1.6p M=.3622 Pb=1 Fc=.5 Cgs=2.414p Kf=9.882E-18 Af=1 mfg=Vishay)

## Equations with optimization (in ohmic and saturation region)

$$\textbf{If } (V_{DS} < (V_{GS} - V_{TO})) \textbf{ , then: }\beta' * V_{DS} *(2 * (V_{GS} - V_{TO}') - V_{DS})*(1 + \lambda * V_{DS})$$


$$\textbf{else: } I_{D} = \beta' * (V_{GS} - V_{TO}')^2 * (1 + \lambda * V_{DS})$$


$$V_{DS} = V_{DD} - I_{D} * (R_{D} + R_{S})$$


$$\beta' = \beta * (1 + \beta_{tce} * (T - T_{ref}))$$


$$V_{TO}' = V_{TO} * V_{TOtc} * (T - T_{ref})$$


## Example images

This image was created using Excel from raw output data.

![alt text2](http://www.vargalaszlo.com/images/out/JFET_output_characteristics-01.jpg)
