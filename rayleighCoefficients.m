function [alpha, beta] = rayleighCoefficients(freq1,damping1,freq2,damping2)

%Frequency in hertz

T1=1/freq1;
T2=1/freq2;

alpha=4*pi*(damping1*T1-damping2*T2)/(T1^2-T2^2);
beta=(T1*T2*(damping1*T1-damping2*T2))/(pi*(T1^2-T2^2));

%Alpha and Beta Coefficient formulas from:
%'Dynamic Time-History Elastic Analysis of Steel Frames Using One Element
%Per Member' doi http://dx.doi.org/10.1016/j.istruc.2016.05.006