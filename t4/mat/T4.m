%Gain Stage%

VT = 25e-3;
BFN = 178.7;
VAFN = 69.7;
RE1 = 100;
RC1 = 2500;
RB1 = 20000;
RB2 = 1500;
VBEON = 0.7;
VCC = 12;
RS = 100;
Ci = 300e-06;
CB = 300e-06;
Co = 200e-6;
Vin = 0.01;
f = logspace(1,8, 100);
w = 2*pi*f;
Zci = 1./(j*w*Ci);
ZcB = 1./(j*w*CB);
RB=1/(1/RB1+1/RB2);
VEQ=RB2/(RB1+RB2)*VCC;
IB1=(VEQ-VBEON)/(RB+(1+BFN)*RE1);
IC1=BFN*IB1;
IE1=(1+BFN)*IB1;
VE1=RE1*IE1;
VO1=VCC-RC1*IC1;
VCE=VO1-VE1;
printf("valoresiniciais_TAB\n");
printf("R1 = %e \n", RB1);
printf("R2 = %e \n", RB2);
printf("RC = %e \n", RC1);
printf("Re = %e \n", RE1);
printf("Rout = %e \n", RE1);
printf("Ci = %e \n", Ci);
printf("CB = %e \n", CB);
printf("Co = %e \n", Co);
printf("valoresiniciais_END\n\n");
gm1=IC1/VT;
rpi1=BFN/gm1;
ro1=VAFN/IC1;


for k=1:length(w)
R_sum = RB+RS;
Z_sum = Zci(k)+ZcB(k); 
G1 = [RE1 + ro1 + RC1, 0, -ro1, -RE1, 0;
    0, -RB, 0, 0, R_sum;
    -RE1, -ZcB(k), 0, RE1 + ZcB(k), 0; 
    0, RB+Z_sum+rpi1, 0, -ZcB(k), -RB
    0, gm1*rpi1, 1, 0, 0;];
G2 = [0; Vin; 0; 0; 0];
finalg = G1\G2;

Q2E1(k) = abs(RC1 * finalg(1));
AV1(k) = Q2E1(k)/Vin;
endfor
ZI1 = 1/(1/RB + 1/rpi1);
ZO1 = 1/(1/ro1+1/RC1);

GS = fopen ("GS.tex", "w");
fprintf(GS, "Input impedance & %e Ohm \\\\ \\hline \n \n", ZI1);
fprintf(GS, "Output impedance & %e Ohm \\\\ \\hline \n \n", ZO1);
fprintf(GS, "Gain & %e V \\\\ \\hline \n \n", mean(abs(AV1)));
fclose (GS);


%Output Stage%
BFP = 227.3;
VAFP = 37.2;
RE2 = 100;
VEBON = 0.7;
VI2 = VO1;
IE2 = (VCC-VEBON-VI2)/RE2;
IC2 = BFP/(BFP+1)*IE2;
VO2 = VCC - RE2*IE2;

gm2 = IC2/VT;
go2 = IC2/VAFP;
gpi2 = gm2/BFP;
ge2 = 1/RE2;
RL=8;

I2 = zeros(length(w),3);
vo2 = zeros(1,length(w));
AV2 = zeros(1,length(w));

for k=1:length(w)
O1 = [1/gpi2+RE2,-RE2, 0;
      gm2*1/gpi2, 0, 1; 
      -RE2, RE2+1/go2, -1/go2];

O2 = [Q2E1(k); 0; 0];
finalo = O1\O2;

Q2E2(k) = abs((finalo(1)-finalo(2))*RE2);
AV2(k) = Q2E2(k)/Q2E1(k);
endfor

gB = 1/(1/gpi2+ZO1)
ZI2 = (gm2+gpi2+go2+ge2)/gpi2/(gpi2+go2+ge2);
ZO2 = 1/(gm2+gpi2+go2+ge2);
ZO = 1/(go2+gm2/gpi2*gB+ge2+gB)

AVtot = (gB+gm2/gpi2*gB)/(gB+ge2+go2+gm2/gpi2*gB)*AV1

fl = log10(f);
gB = 1/(1/gpi2+ZO1)
ZI=ZI1
ZOt=1/(go2+gm2/gpi2*gB+ge2+gB)

OS = fopen ("OS.tex", "w");
fprintf(OS, "Input impedance & %e Ohm \\\\ \\hline \n \n", ZI2);
fprintf(OS, "Output impedance & %e Ohm \\\\ \\hline \n \n", ZO2);
fprintf(OS, "Gain & %e V \\\\ \\hline \n \n", mean(abs(AV2)));
fclose (OS);

graph = figure();
plot(fl,20*log10(AVtot),'r');
xlabel ("Frequency (Hz)");
ylabel ("Gain");
print(graph, "gain.eps", "-color");

TOTAL = fopen ("TOTAL.tex", "w");
fprintf(TOTAL, "Input impedance & %e Ohm \\\\ \\hline \n \n", ZI);
fprintf(TOTAL, "Output impedance & %e Ohm \\\\ \\hline \n \n", ZO);
fprintf(TOTAL, "Gain & %e V \\\\ \\hline \n \n", mean(abs(AVtot)));
fprintf(TOTAL, "Gain & %e dB \\\\ \\hline \n \n", 20*log10(mean(abs(AVtot))));
fclose (TOTAL);

