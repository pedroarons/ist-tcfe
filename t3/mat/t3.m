close all
clear all

R1 = 6000;
C = 5e-6;
f=50;
Vi = 230;
R2 = 24000;
n = 12;
Vs = Vi/n;
t=linspace(0, 10/f, 50000);
w=2*pi*f;
vS1 = Vs * cos(w*t);
vO2 = t*0;
vO1 = t*0;

tOFF = 1/w * atan(1/w/R1/C);

vOn = Vs*cos(w*tOFF)*exp(-(t-tOFF)/R1/C);



vO2 = abs(vS1);


for i=1:length(t)
  if t(i) < tOFF
    vO1(i) = vO2(i);
  elseif vOn(i) > vO2(i)
    vO1(i) = vOn(i);
  else
    tOFF = tOFF + 1/(f*2);
    vOn = Vs*abs(cos(w*tOFF))*exp(-(t-tOFF)/R1/C);
    vO1(i) = vO2(i);
  endif
endfor

average_1 = mean(vO1);
ripple_1= max(vO1) - min(vO1);


n_diodes = 20;
v_on = 12/n_diodes;

vO2 = t*0;
vO2d = 0;
vO2a = t*0;


%dc component regulator ----------------
if average_1 >= v_on*n_diodes
  vO2d = v_on*n_diodes;
else
  vO2d = average_1;
endif

%ac component regulator -----------------
vt = 26e-3;
Is = 1e-14;
eta = 1;

rd = eta*vt/(Is*exp(v_on/(eta*vt)))

% ac regulator
for i = 1:length(t)
  if vO1(i) >= n_diodes*v_on
    vO2a(i) = n_diodes*rd/(n_diodes*rd+R2) * (vO1(i)-average_1);
  else
    vO2a(i) = vO1(i)-average_1;
  endif
endfor

vO2 = vO2d + vO2a;



hfa = figure(1);
title('Regulator and envelope output voltage v_o(t)')
plot (t*1000, vS1, ";vs_{transformer}(t);", t*1000,vO1, ";vo_{envelope}(t);", t*1000,vO2, ";vo_{regulator}(t);");
xlabel ("t[ms]")
ylabel ("v_O [Volts]")
legend('Location','northeast');
print (hfa, "all_vout.eps", "-depsc");

hfb = figure(2);
title('Deviations from desired DC voltage')
plot (t*1000,(vO2-12)*1000, ";vo_{regulator}-12(t);");
xlabel ("t[ms]")
ylabel ("v_O [mV]")
legend('Location','northeast');
print (hfb, "deviation.eps", "-depsc");


diary result_octave.txt
max(vO2)-min(vO2) 
diary on
average_2 = mean(vO2)
max(vO2)
min(vO2)
ripple_2 = max(vO2)-min(vO2) 
diary off

cost = R1/1000 + R2/1000 + C*1e6 + n_diodes*0.1*2; 


cost = cost + 0.4;


MERIT = 1/(cost*(ripple_2 + abs(average_2 - 12) + 1e-6));
v4t = '"sim_TAB"'
v4e = '"sim_END"'
v5t = '"v5_TAB"'
v5e = '"v5_END"'
mt = '"merit_TAB"'
me = '"merit_END"'
file4sim = fopen ("t3.cir", "w");
fprintf(file4sim,"Ngspice T3 \nVs 1 0 0 sin(0 230 50 0 0 90) \nD1 2 4 Default \nD2 GND 2 Default \nD3 3 4 Default \nD4 GND 3 Default \nR1 GND 4 6k \nR2 5 4 24k \nC1 4 GND 5u \nD5 5 GND Default2 \nF1 1 0 E1 0.083333333333333333333 \nE1 3 2 1 0 0.083333333333333333333 \n.model Default D \n.model Default2 D (n=20) \n.end \n.control \nset hcopypscolor=0 \nset color0=white \nset color1=black \nset color2=red \nset color3=blue \nset color4=violet \nset color5=rgb:3/8/0 \nset color6=rgb:4/0/0 \nop \ntran 0.0002 400m 200m \n");
fprintf(file4sim,"let RippleEnvelope = maximum(v(4))-minimum(v(4)) \n");
fprintf(file4sim,"let AverageEnvelope = mean(v(4)) \n");
fprintf(file4sim,"let RippleRegulator = maximum(v(5))-minimum(v(5)) \n");
fprintf(file4sim,"let AverageRegulator = mean(v(5)) \n");
fprintf(file4sim,"echo %s \n", v4t);
fprintf(file4sim,"print RippleEnvelope \nprint AverageEnvelope \n");
fprintf(file4sim,"print RippleRegulator \nprint AverageRegulator \n");
fprintf(file4sim,"echo %s \n", v4e);
fprintf(file4sim,"let Merit = 1/ (37.4* ((maximum(v(5))-minimum(v(5))) + abs(mean(v(5)-12)) + 10e-6)) \n");
fprintf(file4sim,"echo %s \n", mt);
fprintf(file4sim,"print Merit \n");
fprintf(file4sim,"echo %s \n", me);
fprintf(file4sim,"let v4 = v(4)-2*v(2)[0] \nhardcopy comp.ps v(3)-v(2) v4 v(5) \necho comp_FIG \nhardcopy dev.ps v(5)-12 \necho dev_FIG \nquit \n.endc");
fclose (file4sim);

RipAvg = fopen ("RipAvg.tex", "w");
fprintf(RipAvg, "Ripple Envelope & %e \\\\ \\hline \n", ripple_1);
fprintf(RipAvg, "Average Envelope & %e \\\\ \\hline \n", average_1);
fprintf(RipAvg, "Ripple Regulator & %e \\\\ \\hline \n", ripple_2);
fprintf(RipAvg, "Average Regulator & %e \\\\ \\hline \n", average_2);
fclose (RipAvg);

MeritTable = fopen ("MeritTable.tex", "w");
fprintf(MeritTable, "Merit & %e \\\\ \\hline \n", MERIT);
fclose (MeritTable);
