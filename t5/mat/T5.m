C1 = 220e-9
C2 = 180.33e-9
R1 = 1000
R2 = 1000
R3 = 150000
R4 = 1000
%%%
f1 = 1000
w1 = 2*pi*f1
Zin = (1./(j*w1*C1)+R1)
Zout = (1/(1/R2+1/(1./(j*w1*C2))))
AbsZin = abs(1./(j*w1*C1)+R1)
AbsZout = abs(1/(1/R2+1/(1./(j*w1*C2))))
Va = (1+R3/R4)*(R1/(R1+1./(j*w1*C1)))
Gain = abs((1./(j*w1*C2))/((1./(j*w1*C2))+R2)*Va)
Gaindb = 20*log10(Gain)
f2 = logspace(1, 8, 1000);
w2 = 2*pi*f2;
Va = (1+R3/R4).*(R1./(R1+1./(j*w2*C1)));
G = (1./(j*w2*C2))./((1./(j*w2*C2))+R2).*Va;
AbsG = abs(G);
Gdb = 20*log10(AbsG);
Phase = angle(G)*180/pi;
m = max(Gdb)
  fl = 0;
for i=1:length(Gdb)
	if (Gdb(i) >= m-3 && !fl)
	  lw = f2(i-1)
	    fl = 1;
	endif
	if (Gdb(i) <= m-3 && fl)
	  hw = f2(i-1)
	    fl = 0;
  endif
endfor
cw = sqrt(lw*hw)	
Val = fopen ("Val.tex", "w");
fprintf(Val, "$C_{1}$ & %d $\\mu F$ \\\\ \\hline \n \n", C1);
fprintf(Val, "$C_{2}$ & %d $\\mu F$ \\\\ \\hline \n \n", C2);
fprintf(Val, "$R_{1}$ & %d $\\Omega$ \\\\ \\hline \n \n", R1);
fprintf(Val, "$R_{2}$ & %d $\\Omega$ \\\\ \\hline \n \n", R2);
fprintf(Val, "$R_{3}$ & %d $\\Omega$ \\\\ \\hline \n \n", R3);
fprintf(Val, "$R_{4}$ & %d $\\Omega$ \\\\ \\hline \n \n", R4);
fclose (Val);
TA = fopen ("TA.tex", "w");
fprintf(TA, "$|Z_{in}|$ & %f $\\Omega$ \\\\ \\hline \n \n", AbsZin);
fprintf(TA, "$|Z_{out}|$ & %f $\\Omega$ \\\\ \\hline \n \n", AbsZout);
fprintf(TA, "Gain & %f dB \\\\ \\hline \n \n", Gaindb);
fclose (TA);
Comp = fopen ("Comp.tex", "w");
fprintf(Comp, "$Z_{in}$  & %d + %di $\\Omega$ \\\\ \\hline \n \n", real(Zin), imag(Zin));
fprintf(Comp, "$|Z_{in}|$ & %f $\\Omega$ \\\\ \\hline \n \n", AbsZin);
fprintf(Comp, "$Z_{out}$ & %d + %di $\\Omega$ \\\\ \\hline \n \n", real(Zout),imag(Zout));
fprintf(Comp, "$|Z_{out}|$ & %f $\\Omega$  \\\\ \\hline \n \n", AbsZout);
fprintf(Comp, "Gain & %f dB \\\\ \\hline \n \n", Gaindb);
fprintf(Comp, "Low Frequency & %f Hz \\\\ \\hline \n \n", lw);
fprintf(Comp, "High Frequency & %f Hz \\\\ \\hline \n \n", hw);
fprintf(Comp, "Central Frequency & %f Hz \\\\ \\hline \n \n", cw);
fclose (Comp);
gain = figure ();
plot(log10(f2),Gdb,"r");
title("Gain Frequency Response");
xlabel ("Frequency [Hz]");
ylabel ("Gain [dB]");
legend("v_o(f)/v_i(f)");
print (gain, "gain.eps", "-depsc");
phase = figure();
plot(log10(f2),Phase,"g");
title("Phase Frequency Response");
xlabel("Frequency [Hz]");
ylabel("Phase [Deg]");
print(phase, "phase.eps", "-depsc");
