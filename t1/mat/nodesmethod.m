R1 = 1.01080769792 
R2 = 2.07664633274 
R3 = 3.12595649013 
R4 = 4.18722214507 
R5 = 3.0841699201 
R6 = 2.00179338129 
R7 = 1.04556537884 
Va = 5.00120775651 
Id = 1.03136220214 
Kb = 7.12593545434 
Kc = 8.24048597287

A = [1 0 0 0 -1 0 0 ; 0 1/R2 0 0 0 Kb (-1/R2)-Kb ; 0 0 1/R5 0 0 (-1/R5)-Kb Kb ; 0 0 0 (1/R6)+(1/R7) -1/R6 0 0 ; -1/R1 -1/R2 0 0 0 -1/R3 (1/R1)+(1/R2)+(1/R3) ; 1/R1 0 0 -1/R6 (1/R4)+1/R6 -1/R4 -1/R1 ; 0 0 0 Kc/R6 -Kc/R6 1 0]
B = [Va ; 0 ; Id ; 0 ; 0 ; 0 ; 0]
V=A\B
Ia = (V(1)-V(7))/R1
Ib = (V(6)-V(7))*Kb
Ic = (V(5)-V(4))/R6
disp(V)
disp(Ia)
disp(Ib)
disp(Ic)