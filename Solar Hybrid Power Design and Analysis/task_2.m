function [P, Eta,Q_pv, V_L,I_L] = task_2(R_L, Ly, Lz, Id,rc,T, Vg)  
%same as task_1
T_p=[293;303;313;333;353];
V_oc=[0.39;0.37;0.35;0.32;0.27];
T_pv = horzcat(T_p,ones(5,1));
co = T_pv\V_oc;
V_oc = [(T+273), 1]*co;
%constants
Kb= 1.38*10^-23;
qe= 1.6*10^-19;
L=V_oc*qe/(Kb*(T+273));
I0_Iv = 1/(exp(L)-1);
A_pv = (Ly*Lz)/10000; 
%I_L = V_L / R_L; 
%Eta = V_L * I_L / (p * A_pv);
% finding I_L
P_0= Id*rc*A_pv;

phi = (Id*rc*2.404*15)/(pi^4*Kb*(5778));

syms I_L

fun= @(x) (x.^2)./(exp(x)-1);

xmin = qe*Vg/(Kb* 5778);   

inte = integral(fun,xmin,Inf);

phi_g = phi * 0.416 * inte;

Iv = phi_g * qe * A_pv;

I0 =Iv * I0_Iv;
%solving for I_L
eqn= (Kb*(T+273)/qe)*log(((Iv- I_L)/I0)+1)- I_L*R_L==0;
I_L_sol = solve(eqn,I_L);
I_L = double(subs(I_L_sol));
%solving for V_L and then Efficiency
V_L = I_L*R_L;
Eta = V_L*I_L/ (Id*rc*A_pv);
P = I_L*V_L;
Q_pv = -V_L*I_L + P_0;
end 
