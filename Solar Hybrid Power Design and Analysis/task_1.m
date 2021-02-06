T = 30;
%linear curve fit to find the relationship b/t T and V_oc, 5 points are
%estimated based on the graph provided 

T_p=[293;303;313;333;353];
V_oc=[0.39;0.37;0.35;0.32;0.27];
T_pv = horzcat(T_p,ones(5,1));
co = T_pv\V_oc;
V_L = [(T+273), 1]*co;

%constants
Kb= 1.38*10^-23;
qe= 1.6*10^-19;
x=V_L*qe/(Kb*(T+273));
I0_Iv = 1/(exp(x)-1);
I0_Iv
