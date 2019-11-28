rc = 15;
Ly = 10 ;
Lz = 10 ;
Vg = 1.1;
fg = 265 ;
Rs = 0;
Id_i = linspace(720,1080,12);
R_L_pv = 0.0070;
T_sat = 90.2 ;
errr = 1;
alpha= 0.0017;
lam_A= 0.032; 
lam_B= 0.021;
rho_A=0.0020;
rho_B=0.0030;
n=12;
R_L=0.1; 
i =1;
a1 = 150;
b1= 250;
lamb_cont= 300;
for l = 1:length(Id_i)
    Id = Id_i(l);
    disp(l);
    i=1;
    errr = 1;
    a1 = T_sat+1;
    b1= 383;
while errr > 10^(-5)
c1 = (a1+b1)/2; 
T_pv_i(i,:)= [a1,c1,b1];
[P(i,1), Eta(i,1),Q_pv(i,1), V_L(i,1),I_L(i,1)]= task_2(R_L_pv, Ly,Lz, Id, rc,a1,Vg);
[P(i,2), Eta(i,2),Q_pv(i,2), V_L(i,2),I_L(i,2)]= task_2(R_L_pv, Ly,Lz, Id, rc,c1,Vg);
[P(i,3), Eta(i,3),Q_pv(i,3), V_L(i,3),I_L(i,3)]= task_2(R_L_pv, Ly,Lz, Id, rc,b1,Vg);
T_H_i(i,:) = T_pv_i(i,:) - Q_pv(i,:)./lamb_cont;
[Z(i,1), Rlam_min(i,1),  T_C(i,1),Eta_tee(i,1), W(i,1),Q_H(i,1), Q_C(i,1), V_L_te(i,1), I_L_te(i,1) ] = task_3_2(alpha, lam_A, lam_B, rho_A, rho_B, n, R_L, T_H_i(i,1), T_sat);
[Z(i,2), Rlam_min(i,2),  T_C(i,2),Eta_tee(i,2), W(i,2),Q_H(i,2), Q_C(i,2), V_L_te(i,2), I_L_te(i,2) ] = task_3_2(alpha, lam_A, lam_B, rho_A, rho_B, n, R_L, T_H_i(i,2), T_sat);
[Z(i,3), Rlam_min(i,3),  T_C(i,3),Eta_tee(i,3), W(i,3),Q_H(i,3), Q_C(i,3), V_L_te(i,3), I_L_te(i,3) ] = task_3_2(alpha, lam_A, lam_B, rho_A, rho_B, n, R_L, T_H_i(i,3), T_sat);
E_Q(i,:) = Q_pv(i,:) - Q_H(i,:); 
if E_Q(i,1)>0 &&  E_Q(i,2) < 0 
        b1=c1;
    else 
        a1=c1;
    end 
errr = abs(E_Q(i,2));
i = i+1;
end
T_H(1,l) = T_H_i(i-1,2); 
T_C(1,l)= T_C(i-1,2);
T_pv(1,l) = T_pv_i(i-1,2);
Eta_te(1,l) = Eta_tee(i-1,2);
Eta_pv(1,l)= Eta(i-1,2);
P_pv(1,l) = P(i-1,2);
W_te(1,l)= W(i-1,2);
V_L(1,l)= V_L(i-1,2);
I_L(1,l) = I_L(i-1,2);
P_tot(1,l) = P_pv(1,l)+W_te(1,l);
Eta_tot(1,l) = P_tot(1,l)/ (Id*rc*(0.1^2));
V_L_te(1,l) = V_L_te(i-1,2);
I_L_te(1,l) = I_L_te(i-1,2);
end
subplot(2,2,1)
plot(Id_i,transpose(Eta_tot),'.');
title('Subplot 1: Total Efficiency vs Variations in Incident Radiation');
xlabel('Incident Radiation (W/(m^2))')
subplot(2,2,2)
plot(Id_i,transpose(Eta_te),'.');
title('Subplot 2: Thermal Efficiency vs Variations in Incident Radiation');
subplot(2,2,3)
plot(Id_i,transpose(Eta_pv),'.');
title('Subplot 3: PV Efficiency vs Variations in Incident Radiation');
subplot(2,2,4)
plot(Id_i,transpose(P_tot), '.');

title('Subplot 4: Total Power vs Variations in Incident Radiation');
function [P, Eta,Q_pv, V_L,I_L] = task_2(R_L, Ly, Lz, Id,rc,T, Vg)  
T_p=[293;303;313;333;353];
V_oc=[0.39;0.37;0.35;0.32;0.27];
T_pv = horzcat(T_p,ones(5,1));
co = T_pv\V_oc;
V_oc = [(T), 1]*co;

%constants
Kb= 1.38*10^-23;
qe= 1.6*10^-19;
L=V_oc*qe/(Kb*(T));
I0_Iv = 1/(exp(L)-1);
A_pv = (Ly*Lz)/10000; 
%I_L = V_L / R_L; 
%Eta = V_L * I_L / (p * A_pv);
P_0= Id*rc*A_pv;

phi = (Id*rc*2.404*15)/(pi^4*Kb*(5778));

syms I_L

fun= @(x) (x.^2)./(exp(x)-1);

xmin = qe*Vg/(Kb* 5778);   

inte = integral(fun,xmin,Inf);

phi_g = phi * 0.416 * inte;

Iv = phi_g * qe * A_pv;

I0 =Iv * I0_Iv;

eqn= (Kb*(T)/qe)*log(((Iv- I_L)/I0)+1)- I_L*R_L==0;
I_L_sol = solve(eqn,I_L);
I_L = double(subs(I_L_sol));
V_L = I_L*R_L;
Eta = V_L*I_L/ (Id*rc*A_pv);
P = I_L*V_L;
Q_pv = -V_L*I_L + P_0;
end 
 
function [Z, Rlam_min,  T_C,Eta_te, W,Q_H, Q_C, V_L_te, I_L_te] = task_3_2(alpha, lam_A, lam_B, rho_A, rho_B, n, R_L, T_H, T_sat)
%%%%%add I_l and v_l for this question 
Rlam_min = ((lam_A*rho_A)^0.5 + (lam_B*rho_B)^0.5)^2;
Z= alpha^2 / Rlam_min;
A_pv = 0.1^2;
k_w = 6.5;
t_w = 0.4*10^-2;
U_ev = 25;

% guess T_C = (T_H + T_sat)/2 = (323+95)/2 = 209 K 
a = T_sat;
b = T_H+1;
m_max =@(T_C) (1 + 0.5*(T_H+T_C)*Z)^0.5;
R =@(m_max) R_L/(n * m_max);
lam = @(R) Rlam_min / R;
fQ_C =@(T_C) A_pv* (k_w/t_w + U_ev) * (T_C-T_sat);
Eta =@(T_C, m_max) ((T_H-T_C) / T_H) *(((1+m_max)^2 / m_max) * (Z*T_H)^(-1) + 1 + (1/(2*m_max))*(1+T_C/T_H))^(-1);
R_batt = @(R) n* R;
fW =@(T_C, R_batt) (n^2 * alpha^2 * (T_H-T_C)^2 * R_L )/ ((R_L + R_batt)^2);
fQ_H =@(W,Eta) W/Eta;
Q_CTE = @(Q_H, W)Q_H - W; 
E_QC =@(Q_CTE, Q_C) Q_CTE - Q_C; 
err = 1;
while err > 10^(-5)
   c=(a+b)/2; 
    m_max_a =m_max(a);
    m_max_b =m_max(b);
    m_max_c =m_max(c);
    R_a =R(m_max_a); 
    R_b =R(m_max_b); 
    R_c =R(m_max_c);
    lam_a = lam(R_a); 
    lam_b = lam(R_b); 
    lam_c = lam(R_c);
    Q_C_a =fQ_C(a);
    Q_C_b =fQ_C(b);
    Q_C_c =fQ_C(c);
    Eta_a =Eta(a, m_max_a);
    Eta_b =Eta(b, m_max_b);
    Eta_c =Eta(c, m_max_c);
    R_batt_a = R_batt(R_a);
    R_batt_b = R_batt(R_b);
    R_batt_c = R_batt(R_c);
    W_a =fW(a, R_batt_a);
    W_b =fW(b, R_batt_b);
    W_c =fW(c, R_batt_c);
    Q_H_a =fQ_H(W_a,Eta_a);
    Q_H_b =fQ_H(W_b,Eta_b);
    Q_H_c =fQ_H(W_c,Eta_c);
    Q_CTE_a = Q_CTE(Q_H_a, W_a);
    Q_CTE_b = Q_CTE(Q_H_b, W_b);
    Q_CTE_c = Q_CTE(Q_H_c, W_c);
    E_QC_a =E_QC(Q_CTE_a, Q_C_a);
    E_QC_b =E_QC(Q_CTE_b, Q_C_b);
    E_QC_c =E_QC(Q_CTE_c, Q_C_c); 
    if E_QC_a>0 &&  E_QC_c < 0 
        b=c;
    else 
        a=c;
    end 
    err = abs(E_QC_c);
end 

T_C = c;
W = W_c;
Q_H = Q_H_c;
Q_C = Q_C_c;
lam_batt = n*lam_c;
R_batt_2 = R_batt_c;
Eta_te = Eta_c;
syms I_L_te V_L_te
eqn1=V_L_te- (n*alpha*(T_H-T_C))-R_batt_2*I_L_te==0;
eqn2= Q_H- n*alpha*I_L_te*T_H- lam_batt*(T_H-T_C)-0.5*((I_L_te)^2)*R_batt_2==0;
%%%%%need to add task_2 
[V_L_sol, I_L_sol] = solve([eqn1,eqn2],V_L_te,I_L_te,'IgnoreAnalyticConstraints', true);
V_L_2= double(subs(V_L_sol))>0;
V_L_te = double(V_L_sol(double(subs(V_L_sol))>0));
I_L_te = double(I_L_sol(double(subs(I_L_sol))>0));

end 