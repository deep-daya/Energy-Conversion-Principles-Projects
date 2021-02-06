function [Z, Rlam_min,  T_C,Eta_te, W,Q_H, Q_C, V_L, I_L ] = task_3_2(alpha, lam_A, lam_B, rho_A, rho_B, n, R_L, T_H, T_sat)
Rlam_min = ((lam_A*rho_A)^0.5 + (lam_B*rho_B)^0.5)^2;
Z= alpha^2 / Rlam_min;
A_pv = 0.1^2;
k_w = 6.5;
t_w = 0.4*10^-2;
U_ev = 25;
a = T_sat;
b = T_H+1;
m_max =@(T_C) (1 + 0.5*(T_H+T_C)*Z)^0.5;
R =@(m_max) R_L/(n * m_max);
lam = @(R) Rlam_min / R;
fQ_C =@(T_C)  A_pv*(k_w/t_w + U_ev) * (T_C-T_sat);
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
syms I_L V_L
eqn1=V_L- (n*alpha*(T_H-T_C))-R_batt_2*I_L==0;
eqn2= Q_H- n*alpha*I_L*T_H- lam_batt*(T_H-T_C)-0.5*((I_L)^2)*R_batt_2==0;
[V_L_sol, I_L_sol] = solve([eqn1,eqn2],V_L,I_L,'IgnoreAnalyticConstraints', true);
V_L_2= double(subs(V_L_sol))>0;
V_L = double(V_L_sol(double(subs(V_L_sol))>0));
I_L = double(I_L_sol(double(subs(I_L_sol))>0));

end 