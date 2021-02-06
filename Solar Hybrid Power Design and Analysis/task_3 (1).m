%function [Z, Rlam_min,  T_C, W,Q_H, Q_C ] = task_3(alpha, lam_A, lam_B, rho_A, rho_B, n, R_L, T_H, T_sat)
%%%%%add I_l and v_l for this question 
alpha= 0.0017;
lam_A= 0.032; 
lam_B= 0.021;
rho_A=0.0020;
rho_B=0.0030;
n=12;
R_L=0.1; 
T_H= 323; 
T_sat= 95;
%part a and b
Rlam_min = ((lam_A+rho_A)^0.5 + (lam_B+rho_B)^0.5)^2;
Z= alpha^2 / Rlam_min;
A_pv = .10^2;
k_w = 6.5;
t_w = 0.4;
U_ev = 25;

% guess T_C = (T_H + T_sat)/2 = (323+95)/2 = 209 K 
a = T_sat;
b = T_H+1;
m_max =@(T_C) (1 + 0.5*(T_H+T_C)*Z)^0.5;
R =@(m_max) R_L / (n * m_max);
lam = @(R) Rlam_min / R;
fQ_C =@(T_C) A_pv * (k_w/t_w + U_ev) * (T_C-T_sat);
Eta =@(T_C, m_max) ((T_H-T_C) / T_H) *((((1+m_max)^2 / m_max) * ((Z*T_H)^(-1))) + 1 + ((1/(2*m_max))*(1+(T_C/T_H))))^(-1);
R_batt = @(R) n* R;
fW =@(T_C, R_batt) (n^2 * alpha^2 * (T_H-T_C)^2 * R_L )/ ((R_L + R_batt)^2);
fQ_H =@(W,Eta) W/Eta;
Q_CTE = @(Q_H, W)Q_H - W; 
E_QC =@(Q_CTE, Q_C) Q_CTE - Q_C; 
err = 1;
while err > 10^(-5)
    %if a1 == a && b1 == b
        c=(a+b)/2; 
        m_max_mat(:)= [m_max(a), m_max(b), m_max(c)];
        R_mat(:) = [R(m_max_mat(1)),R(m_max_mat(2)),R(m_max_mat(3))];
        lam_mat(:) = [lam(R_mat(1)),lam(R_mat(2)),lam(R_mat(3))];
        Q_C_mat(:) = [fQ_C(a),fQ_C(b),fQ_C(c)]; 
        Eta_mat(:)= [Eta(a,m_max_mat(1)),Eta(b,m_max_mat(2)),Eta(c,m_max_mat(3))];
        R_batt_mat(:) = [R_batt(R(1)),R_batt(R(2)),R_batt(R(3))];
        W_mat(:) = [fW(a, R_batt_mat(1)),fW(b, R_batt_mat(2)),fW(c, R_batt_mat(3))];
        Q_H_mat(:)=[fQ_H(W_mat(1),Eta_mat(1)),fQ_H(W_mat(2),Eta_mat(2)),fQ_H(W_mat(3),Eta_mat(3))]; 
        Q_CTE_mat(:) = [Q_CTE(Q_H_mat(1), W_mat(1)),Q_CTE(Q_H_mat(2), W_mat(2)),Q_CTE(Q_H_mat(3), W_mat(3))];
        E_QC_mat(:)= [E_QC(Q_CTE_mat(1), Q_C_mat(1)),E_QC(Q_CTE_mat(2), Q_C_mat(2)),E_QC(Q_CTE_mat(3), Q_C_mat(3))];
        if E_QC_mat(1) * E_QC_mat(3) < 0
            b=c;
        else 
            a=c;
        end 
        err = abs(E_QC_mat(3));
    %else
        
end 

T_C = c;
W = W_c;
Q_H = Q_H_c;
Q_C = Q_C_c;

%%%%%need to add task_2 

%end 