function [Power, Heat, T5_f, alpha_f, W_c, W_t, n_comp,Work_Class, Heat_Class, n_class]= task_4(~)
P1= 101;
T1=298;
P2=500;
T4=1600;
gamma= .25;
alpha_stoich = 4.76*(2+3*gamma);
Qs=0;
eff_comp=0.85;
eff_turb=0.85;
eps_regen=0.75;
i = 1;
T5_f(1,:) =1000;
alpha(1,:) = alpha_stoich;
while i<4
    m_prod(i,:) = 28.014*3.76*alpha(i,:)/4.76 + 15.9999*2 *(alpha(i,:)/4.76 - 3*gamma -2)+18.01528*(2+2*gamma)+44.01*(1+2*gamma);
    m_air = 6;
    n_dot_air(i,:)= m_air/28.97;
    n_prod(i,:) = 3.76*alpha(i,:)/4.76 + (alpha(i,:)/4.76 - 3*gamma -2)+(2+2*gamma)+(1+2*gamma);
    n_dot_prod(i,:)= m_air/(m_prod(i,:)/n_prod(i,:));
    [T2_final(i,:),W(i,:),cp_2(i,:)] = task_2(T1,P1,P2,eff_comp,1900);
    [T5_f(i+1,:), W_5(i,:),cp_5(i,:)] =task_2a(T4,P1,P2, eff_turb,T5_f(i,:), gamma,alpha(i,:));
    T_avg_25= 0.5*(T2_final(i,:)+T5_f(i,:));
    [y_h20(1,:), y_co2(i,:), y_n2(i,:), y_o2(i,:), cp_prod(i,:), cp_air(i,:)] = task_3(T_avg_25, gamma,alpha(i,:));
    ncp_min = min(cp_prod(i,:)*n_dot_prod(i,:), cp_air(i,:)*n_dot_prod(i,:));
    T2_r(i,:) = T2_final(i,:)+ (ncp_min* eps_regen* (T5_f(i,:)- T2_final(i,:)))/(n_dot_air(i,:)*cp_air(i,:));
    [T3, cp_23] = root_finder(Qs,T2_r(i,:),n_dot_air(i,:)); 
    i = i+1;
    alpha(i,:)= task_1(T4, T3, gamma);
end

alpha_f = alpha(3,:);
W_t = n_dot_prod(3,:)*cp_5(3,:)*(T4-T5_f(3,:));
W_c = n_dot_air(3,:)*cp_2(3,:)*(T2_final(3,:)-T1);
Power= n_dot_prod(3,:)*cp_5(3,:)*(T4-T5_f(3,:))- n_dot_air(3,:)*cp_2(3,:)*(T2_final(3,:)-T1);
[~, ~, ~, ~, cp_burn, ~] = task_3(.5*(T3+T4), gamma,alpha(3,:));
Heat =(n_dot_prod(3,:)*cp_burn*(T4-T3))+ Qs;

k=(28.11/(28.11-8.314)+ (28.3369/(28.3369-8.314)))/2;

Work_Class =(eff_turb*(1-(P1/P2)^((k-1)/k))- (1/eff_comp)*(T1/T4)*((P2/P1)^((k-1)/k)-1));
Heat_Class = (1-eps_regen*(1-eff_turb+eff_turb*((P1/P2)^((k-1)/k)))-(1-eps_regen)*(T1/T4)*(1-(1/eff_comp)+(1/eff_comp)*(P2/P1)^((k-1/k))));
n_comp = (n_dot_prod(3,:)*cp_5(3,:)*(T4-T5_f(3,:))- n_dot_air(3,:)*cp_2(3,:)*(T2_final(3,:)-T1))/(n_dot_prod(3,:)*cp_burn*(T4-T3));
n_class = (eff_turb*(1-(P1/P2)^((k-1)/k))- (1/eff_comp)*(T1/T4)*((P2/P1)^((k-1)/k)-1))/(1-eps_regen*(1-eff_turb+eff_turb*((P1/P2)^((k-1)/k)))-(1-eps_regen)*(T1/T4)*(1-(1/eff_comp)+(1/eff_comp)*(P2/P1)^((k-1/k))));
end
function [T2_f, W,cp] =task_2(T1,P1,P2, eta_comp, T_guess)

%guess T2=T2_g to be average, Tf was taken from task 1 
T2_g=T_guess;
R=8.314;

err=1; 
while err >0.2
    T_avg = 0.5*(T1+T2_g);
    cp = 28.11+ 0.1967*10^(-2)*(T_avg) + 0.4802*10^(-5)*(T_avg)^2 - 1.966*10^(-9)*(T_avg)^3;
    cv = cp-R;
    k= cp/cv;
    %first T2 to calculate difference
    T2 = T1* (1+((P2/P1)^((k-1)/k) -1)/eta_comp);
    %try to get second T2 for difference calculation
    T_avg_1 =  0.5*(T1+T2); 
    cp_1 = 28.11+ 0.1967*10^(-2)*( T_avg_1) + 0.4802*10^(-5)*( T_avg_1)^2 - 1.966*10^(-9)*( T_avg_1)^3;
    cv_1 = cp_1-R;
    k_1= cp_1/cv_1;
    T2_1 = T1* (1+((P2/P1)^((k_1-1)/k_1) -1)/eta_comp);
    err= abs(T2-T2_1);  
    if  err <= 0.2
        break 
    else 
        T2_g = T2; 
    end 
end
T2_f = T2_g;
H2= 0 + cp * (T2_f-T1)/(28.97);
%W=H2-H1, and H1=0
W=H2;

end
function [T2_f, W,cp] =task_2a(T1,P1,P2, eta_comp, T_guess, gamma,alpha)
%guess T2=T2_g to be average, Tf was taken from task 1 
T2_g=T_guess;
R=8.314;

err=1; 
while err >0.2
    T_avg = 0.5*(T1+T2_g);
    [~, ~, ~, ~, cp_prod, ~] = task_3(T_avg, gamma,alpha);
    cp = cp_prod;
    cv = cp-R;
    k= cp/cv;
    %first T2 to calculate difference
    T2 = T1* (1+((P1/P2)^((k-1)/k) -1)*eta_comp);
    err= abs(T2-T2_g);  
    if  err <= 0.2
        break 
    else 
        T2_g = T2; 
    end 
end
T2_f = T2_g;
H2= 0 + cp * (T1-T2_f)/(28.97);
%W=H2-H1, and H1=0
W=H2;

end


function [y_h20, y_co2, y_n2, y_o2, cp_prod, cp_air] = task_3(Tp, gamma,alpha)

%based on the equation given on page1 
n_ex = (3.76*alpha/4.76) + (alpha/4.76-3*gamma-2) + (2+2*gamma) + (1+2*gamma);
y_h20 = (2+2*gamma) / n_ex;
y_co2 = (1+2*gamma)/ n_ex;
y_n2 = (3.76*alpha/4.76)/ n_ex;
y_o2 = (alpha/4.76-3*gamma-2) / n_ex;

%cp of each product
cp_o2 = 25.48 + 1.520*10^(-2)*Tp - 0.7155*10^(-5)*(Tp)^2 + 1.312*10^(-9)*(Tp)^3;
cp_h20 = 32.24 + 0.1923*10^(-2)*Tp + 1.0551*10^(-5)*(Tp)^2 - 3.595*10^(-9)*(Tp)^3; 
cp_co2 = 22.26 + 5.981*10^(-2)*Tp - 3.501*10^(-5)*(Tp)^2 + 7.469*10^(-9)*(Tp)^3;
cp_n2 = 28.9 - 0.1571*10^(-2)*Tp + 0.8081*10^(-5)*(Tp)^2 - 2.873*10^(-9)*(Tp)^3;

%cp of whole product 
cp_prod = y_h20*cp_h20 + y_co2*cp_co2 + y_n2*cp_n2 + y_o2*cp_o2; 
cp_air = 28.11+ 0.1967*10^(-2)*( Tp) + 0.4802*10^(-5)*( Tp)^2 - 1.966*10^(-9)*( Tp)^3;
end

function [alpha] = task_1(Tp, Tr, gamma)


%the heat of formation from table, [kJ/kmol]
h0_c3h8 =-103850;
h0_ch4 = -74850; 
h0_o2 = 0; 
h0_h2o = -241820; 
h0_co2 = -393520; 
h0_n2 = 0;

%specific heat from table, [kJ/kmol*K]
T_avg = 0.5 * (Tr+Tp);
cp_c3h8 = -4.04 + 30.48*10^(-2)*T_avg - 15.72*10^(-5)*(T_avg)^2 + 31.74*10^(-9)*(T_avg)^3;
cp_ch4 = 19.89 + 5.024*10^(-2)*T_avg +1.269*10^(-5)*(T_avg)^2 - 11.01*10^(-9)*(T_avg)^3; 
cp_o2 = 25.48 + 1.520*10^(-2)*T_avg - 0.7155*10^(-5)*(T_avg)^2 + 1.312*10^(-9)*(T_avg)^3;
cp_h2o = 32.24 + 0.1923*10^(-2)*T_avg + 1.055*10^(-5)*(T_avg)^2 - 3.595*10^(-9)*(T_avg)^3; 
cp_co2 = 22.26 + 5.981*10^(-2)*T_avg - 3.501*10^(-5)*(T_avg)^2 + 7.469*10^(-9)*(T_avg)^3;
cp_n2 = 28.9 - 0.1571*10^(-2)*T_avg + 0.8081*10^(-5)*(T_avg)^2 - 2.873*10^(-9)*(T_avg)^3;

%enthalpy 
h_c3h8 =h0_c3h8 + cp_c3h8 *(Tr-298.15);
h_ch4 = h0_ch4 + cp_ch4  * (Tr-298.15); 
h_o2_i = h0_o2 + cp_o2 * (Tr-298.15);
h_o2_o = h0_o2 + cp_o2 * (Tp-298.15);
h_h2o = h0_h2o + cp_h2o * (Tp-298.15); 
h_co2 = h0_co2 + cp_co2 * (Tp-298.15);
h_n2_i = h0_n2 + cp_n2 * (Tr-298.15);
h_n2_o = h0_n2 + cp_n2 * (Tp-298.15);

alpha = 4.76*(-gamma*h_c3h8 - (1-gamma)*h_ch4 - (3*gamma+2)*h_o2_o + (2+2*gamma)*h_h2o + (1+2*gamma)*h_co2)/(h_o2_i + 3.76*h_n2_i - 3.76*h_n2_o - h_o2_o);

end 
function [T3_f, cp23] = root_finder(Qs, T2_r, n_dot_air)
syms cp_23 T3
eqn = [cp_23 == 28.11+ 0.1967*10^(-2)*(0.5*(T3+T2_r)) + 0.4802*10^(-5)*(0.5*(T3+T2_r))^2 - 1.966*10^(-9)*(0.5*(T3+T2_r))^3,T3 == Qs/(n_dot_air*cp_23)+T2_r];
v = solve(eqn, [cp_23, T3]);
cp_23= vpa(v.cp_23);
T3 = vpa(v.T3);
p = find(T3 >800& T3 < 1500);
T3_f = double(T3(p));
cp23 = double(cp_23(p));
end

