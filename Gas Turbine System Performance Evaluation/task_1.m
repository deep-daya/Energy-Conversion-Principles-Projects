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
