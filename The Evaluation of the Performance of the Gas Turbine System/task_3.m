function [y_h20, y_co2, y_n2, y_o2, cp_prod, cp_air] = task_3(Tp, gamma,alpha)

%based on the equation given on page1 
n_ex = (3.76*alpha/4.76) + (alpha/4.76-3*gamma-2) + (2+2*gamma) + (1+2*gamma);
y_h20 = (2+2*gamma) / n_ex;
y_co2 = (1+2*gamma)/ n_ex;
y_n2 = ((3.76*alpha)/4.76)/ n_ex;
y_o2 = ((alpha/4.76)-3*gamma-2) / n_ex;

%cp of each product
cp_o2 = 25.48 + 1.520*10^(-2)*Tp - 0.7155*10^(-5)*(Tp)^2 + 1.312*10^(-9)*(Tp)^3;
cp_h20 = 32.24 + 0.1923*10^(-2)*Tp + 1.0551*10^(-5)*(Tp)^2 - 3.595*10^(-9)*(Tp)^3; 
cp_co2 = 22.26 + 5.981*10^(-2)*Tp - 3.501*10^(-5)*(Tp)^2 + 7.469*10^(-9)*(Tp)^3;
cp_n2 = 28.9 - 0.1571*10^(-2)*Tp + 0.8081*10^(-5)*(Tp)^2 - 2.873*10^(-9)*(Tp)^3;

%cp of whole product 
cp_prod = y_h20*cp_h20 + y_co2*cp_co2 + y_n2*cp_n2 + y_o2*cp_o2; 

cp_air = 28.11+ 0.1967*10^(-2)*( Tp) + 0.4802*10^(-5)*( Tp)^2 - 1.966*10^(-9)*( Tp)^3;
end
