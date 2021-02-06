function [T2_f, W] =task_2a(T1,P1,P2, eta_comp, T_guess, gamma,n_prod,alpha)

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
    %T2 = T1* (1+((P1/P2)^((k-1)/k) -1)/eta_comp);
    T2= T1*(1+ eta_comp*((P1/P2)^((k-1)/k)-1));
    err= abs(T2-T2_g);
    if  err <= 0.2
        break 
    else 
        T2_g = T2; 
    end 
end
T2_f = T2_g;
H2= 0 - cp * (T2_f-T1)*n_prod;
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