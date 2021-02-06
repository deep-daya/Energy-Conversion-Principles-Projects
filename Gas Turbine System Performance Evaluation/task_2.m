function [T2_f, W] =task_2(T1,P1,P2, eta_comp)

%guess T2=T2_g to be average, Tf was taken from task 1 
Tf=1900; 
T2_g=Tf;
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



