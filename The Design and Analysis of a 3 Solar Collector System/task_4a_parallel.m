%Setting up parameters as in previous tasks.
t_0 = 10;
t_i = 16;
T_i = 16;
d =120;
Ac = 3.25;
alpha_c =0.85;
tau_g = 0.89;
h_convo = 7;
del_g = 0.007;
kg = 1.3;
h_convi = 3.1;
del_ins = .06;
k_ins = 0.045;
cp = 4186;
%Setting up data points for different m measurements.
m_min = 0.0267;   % mini m kg/s 
m_max= (350*3.78541)/21600;  % max m kg/s if the total flow is 350 gallons over 6 hrs
m_range= linspace(m_min,m_max,19);
%Setting up data points for Zeta, Epsilon for loops
zet = (0:20:360);
eps = (0:10:180);
T_out= [];
comb = [];
%parallel configuration
%Three nested for loops to determine the T_out of the third collector
for i = 1:length(m_range)
    x(i)=m_range(i)/3;
    for j= 1:length(zet)
        y(j)=zet(j);
        for k=1:length(eps)
            z(k)=eps(k);
            [T_out(:,end+1),Id_ncol] = task_4(t_0,t_i,d,y(j),z(k),T_i,x(i));
            comb(:,end+1) = [y(j),z(k),x(i),Id_ncol];
        end
    end
end
%finding the average temperatures
for i = 1:length(T_out(1,:))
T_avg(i) = sum(T_out(:,i))/length(T_out);
end;
%Finding the closest temp to 65C
I = find(T_avg<66 & T_avg>64);
[M,L] = min(abs(T_avg(I)-65));
%Finding total mass flow rate, Zeta, Epsilon for the parallel configuration
m_total=(3*comb(3,I(L))*21600)/3.78541;
zeta_parallel = comb(1,I(L));
eps_paralled = comb(2,I(L));
%Plotting the Average Temp vs the 3 parameters changed
subplot(2,2,1);
plot(comb(1,:),T_avg,'.k')
title('Average Temperature vs Zeta');
subplot(2,2,2);
plot(comb(2,:),T_avg, '.k');
title('Average Temperature vs Epsilon');
subplot(2,2,[3,4]);
plot(comb(3,:),T_avg,'.k');
title('Average Temperature vs Mass flow rate per collector');
