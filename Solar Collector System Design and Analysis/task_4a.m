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
T_out_ser = [];
T_out_2 = [];
T_out_3 = [];
comb_2 = [];
%series arrangement
%Three nested for loops to determine the T_out of the third collector
for i = 1:length(m_range)
    x(i)=m_range(i);
    for j= 1:length(zet)
        y(j)=zet(j);
        for k=1:length(eps)
            z(k)=eps(k);
               [T_out_ser(:,end+1),Id_ncol] = task_4(t_0,t_i,d,y(j),z(k),T_i,x(i));
               [T_out_2(:,end+1),Id_ncol_2] = task_4(t_0,t_i,d,y(j),z(k),T_out_ser(:,i),x(i));
               [T_out_3(:,end+1),Id_ncol_3] = task_4(t_0,t_i,d,y(j),z(k),T_out_2(:,i),x(i));
                comb_2(:,end+1) = [y(j),z(k),x(i),Id_ncol];
        end;
    end;
end;
%finding the average temperatures
for i = 1:length(T_out_3(1,:))
T_avg_2(i) = sum(T_out_3(:,i))/length(T_out_3);
end;
%Finding temperatures between 64 and 66 and finding the closest one to 65
a= find(T_avg_2<66 & T_avg_2>64);
[m,l] = min(abs(T_avg_2(a)-65));
%Using all the inforamation to find the total mass, zeta, and epsilon of
%series configuration.
m_total_ser=(comb_2(3,a(l))*21600)/3.78541;
zeta_series = comb_2(1,a(l));
eps_series = comb_2(2,a(l));
%Plotting the Average temperature vs the 3 parameters changed 
figure;
subplot(2,2,1);
plot(comb_2(1,:),T_avg_2,'.k')
title('Average Temperature vs Zeta');
subplot(2,2,2);
plot(comb_2(2,:),T_avg_2, '.k');
title('Average Temperature vs Epsilon');
subplot(2,2,[3,4]);
plot(comb_2(3,:),T_avg_2,'.k');
title('Average Temperature vs Mass flow rate per collector');

%end

