%setting up parameters.
t_0 = 10;
t_i = 16;
T_i = 14;
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
m_max= (350*3.78541)/21600;  % max m kg/s if the total flow is 350 gallons over 6 hrs
t=linspace(t_0,t_i,216);
zet = (0:20:360);
eps = (0:10:180);
T_out_ser = [];
T_out_2 = [];
T_out_3 = [];
comb_2 = [];
T_avg_2=[];
%FInding the temp out of the third 
%series arrangement
    for j= 1:length(zet)
        y(j)=zet(j);
        for k=1:length(eps)
            z(k)=eps(k);
               [T_out_ser(:,end+1),Id_ncol] = task_4(t_0,t_i,d,y(j),z(k),T_i,m_max);
               [T_out_2(:,end+1),Id_ncol_2] = task_4(t_0,t_i,d,y(j),z(k),T_out_ser(:,end),m_max);
               [T_out_3(:,end+1),Id_ncol_3] = task_4(t_0,t_i,d,y(j),z(k),T_out_2(:,end),m_max);
                comb_2(:,end+1) = [y(j),z(k),Id_ncol];
        end;
    end;
for i = 1:length(T_out_3(1,:))
T_avg_2(i) = sum(T_out_3(:,i))/length(T_out_3);
end;
%Finding closest temp to 65, mass flow rate and zeta, and epsilon.
a= find(T_avg_2<66 & T_avg_2>64);
[m,l] = min(abs(T_avg_2(a)-65));
m_total_ser=(comb_2(3,a(l))*21600)/3.78541;
zeta_series = comb_2(1,a(l));
eps_series = comb_2(2,a(l));
%returning max temp
max(T_avg_2);

