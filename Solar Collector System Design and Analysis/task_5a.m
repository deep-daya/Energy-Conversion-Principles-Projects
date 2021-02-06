%Setting up parameters,constants
t_0 = 11;
t_i = 16;
T_i = 16;
d =1:14:365;
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
m_in = 0.0286;   % mini m kg/s 
zet = (200);
eps = 30;
T_out_ser = [];
T_out_3 = [];
T_out_2 = [];
T_avg_out = [];
comb_2 = [];
Id_ncol=[];
Id_ncol_1=[];
Id_ncol_2=[];
%Using series to find T_out_3
for i= 1:length(d)
    [T_out_ser(:,end+1),Id_ncol(:,end+1)]= task_4(t_0,t_i,d(i),zet,eps,T_i,m_in);
    [T_out_2(:,end+1),Id_ncol_1(:,end+1)]= task_4(t_0,t_i,d(i),zet,eps,T_out_ser,m_in);
    [T_out_3(:,end+1),Id_ncol_2(:,end+1)]= task_4(t_0,t_i,d(i),zet,eps,T_out_2,m_in);
end
%Finding the total energy and lost energy for the 27 days
for i = 1:length(Id_ncol)
    True_Energy(i)=(Id_ncol(i)+Id_ncol_1(i)+Id_ncol_2(i))*0.648;
    Lost_energy(i) =(Id_ncol(i)+Id_ncol_1(i)+Id_ncol_2(i))*0.352;
end
%Finding total energy lost for the year using trapezoidal rule
Total_lost_energy = 0;
for i =1:length(d)-1
Total_lost_energy = ((365-1)/27)*(Lost_energy(i)+Lost_energy(i+1))/2+Total_lost_energy;
end
%Finding the amount of natural gas amount burned at 90% eff., and T_Avg
Natural_gas = [];
for i = 1:length(T_out_3(1,:))
T_avg_out(i) = sum(T_out_3(:,i))/length(T_out_3);
Natural_gas(end+1)= (163.337*4186*(65-T_avg_out(i)))/(0.9*50050*1000);
end
Total_CH4 = 0;
%Finding total amount burned over year
for i =1:length(d)-1
Total_CH4 = ((365-1)/27)*((Natural_gas(i+1)+Natural_gas(i))/2)+Total_CH4;
end
%Finding total cost of methane burned
Total_cost = Total_CH4*9.45/(28.3168*.656)