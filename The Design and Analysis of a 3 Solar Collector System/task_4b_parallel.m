%Setting up the parameters, constants
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
zet = (0:20:360);
eps = (0:10:180);
T_out= [];
comb = [];
%Finding the T_out for different zeta and epsilon configurations
    for j= 1:length(zet)
        y(j)=zet(j);
        for k=1:length(eps)
            z(k)=eps(k);
            [T_out(:,end+1),Id_ncol] = task_4(t_0,t_i,d,y(j),z(k),T_i,m_max);
            comb(:,end+1) = [y(j),z(k),Id_ncol];
        end
    end
%Finding average temperature
for i = 1:length(T_out(1,:))
T_avg(i) = sum(T_out(:,i))/length(T_out);
end;
%Finding closest avg temp to 65C
I = find(T_avg<66 & T_avg>64);
[M,L] = min(abs(T_avg(I)-65));
%Finding mass flow rate, zeta, eps
m_total=(3*comb(3,I(L))*21600)/3.78541;
zeta_parallel = comb(1,I(L));
eps_paralled = comb(2,I(L));
%returning max avg temp
max(T_avg);