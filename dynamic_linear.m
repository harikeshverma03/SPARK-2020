function[dyn,Alloc_linear,obj_dyn_mu,obj_dyn_sd,obj_lin_mu,obj_lin_sd,dem_mean,sol_mean] = dynamic_linear(n,Beta)
Cost_imp_mu = 12;
Cost_imp_sd = 4;
 Rain_mu = csvread('rainfall_avg.csv')';
 Rain_sd = csvread('rain_sd.csv')';
 
 V1_mu = [627.3 656.99 598.6];
 MDDL = [619.4 646 590.1];
 V1_sd = (V1_mu - MDDL)*0.233826;
 
 E_mu = csvread('solar.csv');
 E_sd = csvread('solar_sd.csv');
[dyn,sum,obj_dyn] = run_dyn(n,Cost_imp_mu,Beta,Rain_mu,Rain_sd,V1_mu,V1_sd,E_mu,E_sd);
obj_dyn_mu = mean(obj_dyn);
obj_dyn_sd = std(obj_dyn);


Day = csvread('Synthetic_Days.csv');
Demand_mu(:,:,1) = xlsread('Demand Data.csv', 'A2:B25');
Demand_mu(:,:,2) = xlsread('Demand Data.csv', 'C2:D25');
Demand_mu(:,:,3) = xlsread('Demand Data.csv', 'E2:F25');
Demand_mu(:,:,4) = xlsread('Demand Data.csv', 'G2:H25');
Demand_mu(:,:,5) = xlsread('Demand Data.csv', 'I2:J25');
Demand_mu(:,:,6) = xlsread('Demand Data.csv', 'K2:L25');
Demand_mu(:,:,7) = xlsread('Demand Data.csv', 'M2:N25');
Demand_mu(:,:,8) = xlsread('Demand Data.csv', 'O2:P25');
Demand_mu(:,:,9) = xlsread('Demand Data.csv', 'Q2:R25');
Demand_mu(:,:,10) = xlsread('Demand Data.csv', 'S2:T25');
Demand_mu(:,:,11) = xlsread('Demand Data.csv', 'U2:V25');
Demand_mu(:,:,12) = xlsread('Demand Data.csv', 'W2:X25');

Demand_sd(:,:,1) = xlsread('Demand Data sd.csv', 'A2:B25');
Demand_sd(:,:,2) = xlsread('Demand Data sd.csv', 'C2:D25');
Demand_sd(:,:,3) = xlsread('Demand Data sd.csv', 'E2:F25');
Demand_sd(:,:,4) = xlsread('Demand Data sd.csv', 'G2:H25');
Demand_sd(:,:,5) = xlsread('Demand Data sd.csv', 'I2:J25');
Demand_sd(:,:,6) = xlsread('Demand Data sd.csv', 'K2:L25');
Demand_sd(:,:,7) = xlsread('Demand Data sd.csv', 'M2:N25');
Demand_sd(:,:,8) = xlsread('Demand Data sd.csv', 'O2:P25');
Demand_sd(:,:,9) = xlsread('Demand Data sd.csv', 'Q2:R25');
Demand_sd(:,:,10) = xlsread('Demand Data sd.csv', 'S2:T25');
Demand_sd(:,:,11) = xlsread('Demand Data sd.csv', 'U2:V25');
Demand_sd(:,:,12) = xlsread('Demand Data sd.csv', 'W2:X25');

Sol_mu(:,:,1) = repmat(xlsread('Solar_mu_daily.csv','A1:A24'),1,2);
Sol_mu(:,:,2) = repmat(xlsread('Solar_mu_daily.csv','B1:B24'),1,2);
Sol_mu(:,:,3) = repmat(xlsread('Solar_mu_daily.csv','C1:C24'),1,2);
Sol_mu(:,:,4) = repmat(xlsread('Solar_mu_daily.csv','D1:D24'),1,2);
Sol_mu(:,:,5) = repmat(xlsread('Solar_mu_daily.csv','E1:E24'),1,2);
Sol_mu(:,:,6) = repmat(xlsread('Solar_mu_daily.csv','F1:F24'),1,2);
Sol_mu(:,:,7) = repmat(xlsread('Solar_mu_daily.csv','G1:G24'),1,2);
Sol_mu(:,:,8) = repmat(xlsread('Solar_mu_daily.csv','H1:H24'),1,2);
Sol_mu(:,:,9) = repmat(xlsread('Solar_mu_daily.csv','I1:I24'),1,2);
Sol_mu(:,:,10) = repmat(xlsread('Solar_mu_daily.csv','J1:J24'),1,2);
Sol_mu(:,:,11) = repmat(xlsread('Solar_mu_daily.csv','K1:K24'),1,2);
Sol_mu(:,:,12) = repmat(xlsread('Solar_mu_daily.csv','L1:L24'),1,2);

Sol_sd(:,:,1)= repmat(xlsread('Solar_sd_daily.csv','A1:A24'),1,2);
Sol_sd(:,:,2) = repmat(xlsread('Solar_sd_daily.csv','B1:B24'),1,2);
Sol_sd(:,:,3) = repmat(xlsread('Solar_sd_daily.csv','C1:C24'),1,2);
Sol_sd(:,:,4) = repmat(xlsread('Solar_sd_daily.csv','D1:D24'),1,2);
Sol_sd(:,:,5) = repmat(xlsread('Solar_sd_daily.csv','E1:E24'),1,2);
Sol_sd(:,:,6) = repmat(xlsread('Solar_sd_daily.csv','F1:F24'),1,2);
Sol_sd(:,:,7) = repmat(xlsread('Solar_sd_daily.csv','G1:G24'),1,2);
Sol_sd(:,:,8) = repmat(xlsread('Solar_sd_daily.csv','H1:H24'),1,2);
Sol_sd(:,:,9) = repmat(xlsread('Solar_sd_daily.csv','I1:I24'),1,2);
Sol_sd(:,:,10) = repmat(xlsread('Solar_sd_daily.csv','J1:J24'),1,2);
Sol_sd(:,:,11) = repmat(xlsread('Solar_sd_daily.csv','K1:K24'),1,2);
Sol_sd(:,:,12) = repmat(xlsread('Solar_sd_daily.csv','L1:L24'),1,2);

reg = zeros(12,24);
hol = zeros(12,24);
obj_lin_mu = zeros(1,12);
obj_lin_sd = zeros(1,12);
dem_mean = zeros(24,2,12);
sol_mean =zeros(24,2,12);
Alloc_linear = zeros(25,17,12);
for i = 1:12
    
    Hydro_avg = sum(i);
   [Alloc_linear(:,:,i),reg(i,:),hol(i,:),obj_lin_mu(i),obj_lin_sd(i),dem_mean(:,:,i),sol_mean(:,:,i)] = run_linear(Day(i,:),Beta,Hydro_avg,Demand_mu(:,:,i),Demand_sd(:,:,i),Cost_imp_mu,Cost_imp_sd,Sol_mu(:,:,i),Sol_sd(:,:,i),n);
   %[Hyd_reg(i,:),Hyd_hol(i,:)] = func.Alloc_Hydro(Alloc_linear(:,:,i),24);
end

x1 = 1:24;
subplot(2,2,[1,2]);
x = 1:12;
plot(x,sum);
subplot(2,2,3); 
    plot(x1,reg(i,:));
    subplot(2,2,4); 
    plot(x1,hol(i,:));
obj_dyn = [obj_dyn_mu,obj_dyn_sd];
obj_lin(1,:) = obj_lin_mu;
obj_lin(2,:) = obj_lin_sd;
csvwrite('Hydro_Alloc_dyn.csv',dyn);
csvwrite('Dynamic_objective.csv',obj_dyn);
csvwrite('Monthly_Alloc.csv',Alloc_linear);
csvwrite('Linear_objective.csv',obj_lin');
csvwrite('Demand_mean.csv',dem_mean);
csvwrite('Solar_mean.csv',sol_mean);

end