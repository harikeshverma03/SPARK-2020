function [dyn,lin] = run()  
%day = 1*2 matrix for no. of holidays and regular days
%cost is a 1*n matrix for per power plant
%demand is a t*2 matrix for demand in holiday and typical day
%cost_imp is a 1*1 matrix for cost of import
%Hydro_cap is a 1*1 matrix for hydro capacity of the hour
%Thermal_cap = 1*n matrix for n power plants
%Hydro_avg 1*1 matrix for hydro capacity of the month
%Ramp_Rate = 1*n matrix for hourly ramp rate of n power plants
%Sol_E = t*2 matrix for solar radiance per hour in hoiday and normal day
%Beta = fraction of solar PV array
%Dem = 1*t matrix for demand in t stages
%Cost = 1*t matrix for cost in t stages
%Sol_E = 1*t matrix for solar radiance in t stages
%Beta = fraction of solar PV array
%Vol_1 = 1*1 matrix for the volume at the start of 1st stage
%Rain = t*p matrix for the incoming water due to rainfall - spillway discharge for p plant in t stage
%MDDL = 1*p matrix for Min. draw down level of p plant
%FRL = 1*p matrix for full reservoir level of p plant
%t = No. of stages for which we want to solve
%p = no. of hydro plants available
%Avg = Average generation in last five years

   Cost = [1 2 3 3 2 4 2 3 2 3 4 2];
   Cost_imp = 6;
   Day = [2 2];
   Demand = [10 15; 11 17; 9 14];
   Hydro_cap = 3000;
   Thermal_cap = [5000 7000 3000];
   Hydro_avg = 30;
   Ramp_Rate = [1 1.2 0.8];
   Sol_E = [1 1.2;1.2 1; 0.9 0.9]*4;
   Cost_1 = [2 5 2.5]; 
     Dem = 300:.1:301.1;
     E = ones(1,12)*4;
     V1 = [60 60 60];
     Rain = ones(12,3)*20;
     Rain(7:10,:) = Rain(7:10,:)+Rain(7:10,:)+Rain(7:10,:)*0.5;
     p = 3;t=3;n=3;t1=12;
     Avg =21600;
     Beta = 0.2;
     FRL = [100 110 90];
     MDDL = [30 10 40];
     Hp = [30 40 120];
     z = func.dynamic_stage(Cost,Dem,E,Beta,V1,Rain,MDDL,FRL,Hp,t1,p,Avg);
     y = func.linear_stage(Day, Cost_1,Cost_imp,Demand,Hydro_cap, Thermal_cap,Hydro_avg,Ramp_Rate,Sol_E,n,t,Beta);
     %disp(z);
     
%dynamic allocation t1 * p
res_dynamic = z.x';
Alloc_dynamic = zeros(t1+1,p+1);
for i = 1:p
    Alloc_dynamic(1,1+i) = i;
end
for i = 1:t1
    Alloc_dynamic(1+i,1) = i;
end
w=2;
for i=1:t1
    for j =1:p
        Alloc_dynamic(1+i,1+j) = res_dynamic(w);
        w=w+1;
    end
    w = 2 + (1+2*p)*i;
end


%linear allocation t*n
     res_linear = y.x';
     Alloc_linear = zeros(t+1,3+2*n);
for i = 1:n
    Alloc_linear(1,1+i) = i;
    Alloc_linear(1,1+i+n) = i;
    
end
Alloc_linear(1,2+2*n) = 1;
Alloc_linear(1,3+2*n) = 1;
for i = 1:t
    Alloc_linear(1+i,1) = i;
end
w = 1;h= 1+n*t;
for i=1:t
    for j =1:n
        Alloc_linear(1+i,1+j) = res_linear(w);
        Alloc_linear(1+i,1+j+n) = res_linear(w+t*(n+2));
        
        w=w+t;
    end
    Alloc_linear(1+i,2+2*n) = res_linear(h);
    Alloc_linear(1+i,3+2*n) = res_linear(h+t*(n+2));
    h=h+1;
    w=1+i;
end
lin = Alloc_linear;
dyn = Alloc_dynamic;
%{
disp(Alloc_dynamic);
disp(Alloc_linear); 
%}
end