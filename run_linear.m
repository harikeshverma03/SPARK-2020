function [Alloc_linear,sum_reg,sum_hol,objval_mu,objval_sd,Dem_mean,sol_mean] = run_linear(Day,Beta,Hydro_avg,Demand_mu,Demand_sd,Cost_imp_mu,Cost_imp_sd,Sol_mu,Sol_sd,p)
    n = 7;
    t = 24;sum1 = 0;objval = 0;Dem_mean=0;sol_mean = 00;
    Hydro_cap = 447;
    Thermal_cap = [180 250 250 250 500 500 150];
    Cost_1 = [2.81 2.88 2.88 4.06 4.15 10.76 14.18];
    Ramp_Rate = [0.3050    0.0610    0.0610    0.0610    0.0610    0.3050    0.6100];
    Sol_E = zeros(24,2,p);
    Demand = zeros(24,2,p);
    Cost_imp = zeros(1,p);
     for i = 1:p
        Sol_E(:,:,i) = Sol_mu + Sol_sd*randn(1,1);
        Demand(:,:,i) = Demand_mu + Demand_sd*randn(1,1);
        Cost_imp(i) = Cost_imp_mu + Cost_imp_sd*randn(1,1);
     end
    %{
    Day = [26,5];
    Demand = csvread('Demand_Jan.csv');
    Cost_imp = 15;
    Hydro_avg = 16.0895;
    Sol_E = csvread('Solar_Jan.csv');
    %}
     count = 0;
     for i = 1:p
         if Demand(:,:,i) - Sol_E(:,:,i) > 0.2
            
            y = func.linear_stage(Day, Cost_1,Cost_imp(i),Demand(:,:,i),Hydro_cap, Thermal_cap,Hydro_avg,Ramp_Rate,Sol_E(:,:,i),n,t,Beta);
            if isfield(y,'x')
                count = count +1; 
                z(i,:) = y.x';
                Dem_mean = Demand(:,:,i) + Dem_mean;
                sol_mean = sol_mean + Sol_E(:,:,i)*Beta;
                objval(1,count) = y.objval;
                sum1 = sum1 + z(i,:);
            else
            end
        end
     end
     disp(count);
    res_linear = sum1/count;
    Dem_mean = Dem_mean/count;
    sol_mean = sol_mean/count;
    objval_mu = mean(objval);
    objval_sd = std(objval);
    fprintf('The value of Objective function is: %e\n', objval_mu);
    [Alloc_linear,sum_reg,sum_hol] = func.Alloc_linear(res_linear,t,n);
end