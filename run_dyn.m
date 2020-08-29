function [dyn,sum,objval] = run_dyn(n,Cost_imp,Beta,Rain_mu,Rain_sd,V1_mu,V1_sd,E_mu,E_sd)
 t1 = 12;
 p = 3;
 z = zeros(n,84);
 objval = 0;
 MDDL = [619.4 646 590.1];
 FRL = [635.2 668 607.1];
 Hp = [51.84 54 216];
 Avg = 1336.232;
 E = zeros(n,t1);
 V1 = zeros(n,p);
 Rain = zeros(t1,p,n);
 
 for i = 1:n
    E(i,:) = E_mu + E_sd*randn(1,1);
    V1(i,:) = V1_mu + V1_sd*randn(1,1);
    Rain(:,:,i) = Rain_mu + Rain_sd*randn(1,1);
 end
 Dem = csvread('Demand.csv')';
 sum1 = 0;count = 0;
 for i = 1:n
     if V1(i,:) - MDDL > 0 
        Dem_cost = Dem - Beta*E_mu;
        Cost = func.Cost(Dem_cost,Cost_imp);
        x = func.dynamic_stage(Cost,Dem,E(i,:),Beta,V1(i,:),Rain(:,:,i),MDDL,FRL,Hp,t1,p,Avg,V1_mu);
        if isfield(x,'x')
        z(i,:) = x.x';
        objval(1,i) = x.objval;
        sum1 = sum1 + z(i,:);
        count = count +1;
        else
        end
     end
 end
res_dynamic = sum1/count;
fprintf('The value of Objective function is: %e\n', mean(objval));      
[dyn,sum] = func.Alloc_dynamic(res_dynamic,t1,p);
x = 1:12;
plot(x,sum);
end