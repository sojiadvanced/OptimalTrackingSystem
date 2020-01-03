clear all
clc
%Adesoji Bello
%Design of a Radar Measurement System for Airport Surveillance

%%
%Measurements of Airborne Target
k = 10000;
theta_target = zeros(k+1,1);                 %Values furnished to the Airborne Target
w = unifrnd(-0.01,0.01,k+1,1);

for i =1:1:k
 	
    theta_target(i+1) =  theta_target(i) + w(i) ;
end


%%
%Radar Measurements

theta_meas = zeros(k+1,1);
n = normrnd(0,sqrt(0.02),k+1,1);
for j =1:1:k+1
    
    theta_meas(j) = theta_target(j) + n(j);
end

%%
%Obtain the estimated values

A= 0:0.01:1;
m = length(A);
theta_estimate = zeros(k+1,m);
for i =1:1:length(A)
    
    for j =1:1:k
 	
        theta_estimate(j+1,i) = A(i)* theta_estimate(j,i) + (1 - A(i))* theta_meas(j) ;
    end  
end


%%
%Computation of estimation error & variance
est_error = zeros(k+1,m);
for i =1:1:m
    
    for j =1:1:k+1
 	
        est_error(j,i) =theta_estimate(j,i) - theta_target(j)  ;
    end  
end



var_error = var(est_error,1,1);
measured_errvar = min(var_error);
Ameasured_val = A(find(var_error ==measured_errvar));

fprintf("The simulated estimation error variance is: %0.5f\n",measured_errvar);
fprintf("The value of A with least simulated estimate error is: %0.2f\n\n",Ameasured_val);


%%
%Plot of True Azimuth values & Estimated Azimuth values

minvar = find(var_error == min(var_error));
optimal_est =  theta_estimate(:,minvar);
figure(1)
plot([0:k],theta_target)
hold on
plot([0:k],optimal_est)
title("Plot of Azimuth True & Estimated Values")
xlabel('Range of samples')
ylabel('Azimuth Values')
legend('True Value', 'Estimated Value')
grid on
hold off


%%
%Impulse response of the system with w[k]


theta_Imp = zeros(k+1,1);                 %Values furnished to the Airborne Target
w_imp = zeros(k+1,1);
w_imp(1) = 1;

for i =1:1:k
 	
    theta_Imp(i+1) =  theta_Imp(i) + w_imp(i) ;
end
theta_estwk = zeros(k+1,m);
imp_estwk = zeros(k+1,m);
for i =1:1:length(A)
   
    for j =1:1:k
 	
        theta_estwk(j+1,i) = A(i)* theta_estwk(j,i) + (1 - A(i))* (theta_Imp(j)) ;
    end  
    for j=1:1:k
         imp_estwk(j,i) = theta_Imp(j,1) - theta_estwk(j,i);
    end
end

imp2_errwk = zeros(1,m);
%Mean Square Value of the impulse response
for i =1:1:m
    
   for j=1:1:k+1
       imp2_errwk(1,i) = imp2_errwk(1,i) + imp_estwk(j,i)^2;
   end
end
var_wk = (0.02)^2 /12;
est_errwknew = imp2_errwk * var_wk;       %Output response based on w[k] only

% optimal_esterr = min(est_errwknew);
% Aoptimal_val = A(find(est_errwknew ==optimal_esterr));
% 
% fprintf("The estimation error variance for w[k] only is: %0.5f\n",optimal_esterr);
% fprintf("The value of A with least measured estimate error is: %0.2f\n\n",Aoptimal_val);
figure(2)
plot(A,est_errwknew)
title("Plot of Estimation Error Variance with w[k] only")
xlabel("A Values")
ylabel("Estimated Error Variance")
grid on
 
figure(3)
plot([0:k],imp_estwk(:,97))
title("Plot of Impulse Response with w[k] only when A =0.96")
grid on
%%
%Impulse response of the system with n[k]
theta_targetnk = zeros(k+1,1);
theta_measimp = zeros(k+1,1);
n_imp = zeros(k+1,1);
n_imp(1) = 1;
imp_estnk = zeros(k+1,m);
for j =1:1:k+1
    
    theta_measimp(j) = theta_targetnk(j) + n_imp(j);
end
theta_estimp = zeros(k+1,m);
for i =1:1:length(A)
    
    for j =1:1:k
 	
        theta_estimp(j+1,i) = A(i)* theta_estimp(j,i) + (1 - A(i))* theta_measimp(j) ;
    end 
    
    for j=1:1:k
         imp_estnk(j,i) =theta_targetnk(j,1) - theta_estimp(j,i);   %Impulse response with n[k] only
    end
end


var_nk = 0.02;
est_errnk = zeros(1,m);
for i =1:1:m
    
   for j=1:1:k+1
       est_errnk(1,i) = est_errnk(1,i) + imp_estnk(j,i)^2;
   end
end
est_errnknew = est_errnk * var_nk;
% optimal_esterr = min(est_errnknew);
% Aoptimal_val = A(find(est_errnk ==optimal_esterr));
% 
% fprintf("The estimation error variance for n[k] only is: %0.5f\n",optimal_esterr);
% fprintf("The value of A with least measured estimate error is: %0.2f\n\n",Aoptimal_val);
figure(4)
plot(A,est_errnknew)
title("Plot of Estimation Error Variance with n[k] only")
xlabel("A Values")
ylabel("Estimated Error Variance")
grid on
% 
figure(5)
plot([0:k],imp_estnk(:,97))
title("Plot of Impulse Response with n[k] only when A = 0.96")
grid on

%%
%Impulse response with both w[k] and n[k]

est_errwknk = est_errnknew + est_errwknew;


figure(6)
plot(A,est_errwknew)
hold on
plot(A,est_errnknew)
plot(A,est_errwknk)
title("Plot of Estimation Error Variance with w[k], n[k] and both")
xlabel("A Values")
ylabel("Estimated Error Variance")
grid on
legend('w[k] only', 'n[k] only','w[k and n[k]')
hold off

optimal_estwknk = min(est_errwknk);
Aoptimal_valwknk = A(find(est_errwknk ==optimal_estwknk));

optimal_esterr = est_errnknew(find(est_errwknk == optimal_estwknk));

optimal_esterrwk = est_errwknew(find(est_errwknk == optimal_estwknk));
fprintf("The estimation error variance for n[k] only is: %0.6f\n",optimal_esterr);
% fprintf("The value of A with least measured estimate error is: %0.2f\n\n",Aoptimal_val);

fprintf("The estimation error variance for w[k] only is: %0.6f\n",optimal_esterrwk);

fprintf("The optimal estimation error variance for both noise is: %0.6f\n",optimal_estwknk);
fprintf("The value of A with least measured estimate error is: %0.2f\n",Aoptimal_valwknk);


