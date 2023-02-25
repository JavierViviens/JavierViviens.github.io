%=========================================================================
%=========================================================================
%                         Econometrics III
%                     EUI PhD programme. 22/23.
%                           Problem Set 3
%                        Javier Viviens Martín
%=========================================================================
%=========================================================================

% Clean and set the seed.
clear; clc;
close all
seed=12345;
rng(seed)

% Quarterly data from 1960:Q1 to 2019:Q4.
% Read data:
datuak = readtable('data.xls');    
inflation = table2array(datuak(1:end,2));      
int = table2array(datuak(1:end,3));
unemployment = table2array(datuak(1:end,4));

%Create variable for time:
time =1960:0.25:2019.75;

%Inflation is clearly not stationary, so I create the growth rate of it,
%which it is.
pi = diff(inflation);

%Take the log of unemployment:
unem = log(unemployment);

%Construct the data for the VAR:
vardata = [pi unem(2:end) int(2:end)] ;

%Estimate VAR(4):
p=4;
constant = true ;
[comp, u, beta, X, Y,c] = var(p, vardata, constant) ;
if abs(eig(comp)) >= 1
    disp('VAR not stable')
else
    disp('stable VAR')
end

Granger_test = granger(X(1:end,2:end),Y)

h = 50; %horizon.
n = size(vardata,2); %number of variables.
irf = cholesky(comp,u,n,h); %Impulse Response Functions.

iter = 1000; %Number of iterations of the bootstrap.
cholboot = cholesky_boot(p, vardata, constant,iter,h) ;

%Plot the IRF at the 68% confidence level:
confidence_level=68; %Confidence level desired.
lb = 50 - confidence_level/2;
ub = 50 + confidence_level/2;
lower_bound = nan(n,n,h);
upper_bound = nan(n,n,h);

%Compute the mean and the confidence interval for each of the variable and
%each of the shocks:
lower_bound(1,1,:) = prctile(cholboot(1,1,:,:),lb,4);
upper_bound(1,1,:) = prctile(cholboot(1,1,:,:),ub,4);

lower_bound(1,2,:) = prctile(cholboot(1,2,:,:),lb,4);
upper_bound(1,2,:) = prctile(cholboot(1,2,:,:),ub,4);

lower_bound(1,3,:) = prctile(cholboot(1,3,:,:),lb,4);
upper_bound(1,3,:) = prctile(cholboot(1,3,:,:),ub,4);

lower_bound(2,1,:) = prctile(cholboot(2,1,:,:),lb,4);
upper_bound(2,1,:) = prctile(cholboot(2,1,:,:),ub,4);

lower_bound(2,2,:) = prctile(cholboot(2,2,:,:),lb,4);
upper_bound(2,2,:) = prctile(cholboot(2,2,:,:),ub,4);

lower_bound(2,3,:) = prctile(cholboot(2,3,:,:),lb,4);
upper_bound(2,3,:) = prctile(cholboot(2,3,:,:),ub,4);

lower_bound(3,1,:) = prctile(cholboot(3,1,:,:),lb,4);
upper_bound(3,1,:) = prctile(cholboot(3,1,:,:),ub,4);

lower_bound(3,2,:) = prctile(cholboot(3,2,:,:),lb,4);
upper_bound(3,2,:) = prctile(cholboot(3,2,:,:),ub,4);

lower_bound(3,3,:) = prctile(cholboot(3,3,:,:),lb,4);
upper_bound(3,3,:) = prctile(cholboot(3,3,:,:),ub,4);
%Plot IRFs:
figure(1)
fill([1:h fliplr(1:h)] ,[squeeze(lower_bound(1,1,:))' fliplr(squeeze(upper_bound(1,1,:))')],...
        [0.4 1 0.8],'EdgeColor','None'); 
hold on;
plot(1:h,squeeze(irf(1,1,:)),'LineWidth',1.5,'Color',[0 0 0])
hold on
plot(1:h,squeeze(lower_bound(1,1,:)), ':','LineWidth',1.5,'Color',[0.3010 0.6350 0.9330])
hold on
plot(1:h,squeeze(upper_bound(1,1,:)),':','LineWidth',1.5,'Color',[0.3010 0.6350 0.9330])
xlabel(['$ horizon $'],'interpreter','latex')
xlim([1,h])
title('Impulse response function: Inflation on Inflation','interpreter','latex')

figure(2)
fill([1:h fliplr(1:h)] ,[squeeze(lower_bound(1,2,:))' fliplr(squeeze(upper_bound(1,2,:))')],...
        [0.4 1 0.8],'EdgeColor','None'); 
hold on;
plot(1:h,squeeze(irf(1,2,:)),'LineWidth',1.5,'Color',[0 0 0])
hold on
plot(1:h,squeeze(lower_bound(1,2,:)),':','LineWidth',1.5,'Color',[0.3010 0.6350 0.9330])
hold on
plot(1:h,squeeze(upper_bound(1,2,:)),':','LineWidth',1.5,'Color',[0.3010 0.6350 0.9330])
xlabel(['$ horizon $'],'interpreter','latex')
xlim([1,h])
title('Impulse response function: Unemployment on Inflation','interpreter','latex')

figure(3)
fill([1:h fliplr(1:h)] ,[squeeze(lower_bound(1,3,:))' fliplr(squeeze(upper_bound(1,3,:))')],...
        [0.4 1 0.8],'EdgeColor','None'); 
hold on;
plot(1:h,squeeze(irf(1,3,:)),'LineWidth',1.5,'Color',[0 0 0])
hold on
plot(1:h,squeeze(lower_bound(1,3,:)),':','LineWidth',1.5,'Color',[0.3010 0.6350 0.9330])
hold on
plot(1:h,squeeze(upper_bound(1,3,:)),':','LineWidth',1.5,'Color',[0.3010 0.6350 0.9330])
xlabel(['$ horizon $'],'interpreter','latex')
xlim([1,h])
title('Impulse response function: Interest rate on Inflation','interpreter','latex')

figure(4)
fill([1:h fliplr(1:h)] ,[squeeze(lower_bound(2,1,:))' fliplr(squeeze(upper_bound(2,1,:))')],...
        [0.4 1 0.8],'EdgeColor','None'); 
hold on;
plot(1:h,squeeze(irf(2,1,:)),'LineWidth',1.5,'Color',[0 0 0])
hold on
plot(1:h,squeeze(lower_bound(2,1,:)), ':','LineWidth',1.5,'Color',[0.3010 0.6350 0.9330])
hold on
plot(1:h,squeeze(upper_bound(2,1,:)),':','LineWidth',1.5,'Color',[0.3010 0.6350 0.9330])
xlabel(['$ horizon $'],'interpreter','latex')
xlim([1,h])
title('Impulse response function: Inflation on Unemployment','interpreter','latex')

figure(5)
fill([1:h fliplr(1:h)] ,[squeeze(lower_bound(2,2,:))' fliplr(squeeze(upper_bound(2,2,:))')],...
        [0.4 1 0.8],'EdgeColor','None'); 
hold on;
plot(1:h,squeeze(irf(2,2,:)),'LineWidth',1.5,'Color',[0 0 0])
hold on
plot(1:h,squeeze(lower_bound(2,2,:)),':','LineWidth',1.5,'Color',[0.3010 0.6350 0.9330])
hold on
plot(1:h,squeeze(upper_bound(2,2,:)),':','LineWidth',1.5,'Color',[0.3010 0.6350 0.9330])
xlabel(['$ horizon $'],'interpreter','latex')
xlim([1,h])
title('Impulse response function: Unemployment on Unemployment','interpreter','latex')

figure(6)
fill([1:h fliplr(1:h)] ,[squeeze(lower_bound(2,3,:))' fliplr(squeeze(upper_bound(2,3,:))')],...
        [0.4 1 0.8],'EdgeColor','None'); 
hold on;
plot(1:h,squeeze(irf(2,3,:)),'LineWidth',1.5,'Color',[0 0 0])
hold on
plot(1:h,squeeze(lower_bound(2,3,:)),':','LineWidth',1.5,'Color',[0.3010 0.6350 0.9330])
hold on
plot(1:h,squeeze(upper_bound(2,3,:)),':','LineWidth',1.5,'Color',[0.3010 0.6350 0.9330])
xlabel(['$ horizon $'],'interpreter','latex')
xlim([1,h])
title('Impulse response function: Interest rate on Unemployment','interpreter','latex')

figure(7)
fill([1:h fliplr(1:h)] ,[squeeze(lower_bound(3,1,:))' fliplr(squeeze(upper_bound(3,1,:))')],...
        [0.4 1 0.8],'EdgeColor','None'); 
hold on;
plot(1:h,squeeze(irf(3,1,:)),'LineWidth',1.5,'Color',[0 0 0])
hold on
plot(1:h,squeeze(lower_bound(3,1,:)), ':','LineWidth',1.5,'Color',[0.3010 0.6350 0.9330])
hold on
plot(1:h,squeeze(upper_bound(3,1,:)),':','LineWidth',1.5,'Color',[0.3010 0.6350 0.9330])
xlabel(['$ horizon $'],'interpreter','latex')
xlim([1,h])
title('Impulse response function: Inflation on Interest rate','interpreter','latex')

figure(8)
fill([1:h fliplr(1:h)] ,[squeeze(lower_bound(3,2,:))' fliplr(squeeze(upper_bound(3,2,:))')],...
        [0.4 1 0.8],'EdgeColor','None'); 
hold on;
plot(1:h,squeeze(irf(3,2,:)),'LineWidth',1.5,'Color',[0 0 0])
hold on
plot(1:h,squeeze(lower_bound(3,2,:)),':','LineWidth',1.5,'Color',[0.3010 0.6350 0.9330])
hold on
plot(1:h,squeeze(upper_bound(3,2,:)),':','LineWidth',1.5,'Color',[0.3010 0.6350 0.9330])
xlabel(['$ horizon $'],'interpreter','latex')
xlim([1,h])
title('Impulse response function: Unemployment on Interest rate','interpreter','latex')

figure(9)
fill([1:h fliplr(1:h)] ,[squeeze(lower_bound(3,3,:))' fliplr(squeeze(upper_bound(3,3,:))')],...
        [0.4 1 0.8],'EdgeColor','None'); 
hold on;
plot(1:h,squeeze(irf(3,3,:)),'LineWidth',1.5,'Color',[0 0 0])
hold on
plot(1:h,squeeze(lower_bound(3,3,:)),':','LineWidth',1.5,'Color',[0.3010 0.6350 0.9330])
hold on
plot(1:h,squeeze(upper_bound(3,3,:)),':','LineWidth',1.5,'Color',[0.3010 0.6350 0.9330])
xlabel(['$ horizon $'],'interpreter','latex')
xlim([1,h])
title('Impulse response function: Interest rate on Interest rate','interpreter','latex')

%=========================================
%% ============ Functions ================
%=========================================
% ------------------------------- 
% --- Granger-Causality Tests ---
% -------------------------------
function [Granger_test] = granger(X,Y)
    Granger_test = zeros(size(Y,2),size(Y,2));
    for i = 1:size(Y,2)
        ols = fitlm(X,Y(:,i)); %Run OLS.
        Granger_test(i,1) = coefTest(ols,[0,1,0,0,1,0,0,1,0,0,1,0,0]); %Test inflation has no effect.
        Granger_test(i,2) = coefTest(ols,[0,0,1,0,0,1,0,0,1,0,0,1,0]); %Test unemployment has no effect.
        Granger_test(i,3) = coefTest(ols,[0,0,0,1,0,0,1,0,0,1,0,0,1]); %Test interest rate has no effect.
    end
end

% ------------------------------- 
% ------------- VAR -------------
% -------------------------------
function [B, u, beta, X, Y,c] = var(p, datuak, constant)
% B: Companion form matrix
% p: VAR order
% constant: Boolean
    Y = datuak(1+p:end,:);
    n = size(datuak,2);
    T = size(datuak,1);
    X = ones(T-p,n*p);
    %SUR representation
    for i = 1:p
        X(:,(i-1)*n+1:i*n) = datuak(p-i+1:end-i,:);
    end
    if constant
        X = [ones(T-p,1) X]; 
        beta = inv(X'*X)*X'*Y;
        B=[beta(2:end,:)';eye(n*(p-1)) zeros(n*(p-1),n)];
        c = [beta(1,:)';zeros(n*(p-1),1)];
    else
        beta = inv(X'*X)*X'*Y;
        B=[beta';eye(n*(p-1)) zeros(n*(p-1),n)];
    end
    u = Y - X*beta;
end
% ------------------------------- 
% ---------- CHOLESKY -----------
% -------------------------------
 function Chol_IRF = cholesky(companion,u,n,h)
    np = size(companion,1);
    CL = nan(np,np,h);
    Chol_IRF = nan(n,n,h);
    Omega = cov(u);
    for i=1:h %For each period in the horizon, we get the IRF
        CL(:,:,i) = companion^(i-1);
        Chol_IRF(:,:,i) = CL(1:n,1:n,i)*chol(Omega,'lower');
    end
 end
% ------------------------------- 
% ------ CHOLESKY BOOTSRAP ------
% -------------------------------
 function cholirfboot = cholesky_boot(p, datuak, constant,iter,h)
    n = size(datuak,2);
    T = size(datuak,1);
    [comp, u, beta, X, Y] = var(p,datuak,constant);
    chol_irf = cholesky(comp,u,n,h);
    cholirfboot = nan(n,n,h,iter);
    for i=1:iter
        synY = nan(T-p,n); %Save memory for Synthetic Y
        if constant
            synX = ones(T-p,n*p+1); %Save memory for Synthetic X
        else
            synX = ones(T-p,n*p);
        end
        for j = 1:size(synY,1)
            if j==1
                synX(1,:) = X(1,:); 
            else
                if constant
                    synX(j,:) = [1 synY(j-1,:) synX(j-1,2:n*(p-1)+1)];
                else
                    synX(j,:) = [synY(j-1,:) synX(j-1, 1:n*(p-1))];
                end
            end
            %We construct the synthetic Y by getting the fitted value of
            %the synthetic X and adding one of the residuals of the VAR
            %estimation randomly picked.
            synY(j,:) = synX(j,:)*beta + u(randi([1,T-p]),:); 
        end
        
         [comp_loop, u_loop, ~, ~, ~] = var(p,synY,constant);
         chol_irf = cholesky(comp_loop,u_loop,n,h);
         cholirfboot(:,:,:,i) = chol_irf;
    end
 end