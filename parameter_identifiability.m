% Supplementary Analysis: Parameter Identifiability
% identifying the correct parameters for a two population model
close all; clear all; clc
% First make a simulated data set with reasonable model parameters

% Find variability as a function of dose
pat = xlsread('../data/subpop_sorting_Cayman.xls');
%% Find variability as a function of dose
dose = pat(:,2);
viability = pat(:,6);
figure(1)
plot(dose, viability, 'o')
xlabel('dose(\muM)')
ylabel('viability')
title('Raw data all mixtures')

[d, I] =sort(dose);
v = viability(I);

uniq_dose = unique(dose);
for i = 1:length(uniq_dose)
    iv = find(ismember(dose, uniq_dose(i)));
    viab = viability(iv);
    variability(i) = std(viab);
    
end

figure(2)
plot(uniq_dose, variability, 'go-', 'LineWidth', 3)
xlabel ('dose (\muM)')
ylabel('Standard deviation of viability')
%title('Standard deviation in viability as a function of dose')



%%

params = [ 185 0.03 35 0.06 0 0.25 0.50 0.75 1]; % make fake params
nsamp = length(params)-4;
Vmaxbymix =  [0.9148    0.9254    0.9122    0.9327    0.9400]; % from actual data
dose = [  0 18 54 90 130 180 240 300 375 450 525 600];   
D = 1:1:max(dose);
model_sim = @(p)two_pop_model(p, Vmaxbymix, dose);
model = model_sim(params);
simmodel = repmat(model,1,3);
simdose = repmat(dose,1,3);
simeta = repmat(variability, 1, 3);
%eta = 0.15;

noise =2*simeta.*(0.5-rand(nsamp, length(simmodel)));
simdata = simmodel + noise;


for i = 1:size(simdata,1)
    for j = 1:size(simdata,2)
        if simdata(i,j) < 0
            simdata(i,j) = 0;
        end
        if simdata(i,j) >1
            simdata(i,j) = 1;
        end
    end
end 
    

Color = {'b'; 'c';'g'; 'm'; 'r'};
figure(1);
hold off
for i = 1:nsamp
    plot(dose, model(i,:), Color{i},'LineWidth',3)
    hold on
end
for i = 1:nsamp
    plot(simdose, simdata(i,:), '*', 'color', Color{i})
    %text(400, model(i, 4)+0.1, ['f_{res}=, ',num2str(round(params(4+i),2)), ', f_{sens}=, ', num2str(round(1-params(4+i),2))], 'HorizontalAlignment','left','VerticalAlignment','top','color', Color{i})
    hold on
end

xlabel( 'Dose (\muM)')
ylabel('Viability')
legend('0% ADR', '25% ADR','50% ADR','75% ADR','100% ADR')
% legend( '0% ADR_{model}','0% ADR_{meas}', '25% ADR_{model}','25% ADR_{meas}' ,'50% ADR_{model}','50% ADR_{meas}',  '75% ADR_{model}', '75% ADR_{meas}','100% ADR_{model}', '100% ADR_{meas}')
% legend boxoff
%title('Two Population Model Mixtures: Simulated Data Set, \eta = f(dose)_{exp}', 'FontSize', 10)

%% Now find parameters

ns = 12;
nreps = 3;
% something is wrong with restaching of data
doselong = repmat(dose, 1, nreps*nsamp);
simdatalong = [];
for i = 1:nsamp
    simdatalong = vertcat(simdatalong, simdata(i,:)');
end

ind = find(doselong == 0);
Vmax = simdatalong(ind);
Vmaxbymixsim = [];
Vmaxallsim = [];
for j = 1:nsamp
    Vmaxbymixsim(j) = mean(Vmax(3*j-2:3*j));
    Vmaxmat = repmat(Vmaxbymixsim(j), ns.*nreps,1);
    Vmaxallsim = vertcat(Vmaxallsim, Vmaxmat);
end 
% make vertical
doselong = doselong';
for i = 1:length(simdatalong)
    if simdatalong(i) >1
        simdatalong(i) = 1;
    end
    if simdatalong(i)< 0
        simdatalong(i) = 0;
    end
end



%% Start with same initial params as in model
params0 = horzcat( [200 0.01 25 0.01], 0.5.*ones(1, nsamp));%  LD50res, sloperes, LD50sens, slopesens, fres
paramslb = zeros( 1, 4+nsamp);
paramsub = horzcat( [ Inf 1 Inf 1], ones(1, nsamp));
options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
    
    ihi = 145:1:180;
    ilo = 1:1:36;


    psingle0 = [ 70 0.1];

[Phisim, resnormhi, residualshi] = lsqnonlin(@fitsinglepop,...
    psingle0,...
    paramslb(1:2),...
    paramsub(1:2),...
    options,...
    doselong(ihi),...
    simdatalong(ihi));
[Plosim, resnorm, residuals] = lsqnonlin(@fitsinglepop,...
    psingle0,...
    paramslb(1:2),...
    paramsub(1:2),...
    options,...
    doselong(ilo),...
    simdatalong(ilo));

[Psim, resnorm, residuals] = lsqnonlin(@fitmixedpops,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    doselong,...
    simdatalong,...
    nsamp,...
    Vmaxallsim);
%% How sensitive is the model to noise at current parameter values??

% What we do: simulate 50 runs for with noise as a function of dose, outputting parameter
% estimates for each run at each noise level.

% From these distributions of models, can determine a few things:
% 1. parameter distribution, can plot these and see when parameters are the
% most senstive to noice
% 2. Need to assess how sensitive model is to noise-- how does noise
% correlate to spread of parameter distribution?

% Make model vector
modellong = [];
for i = 1:nsamp
    modellong = vertcat(modellong, simmodel(i,:)');
end
simetalong = repmat(simeta, 1, 5);

% Add noise and fit in a loop, output parameters
%eta =  0.3;
simmodelong = [];
Psim = [];
Plosim = [];
Phisim = [];

for k = 1:100
    noise = simetalong.*(0.5-rand(1, length(modellong)));
    simmodellong = modellong + noise';
    
    for i = 1:length(simmodellong)
        if simmodellong(i) >1
            simmodellong(i) = 1;
        end
        if simmodellong(i)< 0
            simmodellong(i) = 0;
        end
    end
    
    ind = find(doselong == 0);
    Vmax = simdatalong(ind);
    Vmaxbymixsim = [];
    Vmaxallsim = [];
    for j = 1:nsamp
        Vmaxbymixsim(j) = mean(Vmax(3*j-2:3*j));
        Vmaxmat = repmat(Vmaxbymixsim(j), ns.*nreps,1);
        Vmaxallsim = vertcat(Vmaxallsim, Vmaxmat);
    end
    
    [Phisim(:,k)] = lsqnonlin(@fitsinglepop,...
    psingle0,...
    paramslb(1:2),...
    paramsub(1:2),...
    options,...
    doselong(ihi),...
    simmodellong(ihi));
[Plosim(:,k)] = lsqnonlin(@fitsinglepop,...
    psingle0,...
    paramslb(1:2),...
    paramsub(1:2),...
    options,...
    doselong(ilo),...
    simmodellong(ilo));

[Psim(:,k)] = lsqnonlin(@fitmixedpops,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    doselong,...
    simmodellong,...
    nsamp,...
    Vmaxallsim);

end

% Compute stats and plot distribution of parameters

mean_param_val = mean(Psim,2);
std_param_val = std(Psim,1,2);
lower_param_val = mean_param_val - 1.96*std_param_val;
upper_param_val = mean_param_val + 1.96*std_param_val;

for j = 1:length(mean_param_val)
    sum_err = 0;
    for i = 1:length(Psim)
        err = (params(j) - Psim(j, i))^2;
        sum_err = sum_err + err;
    end
    MSE(j) = sum_err/(length(Psim)-length(params)); % mean squared error
    CV(j) = std_param_val(j)/(mean_param_val(j));
    
end
        
comb_tbl = vertcat(mean_param_val', params, lower_param_val', upper_param_val', std_param_val', MSE, CV);
%
param_distrib = dataset({comb_tbl, 'LD50_res', 'slope_res', 'LD50_sens', 'slope_sens', 'fres1', 'fres2', 'fres3','fres4', 'fres5'}, 'obsnames',{'mean','actual', 'lower', 'upper', 'std_dev', 'MSE', 'CV'});
        save('../out/param_distrib.mat', 'param_distrib')
errorbars = upper_param_val-lower_param_val;
LD50s=[ mean_param_val(1), mean_param_val(3)];
errorbarLD50s = [errorbars(1), errorbars(3)];
slopes=[ mean_param_val(2), mean_param_val(4)];
errorbarslopes= [errorbars(2), errorbars(4)];
fracs = mean_param_val(5:9);
errorbarfracs = errorbars(5:9);

figure(4)
hold off
subplot(2, 2, 1)
hold off
set(gca,'LineWidth',1.2,'FontSize',10)
errorbar(1:1:2, LD50s, errorbarLD50s, 'o')
hold on
plot(1:1:2, vertcat(params(1), params(3)), 'g*', 'LineWidth',2)
xlim ([0 3])
ylim([ 0 225])
ylabel('LD50 (\muM)')

title('LD50 distribution')

subplot(2,2,2)
hold off
set(gca,'LineWidth',1.2,'FontSize',10)
errorbar(1:1:2, slopes, errorbarslopes, 'o')
hold on
plot(1:1:2, vertcat(params(2), params(4)), 'g*','LineWidth',2)
xlim([0 3])
ylabel ('slope')
title('slope distribution')

subplot(2,2,3:4)
hold off
set(gca,'LineWidth',1.2,'FontSize',10)
errorbar(1:1:5, fracs, errorbarfracs, 'o')
hold on
plot(1:1:5, params(5:9), 'g*','LineWidth',2)
title('fraction estimate distribution')
xlim ([0 6])
ylim([-0.1 1.1])
ylabel('fraction resistant')
xlabel('mixtures')



%% How well is the model able to identify proportions very close to one another?
% Add noise and fit in a loop, output parameters
%eta =  0.3;
paramshigh = [185.0000, 0.0300, 35.0000, 0.0600 0.7, 0.8, 0.9,0.95,1];
model= model_sim(paramshigh)
simmodel = repmat(model,1,3)
%%
% Make model vector
modellong = [];
for i = 1:nsamp
    modellong = vertcat(modellong, simmodel(i,:)');
end
simetalong = repmat(simeta, 1, 5);
%%
simmodelong = [];
Psim = [];
Plosim = [];
Phisim = [];

for k = 1:100
    noise = 2.*simetalong.*(0.5-rand(1, length(modellong)));
    simmodellong = modellong + noise';
    
    for i = 1:length(simmodellong)
        if simmodellong(i) >1
            simmodellong(i) = 1;
        end
        if simmodellong(i)< 0
            simmodellong(i) = 0;
        end
    end
    
    ind = find(doselong == 0);
    Vmax = simdatalong(ind);
    Vmaxbymixsim = [];
    Vmaxallsim = [];
    for j = 1:nsamp
        Vmaxbymixsim(j) = mean(Vmax(3*j-2:3*j));
        Vmaxmat = repmat(Vmaxbymixsim(j), ns.*nreps,1);
        Vmaxallsim = vertcat(Vmaxallsim, Vmaxmat);
    end
    
    [Phisim(:,k)] = lsqnonlin(@fitsinglepop,...
    psingle0,...
    paramslb(1:2),...
    paramsub(1:2),...
    options,...
    doselong(ihi),...
    simmodellong(ihi));
[Plosim(:,k)] = lsqnonlin(@fitsinglepop,...
    psingle0,...
    paramslb(1:2),...
    paramsub(1:2),...
    options,...
    doselong(ilo),...
    simmodellong(ilo));

[Psim(:,k)] = lsqnonlin(@fitmixedpops,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    doselong,...
    simmodellong,...
    nsamp,...
    Vmaxallsim);
    
%     subplot(1,1, j)
%     hold off
%     plot(doselong, simmodellong, '*')
%     hold on
%     plot(doselong, modellong, 'ko')
%     
end

% Compute stats and plot distribution of parameters

mean_param_val = mean(Psim,2);
std_param_val = std(Psim,1,2);
lower_param_val = mean_param_val - 1.96*std_param_val;
upper_param_val = mean_param_val + 1.96*std_param_val;

for j = 1:length(mean_param_val)
    sum_err = 0;
    for i = 1:length(Psim)
        err = (params(j) - Psim(j, i))^2;
        sum_err = sum_err + err;
    end
    MSE(j) = sum_err/(length(Psim)-length(params)); % mean squared error
    CV(j) = std_param_val(j)/(mean_param_val(j));
    
end
        
comb_tbl = vertcat(mean_param_val', params, lower_param_val', upper_param_val', std_param_val', MSE, CV);
%
param_distrib = dataset({comb_tbl, 'LD50_res', 'slope_res', 'LD50_sens', 'slope_sens', 'fres1', 'fres2', 'fres3','fres4', 'fres5'}, 'obsnames',{'mean','actual', 'lower', 'upper', 'std_dev', 'MSE', 'CV'});
        save('../out/param_distrib.mat', 'param_distrib')
errorbars = upper_param_val-lower_param_val;
LD50s=[ mean_param_val(1), mean_param_val(3)];
errorbarLD50s = [errorbars(1), errorbars(3)];
slopes=[ mean_param_val(2), mean_param_val(4)];
errorbarslopes= [errorbars(2), errorbars(4)];
fracs = mean_param_val(5:9);
errorbarfracs = errorbars(5:9);

figure(4)
hold off
subplot(2, 2, 1)
hold off
set(gca,'LineWidth',1.2,'FontSize',10)
errorbar(1:1:2, LD50s, errorbarLD50s, 'o')
hold on
plot(1:1:2, vertcat(params(1), params(3)), 'g*', 'LineWidth',2)
xlim ([0 3])
ylim([ 0 225])
ylabel('LD50 (\muM)')

title('LD50 distribution')

subplot(2,2,2)
hold off
set(gca,'LineWidth',1.2,'FontSize',10)
errorbar(1:1:2, slopes, errorbarslopes, 'o')
hold on
plot(1:1:2, vertcat(params(2), params(4)), 'g*','LineWidth',2)
xlim([0 3])
ylabel ('slope')
title('slope distribution')

subplot(2,2,3:4)
hold off
set(gca,'LineWidth',1.2,'FontSize',10)
errorbar(1:1:5, fracs, errorbarfracs, 'o')
hold on
plot(1:1:5, params(5:9), 'g*','LineWidth',2)
title('fraction estimate distribution')
xlim ([0 6])
ylim([-0.1 1.1])
ylabel('fraction resistant')
xlabel('mixtures')

%%
X = Psim;
Xnew = reshape(X, size(Psim,1), size(Psim, 2));
X = Xnew';


figure(2)
hold off
[S,AX,BigAx,H,HAx] = plotmatrix(X(:,1:4));
title ('Distribution of LD50 and slope parameters')
ylabel(AX(1,1),'LD50_{res}')
ylabel(AX(2,1),'slope_{res}')
ylabel(AX(3,1), 'LD50_{sens}')
ylabel(AX(4,1),'slope_{sens}')
xlabel(AX(4,1),'LD50_{res}')
xlabel(AX(4,2),'slope_{res}')
xlabel(AX(4,3), 'LD50_{sens}')
xlabel(AX(4,4),'slope_{sens}')

figure(3)
hold off
[S,AX,BigAx,H,HAx] = plotmatrix(X(:,5:9));
title('Distribution of fraction estimate parameters')
ylabel(AX(1,1),'f_{res1}')
ylabel(AX(2,1),'f_{res2}')
ylabel(AX(3,1), 'f_{res3}')
ylabel(AX(4,1),'f_{res4}')
ylabel(AX(5,1),'f_{res5}')
xlabel(AX(5,1),'f_{res1}')
xlabel(AX(5,2),'f_{res2}')
xlabel(AX(5,3), 'f_{res2}')
xlabel(AX(5,4),'f_{res4}')
xlabel(AX(5,5),'f_{res5}')
%% Perturb parameters: Distance between LD50s and fres1, what is corresponding average error?


MSE_LD50 = [];
MSE_init_frac = [];
%for m = 1 % find MSE in  5 different sets of parameter conditions
params = [];
params(1,:) = [ 50 0.0172 35 0.0501 0 0.25 0.5 0.75 1]; % make fake params
% First just perturb LD50 res and plot corresponding error in key
% parameters
for h = 2:6
params(h,:) = [ 50+50*(h) 0.0172 35 0.0501 0 0.25 0.5 0.75 1]; 
end

dist_rmins = params(:,1)-params(:,3);
%%

for m = 1:6


LD50_dist(m,1) = params(m,1)- params(m,3);
init_frac(m,1) = params(5);
model_sim = @(p)two_pop_model(p, Vmaxbymix, dose);
model = model_sim(params(m,:));
simmodel = repmat(model,1,3);
simdose = repmat(dose,1,3);


% Start with initial parameter distribution and noise from before:
modellong = [];
for i = 1:nsamp
    modellong = vertcat(modellong, simmodel(i,:)');
end

% Add noise and fit in a loop, output parameters
eta =  0.3;
simmodelong = [];
Psim = [];
Plosim = [];
Phisim = [];


for k = 1:100
    noise = eta.*(0.5-rand(1, length(modellong)));
    simmodellong = modellong + noise';
    
    for i = 1:length(simmodellong)
        if simmodellong(i) >1
            simmodellong(i) = 1;
        end
        if simmodellong(i)< 0
            simmodellong(i) = 0;
        end
    end
    
    ind = find(doselong == 0);
    Vmax = simdatalong(ind);
    Vmaxbymixsim = [];
    Vmaxallsim = [];
    for j = 1:nsamp
        Vmaxbymixsim(j) = mean(Vmax(3*j-2:3*j));
        Vmaxmat = repmat(Vmaxbymixsim(j), ns.*nreps,1);
        Vmaxallsim = vertcat(Vmaxallsim, Vmaxmat);
    end
    
    [Phisim(:,k)] = lsqnonlin(@fitsinglepop,...
    psingle0,...
    paramslb(1:2),...
    paramsub(1:2),...
    options,...
    doselong(ihi),...
    simmodellong(ihi));
[Plosim(:,k)] = lsqnonlin(@fitsinglepop,...
    psingle0,...
    paramslb(1:2),...
    paramsub(1:2),...
    options,...
    doselong(ilo),...
    simmodellong(ilo));

[Psim(:,k)] = lsqnonlin(@fitmixedpops,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    doselong,...
    simmodellong,...
    nsamp,...
    Vmaxallsim);

if Psim(1,k) < Psim(3,k)
    Psimold = Psim(:,k);
    Psim(1,k)= Psimold(3);
    Psim(2,k) = Psimold(4);
    Psim(3,k)= Psimold(1);
    Psim(4,k) = Psimold(2);
end
  percent_error(:,k) =abs(100.*( Psim(:,k)-params(j,:)')./params(j,:)');


end

mean_param_val = mean(Psim,2);
std_param_val = std(Psim,1,2);
mean_percent_error = mean(percent_error,2);

lower_param_val = mean_param_val - 1.96*std_param_val;
upper_param_val = mean_param_val + 1.96*std_param_val;

    for j = 1:length(mean_param_val)
        sum_err = 0;
        for i = 1:length(Psim)
            err = (params(j) - Psim(j, i))^2;
            sum_err = sum_err + err;
        end
        MSE(j) = sum_err/(length(Psim)-length(params)); % mean squared error
        CV(j) = std_param_val(j)/(mean_param_val(j));

    
    end

% For this parameter distribution, find the MSE in LD50 estimates and fres
% estimates

MSE_LD50(m,1:2) = horzcat(MSE(1), MSE(3));
MSE_all_params(m,:) = MSE;
MSE_LD50_avg(m) = mean(MSE_LD50(m,:));
STD_LD50(m,1:2) = horzcat(std_param_val(1), std_param_val(3));
STD_LD50_avg(m,1:2) = mean(STD_LD50(m,:));
MSE_init_frac(m) = MSE(5);

mean_percent_error_LD50res(m) = mean_percent_error(1);
mean_percent_error_LD50sens(m) = mean_percent_error(3);

end
%% Visual Displays
Params = {'LD50_{res}', 'slope_{res}', 'LD50_{sens}', 'slope_{sens}', 'f_{res1}', 'f_{res2}', 'f_{res3}', 'f_{res4}', 'f_{res5}'}
MSE_norm = MSE_all_params./params;
figure(7)
hold off
set(gca,'LineWidth',1.5,'FontSize',14)
    for k =1:9 % 4%size(MSE_all_params,2) % for all 6 diferent param perturbations
        semilogy(params(:,1), MSE_norm(:,k),'o-', 'LineWidth', 2);
        text(params(4,1), MSE_norm(4,k), ['MSE in ', Params{k}], 'HorizontalAlignment','left','VerticalAlignment','top')
        hold on
    end
xlabel('LD50_{res} uM')
ylabel( 'Normalized Mean-Squared Error in Parameter')
title('Error in Parameters as a Function of LD50_{res}')
xlim([50 420])

% Make figure of error in LD50 res and sensitive as function of their
% distance
dist_rmins = params(:,1)-params(:,3);
MSE_res = MSE_norm(:,1);
MSE_sens = MSE_norm(:,3);
%%
figure (8)
set(gca,'LineWidth',1.5,'FontSize',14)
hold off
plot(dist_rmins, MSE_res, 'ro-', 'LineWidth', 2)
text(dist_rmins(2), MSE_res(1),'MSE LD50_{res}', 'VerticalAlignment', 'bottom')
hold on
plot(dist_rmins, MSE_sens, 'bo-', 'LineWidth',2)
text(dist_rmins(2), MSE_sens(1),'MSE LD50_{sens}', 'VerticalAlignment', 'bottom')
xlabel('LD50_{res} - LD50_{sens}')
ylabel('Normalized Mean-Squared Error in LD50 estimation')
title('Error vs. Distance between LD50s as LD50_{res} increases')
%% Repeat Same analysis but alter LD50_sens

MSE_LD50 = [];
MSE_init_frac = [];
%for m = 1 % find MSE in  5 different sets of parameter conditions
params2 = [];
params2(1,:) = [ 185 0.0172 100 0.05 0 0.25 0.5 0.75 1]; % make fake params
% First just perturb LD50 res and plot corresponding error in key
% parameters
for h = 2:6
params2(h,:) = [ 150 0.0172 100-(15*h) 0.05 0 0.25 0.5 0.75 1]; 
end
%%

for m = 1:6


LD50_dist(m,1) = params2(m,1)- params2(m,3);
init_frac(m,1) = params2(5);
model_sim = @(p)two_pop_model(p, Vmaxbymix, dose);
model = model_sim(params2(m,:));
simmodel = repmat(model,1,3);
simdose = repmat(dose,1,3);


% Start with initial parameter distribution and noise from before:
modellong = [];
for i = 1:nsamp
    modellong = vertcat(modellong, simmodel(i,:)');
end

% Add noise and fit in a loop, output parameters
eta =  0.3;
simmodelong = [];
Psim = [];
Plosim = [];
Phisim = [];


for k = 1:100
    noise = eta.*(0.5-rand(1, length(modellong)));
    simmodellong = modellong + noise';
    
    for i = 1:length(simmodellong)
        if simmodellong(i) >1
            simmodellong(i) = 1;
        end
        if simmodellong(i)< 0
            simmodellong(i) = 0;
        end
    end
    
    ind = find(doselong == 0);
    Vmax = simdatalong(ind);
    Vmaxbymixsim = [];
    Vmaxallsim = [];
    for j = 1:nsamp
        Vmaxbymixsim(j) = mean(Vmax(3*j-2:3*j));
        Vmaxmat = repmat(Vmaxbymixsim(j), ns.*nreps,1);
        Vmaxallsim = vertcat(Vmaxallsim, Vmaxmat);
    end
    
    [Phisim(:,k)] = lsqnonlin(@fitsinglepop,...
    psingle0,...
    paramslb(1:2),...
    paramsub(1:2),...
    options,...
    doselong(ihi),...
    simmodellong(ihi));
[Plosim(:,k)] = lsqnonlin(@fitsinglepop,...
    psingle0,...
    paramslb(1:2),...
    paramsub(1:2),...
    options,...
    doselong(ilo),...
    simmodellong(ilo));

[Psim(:,k)] = lsqnonlin(@fitmixedpops,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    doselong,...
    simmodellong,...
    nsamp,...
    Vmaxallsim);

if Psim(1,k) < Psim(3,k)
    Psimold = Psim(:,k);
    Psim(1,k)= Psimold(3);
    Psim(2,k) = Psimold(4);
    Psim(3,k)= Psimold(1);
    Psim(4,k) = Psimold(2);
end
  percent_error(:,k) =abs(100.*( Psim(:,k)-params(j,:)')./params(j,:)');


end

mean_param_val = mean(Psim,2);
std_param_val = std(Psim,1,2);
mean_percent_error = mean(percent_error,2);

lower_param_val = mean_param_val - 1.96*std_param_val;
upper_param_val = mean_param_val + 1.96*std_param_val;

    for j = 1:length(mean_param_val)
        sum_err = 0;
        for i = 1:length(Psim)
            err = (params(j) - Psim(j, i))^2;
            sum_err = sum_err + err;
        end
        MSE(j) = sum_err/(length(Psim)-length(params)); % mean squared error
        CV(j) = std_param_val(j)/(mean_param_val(j));

    
    end

% For this parameter distribution, find the MSE in LD50 estimates and fres
% estimates

MSE_LD50(m,1:2) = horzcat(MSE(1), MSE(3));
MSE_all_params(m,:) = MSE;
MSE_LD50_avg(m) = mean(MSE_LD50(m,:));
STD_LD50(m,1:2) = horzcat(std_param_val(1), std_param_val(3));
STD_LD50_avg(m,1:2) = mean(STD_LD50(m,:));
MSE_init_frac(m) = MSE(5);

mean_percent_error_LD50res(m) = mean_percent_error(1);
mean_percent_error_LD50sens(m) = mean_percent_error(3);

end
%% Visual Displays
% Effect on all parameters for perturbation of a single parameter
Params = {'LD50_{res}', 'slope_{res}', 'LD50_{sens}', 'slope_{sens}', 'f_{res1}', 'f_{res2}', 'f_{res3}', 'f_{res4}', 'f_{res5}'}
MSE_norm2 = MSE_all_params./params;

dist_rmins2 = params2(:,1)-params2(:,3);
MSE_res2 = MSE_norm2(:,1);
MSE_sens2 = MSE_norm2(:,3);

figure (9)
set(gca,'LineWidth',1.5,'FontSize',14)
hold off
plot(dist_rmins2, MSE_res2, 'ro-', 'LineWidth', 2)
text(dist_rmins2(3), MSE_res2(3),'MSE LD50_{res}', 'VerticalAlignment', 'bottom')
hold on
plot(dist_rmins2, MSE_sens2, 'bo-', 'LineWidth',2)
text(dist_rmins2(3), MSE_sens2(3),'MSE LD50_{sens}', 'VerticalAlignment', 'bottom')
xlabel('LD50_{res} - LD50_{sens}')
ylabel('Normalized Mean-Squared Error in LD50 estimation')
title('Error vs. Distance between LD50s as LD50_{sens} decreases')

%%

figure(8)
hold off
set(gca,'LineWidth',1.5,'FontSize',14)
    for k =1:9 % 4%size(MSE_all_params,2) % for all 6 diferent param perturbations
        semilogy(params(:,3), MSE_norm2(:,k),'o-', 'LineWidth', 2);
        text(params(4,3), MSE_norm2(4,k), ['MSE in ', Params{k}], 'HorizontalAlignment','left','VerticalAlignment','top')
        hold on
    end
xlabel('LD50_{sens} uM')
ylabel( 'Normalized Mean-Squared Error in Parameter')
title('Error in Parameters as a Function of LD50_{sens}')



%%
%Display MSE in parameter estimation as a function of distance between LD50s
figure(8)
hold off
semilogy(LD50_dist, MSE_LD50(:,1)./params(:,1), 'ro-', 'LineWidth',2)
hold on 
semilogy(LD50_dist, MSE_LD50(:,2)./params(:,3), 'bo-', 'LineWidth',2)
%plot(LD50_dist, MSE_LD50_avg, 'ko')
legend( 'resistant LD50 MSE', 'sensitive LD50 MSE')
legend boxoff
title ('Error in LD50 identifiability as a function of distance between LD50s')
ylim([ 0 1e5])
ylabel('Normalized MSE in LD50 parameter estimation')
xlabel('LD50_{res}-LD50_{sens} (uM)')
%%
figure(9)
subplot(1,2,1)
hold off
plot(params(:,1), MSE_LD50(:,1), 'r*')
title ('MSE in parameter estimation of LD50')
ylabel('MSE in LD50_{res} parameter estimation')
xlabel('LD50_{res}(uM)')
subplot(1,2,2)
hold off
plot(params(:,3), MSE_LD50(:,2), 'b*')
title ('MSE in parameter estimation of LD50')
ylabel('MSE in LD50_{sens} parameter estimation')
xlabel('LD50_{sens}(uM)')

figure(10)
hold off
plot(LD50_dist, STD_LD50(:,1), 'r*')
hold on 
plot(LD50_dist, STD_LD50(:,2), 'b*')
plot(LD50_dist, STD_LD50_avg, 'ko')
legend( 'resistant LD50 STD', 'sensitive LD50 STD')
legend boxoff
title ('STD in parameter estimation of LD50s as a function of their difference')
ylabel('STD in LD50 parameter estimation')
xlabel('LD50_{res}-LD50_{sens} (uM)')

figure(11)
subplot(1,2,1)
hold off
plot(params(:,1), STD_LD50(:,1), 'r*')
title ('STD in parameter estimation of LD50')
ylabel('STD in LD50_{res} parameter estimation')
xlabel('LD50_{res}(uM)')
subplot(1,2,2)
hold off
plot(params(:,3), STD_LD50(:,2), 'b*')
title ('STD in parameter estimation of LD50')
ylabel('STD in LD50_{sens} parameter estimation')
xlabel('LD50_{sens}(uM)')


figure(12)
hold off
plot(params(:,1), mean_percent_error_LD50res, 'r*')
hold on
plot(params(:,3), mean_percent_error_LD50sens,'b*')
plot(LD50_dist, mean(horzcat(mean_percent_error_LD50res, mean_percent_error_LD50sens)), 'ko')
