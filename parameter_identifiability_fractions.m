close all; clear all; clc
% Find variability as a function of dose
pat = xlsread('../data/subpop_sorting_Cayman.xls');
% Find variability as a function of dose
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

params = [ 185 0.03 35 0.06 0 0.7 0.8 0.9 1]; % high 
%params = [ 185 0.03 35 0.06 0 0.1 0.2 0.3 1]; % low
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

subplot(1,2,1)
hold off
for i = 1:nsamp
    plot(dose, model(i,:), Color{i},'LineWidth',3)
    hold on
end
for i = 1:nsamp
    plot(simdose, simdata(i,:), '*', 'color', Color{i})
    %text(400, model(i, 4)+0.1, ['f_{res}=, ',num2str(round(params(4+i),2)), ', f_{sens}=, ', num2str(round(1-params(4+i),2))], 'HorizontalAlignment','left','VerticalAlignment','top','color', Color{i})
end
legend('0% ADR', '25% ADR','50% ADR','75% ADR','100% ADR')
title('a', 'FontSize', 18)
xlabel( 'Dose (\muM)')
ylabel('Viability')
% legend( '0% ADR_{model}','0% ADR_{meas}', '25% ADR_{model}','25% ADR_{meas}' ,'50% ADR_{model}','50% ADR_{meas}',  '75% ADR_{model}', '75% ADR_{meas}','100% ADR_{model}', '100% ADR_{meas}')
legend boxoff
%title('Two Population Model Mixtures: Low Resistance Simulated Data Set, \eta = f(dose)_{exp}', 'FontSize', 10)

% Find variability as a function of dose
pat = xlsread('../data/subpop_sorting_Cayman.xls');
% Find variability as a function of dose
dose = pat(:,2);
viability = pat(:,6);




subplot(1,2,2)
hold off
plot(uniq_dose, variability, 'go-', 'LineWidth', 3)
title('b', 'FontSize', 18)
xlabel ('Dose (\muM)')
ylabel('Standard deviation of viability')
%title('Standard deviation in viability as a function of dose')

%%

doselong = dose;
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
ylim([0 .1])
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
%% run to get Psiml
Psiml = Psim;
%% Run to get Psimh
Psimh = Psim;
%%
Color = {'b'; 'c';'g'; 'm'; 'r'};
colorset= varycolor(13)
Colorl =  [colorset(1:4, :); colorset(10, :)];
Colorh =  [colorset(1, :); colorset(7:10, :)];
%%
figure(5)
hold off
subplot(2,1,1)
hold off
for i = 5:9
histogram(Psiml(i,:),'BinWidth',0.02, 'FaceColor', Colorl(i-4, :))
hold on
end
xlabel ('Parameter estimation of fraction resistant')
ylabel('Frequency')
legend('0% ADR', '10% ADR','20% ADR','30% ADR','100% ADR')
title('a                        ', 'FontSize', 18)
xlim ([0 1])
subplot(2,1,2)
hold off
for i = 5:9
histogram(Psimh(i,:),'BinWidth',0.02, 'FaceColor', Colorh(i-4, :))
hold on
end
xlabel ('Parameter estimation of fraction resistant')
ylabel('Frequency')
legend('0% ADR', '70% ADR','80% ADR','90% ADR','100% ADR')
xlim ([0 1])
title('b                        ', 'FontSize', 18)


%% MSE

figure(3)

bar(MSE(5:9))
ylim ([ 0 6e-3])
ylabel('Mean-Squared Error')
title('MSE for each fraction estimate')
%% Determine if parameter distributions are statistically signficantly different
[ h, p]= kstest2(Psiml(6,:), Psim(7,:))


Params = {Psim(5,:)', Psim(6,:)', Psim(7,:)'};
Fracs = {1, 2, 3};

[pl,t, statsl] = anova2(Psiml(5:9,:)',5,'off')

[ph,t, statsh] = anova2(Psimh(5:9,:)',5,'off')

cl= multcompare(statsl)
ch= multcompare(statsh)

