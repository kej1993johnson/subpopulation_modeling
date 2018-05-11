% Subpopulation sorting modeling
% This script loads in the data from the dose-response experiments of
% mixed resistant and naive MCF-7 breast cancer cells and fits them to a
% single, and double population model in which the model performs parameter
% estimation of the center and slope of each population and the proportion
% of cells in each subpopulation. These will be used as model validation to
% compare to  the known experimental LD50s and subpopulation proportions

clear all, close all, clc

%pat = xlsread('../data/subpop_sorting_90_10.xls'); % replace 75 ADR 25 WT with 90:10
%pat = xlsread('../data/subpop_sorting_oldADR.xls');
pat = xlsread('../data/subpop_sorting_Cayman.xls');

%%
ns = 12;
nreps = 3;
% for replicates from 6/30/17 use nsamp = 4
% for replicates including both (6/30/17 + 8/4/17) use nsamp = 7, will
% change this when Grant reruns all mixtures to be nsamp = 12
nsamp = 5;
ind = pat(:,7)<= nsamp;

dose = pat(ind,2);
viability = pat(ind,6);

ind = find(dose == 0);
Vmax = viability(ind);
Vmaxbymix = [];
Vmaxall =[];
for j = 1:nsamp
    Vmaxbymix(j) = mean(Vmax(3*j-2:3*j));
    Vmaxmat = repmat(Vmaxbymix(j), ns.*nreps,1);
    Vmaxall = vertcat(Vmaxall, Vmaxmat);
end

%%
params0 = horzcat( [200 0.01 25 0.01], 0.5.*ones(1, nsamp));%  LD50res, sloperes, LD50sens, slopesens, fres
paramslb = zeros( 1, 4+nsamp);
paramsub = horzcat( [ Inf 1 Inf 1], ones(1, nsamp));
options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
               

    ihi = find(pat(:,7) ==5); % 1st sample is pure ADR
    ilo = find(pat(:,7) == 1); % 4th sample is pure WT
    psingle0 = [ 70 0.1];

[Phi, resnormhi, residualshi] = lsqnonlin(@fitsinglepop,...
    psingle0,...
    paramslb(1:2),...
    paramsub(1:2),...
    options,...
    dose(ihi),...
    viability(ihi),...
    Vmaxall(ihi));
[Plo, resnorm, residualslo] = lsqnonlin(@fitsinglepop,...
    psingle0,...
    paramslb(1:2),...
    paramsub(1:2),...
    options,...
    dose(ilo),...
    viability(ilo),...
    Vmaxall(ilo));

[P, resnorm, residuals] = lsqnonlin(@fitmixedpops,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    viability,...
    nsamp,...
    Vmaxall);

D = 1:1:max(dose);

% Find CI on parameter estimates

[lowerlim, upperlim] = finderrorBS(residuals,P, dose, nsamp, nreps, Vmaxbymix, Vmaxall);
%%
[lowerlimlo, upperlimlo, betaboot] = finderrorBSsingle(residualslo, Plo, dose(ilo), 1, nreps, Vmaxbymix(1), Vmaxall(ilo))
[lowerlimhi, upperlimhi, betaboothi] = finderrorBSsingle(residualshi, Phi, dose(ihi), 1, nreps, Vmaxbymix(5), Vmaxall(ihi))
%%
% find thwaw for single pop


%Pactual = [ 0.07 0.28 0.51 0.75 0.96]; old ADR old dox

Pactual = [0.07 0.30 0.54 0.83 1]; % with Cayman dox(% ADR)


stddevs = [ 0.0054 0.0227 0.0258 0.0264 0.0127];
CImeas = stddevs*1.96;
%%
CImodel = (upperlim-lowerlim)./2
CImeasres = (upperlimhi-lowerlimhi)./2
CImeassens= (upperlimlo - lowerlimlo)./2

%%
switch nsamp
    case 4 %6/30/17 samples
        Params = [ Phi, Plo, P];
        param_table = dataset({Params,'LD50ADR', 'slopeADR', 'LD50WT', 'slopeWT', 'LD50_res', 'slope_res', 'LD50_sens', 'slope_sens', 'fres1', 'fres2', 'fres3','fres4'});
        save('../out/param_table4.mat', 'param_table')
        case 5 %
        Params = [ P];
        Meas = [Phi, Plo, Pactual];
        param_table = dataset({Params, 'LD50_res', 'slope_res', 'LD50_sens', 'slope_sens', 'fres1', 'fres2', 'fres3','fres4', 'fres5'});
        save('../out/param_table4.mat', 'param_table')
        meas_table = dataset({Meas, 'LD50ADR', 'slopeADR', 'LD50WT', 'slopeWT','fres1', 'fres2', 'fres3','fres4', 'fres5'});
    case 7 % 6/30/17 samples + 8/4/17 samples
        param_table = dataset({P, 'LD50_res', 'slope_res', 'LD50_sens', 'slope_sens', 'fres1', 'fres2', 'fres3','fres4', 'fres5', 'fres6', 'fres7'});
        save('../out/param_table7.mat', 'param_table')

end

%% Plot the results against the raw data

% Find maximum viability for each dose-response curve

%[orderedP, ind]= sort(P(5:9)); % puts the fraction parameters in order from lowest to highest
%
Color = {'b'; 'c';'g'; 'm'; 'r'};
hold off
n = length(viability);
ns = n/nsamp;
figure(1)
hold off
subplot(1,3,1:2)
set(gca,'LineWidth',1.2,'FontSize',10)
for i = 1:nsamp
model(i,:) = Vmaxbymix(i).*(((P(4+i)./( 1 + exp(P(2).*(D - P(1))))) + ((1-P(4+i))./(1 + exp(P(4).*(D - P(3)))))));

plot(D, model(i,:), Color{i},'LineWidth',3)
hold on
end
%text(600, model(i, 4), ['f_{res}=, ',num2str(round(P(4+i),2)), ', f_{sens}=, ', num2str(round(1-P(4+i),2))], 'HorizontalAlignment','left','VerticalAlignment','top','color', Color{i})
hold on
for i = 1:nsamp
scatter(dose(ns*(i-1)+1:ns*i), viability(ns*(i-1)+1:ns*i), Color{i}, 'LineWidth', 1)
end
%title ('Subpopulation Sorting')
xlabel( 'Dose (\muM)')
ylabel('Viability')

legend( '0% ADR', '25% ADR','50% ADR',  '75% ADR', '100% ADR')
legend boxoff
title('a', 'FontSize', 18)
Pactual = [0.07 0.31 0.54 0.83 1];

stddevs = [ 0.01 0.02 0.01 0.01 0.02];
CImeas = stddevs*1.96;

errorparams = upperlim-lowerlim;
errorfracs = errorparams(5:9);

% Find R-squared value 

SStot = sum((Pactual - mean(Pactual)).^2)
SSres = sum((Pactual - P(5:end)).^2)
Rsq = 1-SSres/SStot


x = 1:5:100;
y = x;
subplot(1,3,3)
set(gca,'LineWidth',1.2,'FontSize',10)
hold off
for i = 1:5
errorbar(100.*Pactual(i), 100.*P(4+i), 100.*errorfracs(i)/2,'o','Color', Color{i}, 'LineWidth',1.2)
hold on
end
text(50, 50, ['R-sq =' num2str(round(Rsq, 3))])
hold on
plot(x,y, 'k--')
xlabel('Measured Percent ADR')
ylabel('Model Estimated Percent ADR')
title ('b', 'FontSize', 18)
legend('0% ADR', '25% ADR','50% ADR','75% ADR','100% ADR', 'line of unity')
legend boxoff
title ('b', 'FontSize', 18)

%% Find R-squared value 

SStot = sum((Pactual - mean(Pactual)).^2)
SSres = sum((Pactual - P(5:end)).^2)
Rsq = 1-SSres/SStot






%% Figure of 95% CI overlaid on data
Color = {'b'; 'r'; 'k'; 'g'; 'm'};
figure;
hold off
for i = 1:nsamp
model(i,:) = (Vmaxbymix(i).*((P(4+i)./( 1 + exp(P(2).*(D - P(1))))) + ((1-P(4+i))./(1 + exp(P(4).*(D - P(3)))))));
model_low (i,:) = (Vmaxbymix(i).*((lowerlim(4+i)./( 1 + exp(P(2).*(D - lowerlim(1))))) + ((1-lowerlim(4+i))./(1 + exp(P(4).*(D - lowerlim(3)))))));
model_high (i,:) = (Vmaxbymix(i).*((upperlim(4+i)./( 1 + exp(P(2).*(D - upperlim(1))))) + ((1-upperlim(4+i))./(1 + exp(P(4).*(D - upperlim(3)))))));

plot(D, model(i,:), Color{i},'LineWidth',3)
hold on
plot(D, model_low(i,:), Color{i},'LineWidth',1)
plot(D, model_high(i,:), Color{i},'LineWidth',1)

scatter(dose(ns*(i-1)+1:ns*i), viability(ns*(i-1)+1:ns*i), Color{i}, 'LineWidth', 1)
title (['Subpopulation Sorting Sample ', num2str(i)])


end
%%
switch nsamp
    case 4
comb_tbl = vertcat(lowerlim, P, upperlim);
param_table_CI = dataset({comb_tbl, 'LD50_res', 'slope_res', 'LD50_sens', 'slope_sens', 'fres1', 'fres2', 'fres3','fres4'});
        save('../out/param_table_CI4.mat', 'param_table_CI')
end



