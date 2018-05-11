function[lowerlim, upperlim] = finderrorBS_ADRenr(residuals, P, dose, nsamp, nreps, Vmaxbymix, Vmaxall)

% Find model values at each dose
for i = 1:nsamp
model_val(i,:) = Vmaxbymix(i).*((((P(4+i)./( 1 + exp(P(2).*(dose(1:12) - P(1))))) + ((1-P(4+i))./(1 + exp(P(4).*(dose(1:12) - P(3))))))));
end
model_val(6,:)= model_val(5,:);
model_vallong = repmat(model_val, 1, nreps);

model_new= reshape(model_vallong', [length(dose), 1]); % puts into single vector


nboot = 500;
[~, bootIndices] = bootstrp(nboot, [], residuals); % randomly generates indices
bootResiduals = residuals(bootIndices); % uses indices to sample from residuals with replacement
viaBoot = repmat(model_new,1,nboot) + bootResiduals; % creates simulated data sets with randomly added residuals
for j = 1: size(viaBoot,1)
    for k = 1:size(viaBoot,2)
        if viaBoot(j,k) > 1 
            viaBoot(j,k) = 1;
        end
         if viaBoot(j,k) < 0
             viaBoot(j,k) = 0;
         end
             
    end
end
% build up the bootstrap data set
options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
params0 = horzcat( [180 0.1 35 0.1], 0.5.*ones(1, nsamp));%  LD50res, sloperes, LD50sens, slopesens, fres
paramslb = zeros( 1, 4+nsamp);
paramsub = horzcat( [ Inf 1 Inf 1], ones(1, nsamp));
 for i = 1:nboot
     viasim = viaBoot(:,i);
%    switch nsamp
%          case 4
        betaboot(:,i) = lsqnonlin(@fitmixedpops_ADRenr,...
        params0,...
        paramslb,...
        paramsub,...
        options,...
        dose,...
        viasim,...
        nsamp,...
        Vmaxall);
%     end 

end

 bootCI = prctile(betaboot', [2.5 97.5]);
 lowerlim = bootCI(1,:);
 upperlim = bootCI(2,:);
end