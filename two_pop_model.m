function [ model_sim] = two_pop_model( P, Vmaxbymix, D)
%This function takes model parameters and outputs the results for a two
%population model
% Here the first 4 params are LD50 and slope of the resistant and sensitive
% populations respectively, and the rest of the parameters are the fraction
% of sensitive cells

nsamp = length(P) -4;
%   Detailed explanation goes here

for i = 1:nsamp
% model_sim(i,:) = Vmaxbymix(i).*(((P(4+i)./( 1 + exp(P(2).*(D - P(1))))) + ((1-P(4+i))./(1 + exp(P(4).*(D - P(3)))))));
model_sim(i,:) = Vmaxbymix(i).*(((P(4+i)./( 1 + exp(P(2).*(D - P(1))))) + ((1-P(4+i)))./(1 + exp(P(4).*(D - P(3))))));
end
end

