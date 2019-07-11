function [fluxes, rxnIdx] = getMetProduction(model,metName,solVector,nonZero)
% Which fluxes produce/consume a metabolite (excluding transport
% reactions).
% nonZero   only include non zero fluxes, default false

metLoc          = find(ismember(model.metNames,metName));
[metIdx,rxnIdx] = find(model.S(metLoc,:));
coeff           = full(model.S(metLoc(metIdx),rxnIdx));

for i=1:size(solVector,2)
    fluxes(i,:) = sum(solVector(rxnIdx,i) .* transpose(coeff),2);
end

if nonZero
    nonZero     = sum(fluxes,1)~=0;
    fluxes      = transpose(fluxes(:,nonZero));
    rxnIdx      = rxnIdx(nonZero);
else
    fluxes      = transpose(fluxes);
end
end