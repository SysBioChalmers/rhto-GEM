function model = setGAM(model,GAM)

xr_pos = getIndexes(model','r_4041','rxns');
for i = 1:length(model.mets)
    S_ix  = model.S(i,xr_pos);
    isGAM = sum(strcmp({'ATP','ADP','H2O','H+','phosphate'},model.metNames{i})) == 1;
    if S_ix ~= 0 && isGAM
        model.S(i,xr_pos) = sign(S_ix) * GAM;
    end
end
end