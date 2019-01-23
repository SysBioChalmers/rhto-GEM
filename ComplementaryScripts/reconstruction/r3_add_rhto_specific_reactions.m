%% Rhodosporidium specific reactions
load('../../scrap/model_r2.mat','model');

metsToAdd.metNames={'15-cis-phytoene';'3,4-dehydrolycopene';'O-docosanoylcarnitine';'O-docosanoylcarnitine';'O-hexacosanoylcarnitine';'O-hexacosanoylcarnitine';'O-icosanoylcarnitine';'O-icosanoylcarnitine';'O-lauroylcarnitine';'O-lauroylcarnitine';'O-myristoylcarnitine';'O-myristoylcarnitine';'O-oleoylcarnitine';'O-oleoylcarnitine';'O-palmitoleoylcarnitine';'O-palmitoleoylcarnitine';'O-palmitoylcarnitine';'O-palmitoylcarnitine';'O-stearoylcarnitine';'O-stearoylcarnitine';'O-tetracosanoylcarnitine';'O-tetracosanoylcarnitine';'beta-carotene';'docosanoyl-CoA';'gamma-carotene';'hexacosanoyl-CoA';'icosanoyl-CoA';'lauroyl-CoA';'lycopene';'myristoyl-CoA';'neurosporene';'oleoyl-CoA';'palmitoleoyl-CoA(4-)';'palmitoyl-CoA';'stearoyl-CoA';'tetracosanoyl-CoA';'torulene';'torularhodin'};
metsToAdd.compartments={'c';'c';'c';'m';'c';'m';'c';'m';'c';'m';'c';'m';'c';'m';'c';'m';'c';'m';'c';'m';'c';'m';'c';'m';'c';'m';'m';'m';'c';'m';'c';'m';'m';'m';'m';'m';'c';'c'};
metsToAdd.mets=generateNewIds(model,'mets','st_',length(metsToAdd.metNames));
%metsToAdd.metFormulas
model=addMets(model,metsToAdd); clear metsToAdd;

rxnsToAdd.equations={'geranylgeranyl diphosphate[c] <=> diphosphate[c] + 15-cis-phytoene[c]';'3 NAD[c] + 15-cis-phytoene[c] <=> 3 H+[c] + 3 NADH[c] + neurosporene[c]';'H+[c] + NADH[c] + oxygen[c] + neurosporene[c] <=> 2 H2O[c] + NAD[c] + lycopene[c]';'lycopene[c] <=> gamma-carotene[c]';'gamma-carotene[c] <=> beta-carotene[c]';'H+[c] + NADH[c] + oxygen[c] + lycopene[c] <=> 2 H2O[c] + NAD[c] + 3,4-dehydrolycopene[c]';'3,4-dehydrolycopene[c] <=> torulene[c]';'oxygen[c] + torulene[c] <=> torularhodin[c]';'(R)-carnitine[c] + palmitoyl-CoA[c] => coenzyme A[c] + O-palmitoylcarnitine[c]';'(R)-carnitine[m] + O-palmitoylcarnitine[c] <=> (R)-carnitine[c] + O-palmitoylcarnitine[m]';'coenzyme A[m] + O-palmitoylcarnitine[m] => (R)-carnitine[m] + palmitoyl-CoA[m]';'(R)-carnitine[c] + icosanoyl-CoA[c] => coenzyme A[c] + O-icosanoylcarnitine[c]';'(R)-carnitine[m] + O-icosanoylcarnitine[c] <=> (R)-carnitine[c] + O-icosanoylcarnitine[m]';'coenzyme A[m] + O-icosanoylcarnitine[m] => (R)-carnitine[m] + icosanoyl-CoA[m]';'(R)-carnitine[c] + docosanoyl-CoA[c] => coenzyme A[c] + O-docosanoylcarnitine[c]';'(R)-carnitine[m] + O-docosanoylcarnitine[c] <=> (R)-carnitine[c] + O-docosanoylcarnitine[m]';'coenzyme A[m] + O-docosanoylcarnitine[m] => (R)-carnitine[m] + docosanoyl-CoA[m]';'(R)-carnitine[c] + tetracosanoyl-CoA[c] => coenzyme A[c] + O-tetracosanoylcarnitine[c]';'(R)-carnitine[m] + O-tetracosanoylcarnitine[c] <=> (R)-carnitine[c] + O-tetracosanoylcarnitine[m]';'coenzyme A[m] + O-tetracosanoylcarnitine[m] => (R)-carnitine[m] + tetracosanoyl-CoA[m]';'(R)-carnitine[c] + hexacosanoyl-CoA[c] => coenzyme A[c] + O-hexacosanoylcarnitine[c]';'(R)-carnitine[m] + O-hexacosanoylcarnitine[c] <=> (R)-carnitine[c] + O-hexacosanoylcarnitine[m]';'coenzyme A[m] + O-hexacosanoylcarnitine[m] => (R)-carnitine[m] + hexacosanoyl-CoA[m]';'(R)-carnitine[c] + oleoyl-CoA[c] => coenzyme A[c] + O-oleoylcarnitine[c]';'(R)-carnitine[m] + O-oleoylcarnitine[c] <=> (R)-carnitine[c] + O-oleoylcarnitine[m]';'coenzyme A[m] + O-oleoylcarnitine[m] => (R)-carnitine[m] + oleoyl-CoA[m]';'(R)-carnitine[c] + palmitoleoyl-CoA(4-)[c] => coenzyme A[c] + O-palmitoleoylcarnitine[c]';'(R)-carnitine[m] + O-palmitoleoylcarnitine[c] <=> (R)-carnitine[c] + O-palmitoleoylcarnitine[m]';'coenzyme A[m] + O-palmitoleoylcarnitine[m] => (R)-carnitine[m] + palmitoleoyl-CoA(4-)[m]';'(R)-carnitine[c] + stearoyl-CoA[c] => coenzyme A[c] + O-stearoylcarnitine[c]';'(R)-carnitine[m] + O-stearoylcarnitine[c] <=> (R)-carnitine[c] + O-stearoylcarnitine[m]';'coenzyme A[m] + O-stearoylcarnitine[m] => (R)-carnitine[m] + stearoyl-CoA[m]';'(R)-carnitine[c] + lauroyl-CoA[c] => coenzyme A[c] + O-lauroylcarnitine[c]';'(R)-carnitine[m] + O-lauroylcarnitine[c] <=> (R)-carnitine[c] + O-lauroylcarnitine[m]';'coenzyme A[m] + O-lauroylcarnitine[m] => (R)-carnitine[m] + lauroyl-CoA[m]';'(R)-carnitine[c] + myristoyl-CoA[c] => coenzyme A[c] + O-myristoylcarnitine[c]';'(R)-carnitine[m] + O-myristoylcarnitine[c] <=> (R)-carnitine[c] + O-myristoylcarnitine[m]';'coenzyme A[m] + O-myristoylcarnitine[m] => (R)-carnitine[m] + myristoyl-CoA[m]';'(S)-malate[c] + NADP(+)[c] => carbon dioxide[c] + NADPH[c] + pyruvate[c]'};
rxnsToAdd.rxns=generateNewIds(model,'rxns','t_',length(rxnsToAdd.equations));
rxnsToAdd.rxnNames={'phytoene synthase';'phytoene dehydrogenase';'phytoene dehydrogenase';'phytoene synthase';'phytoene synthase';'phytoene dehydrogenase';'phytoene synthase';'torularhodin synthase';'cytoplasmatic carnitine acyltransferase (C16:0)';'fatty acylcarnitine transport (C16:0)';'mitochondrial carnitine acyltransferase (C16:0)';'cytoplasmatic carnitine acyltransferase (C20:0)';'fatty acylcarnitine transport (C20:0)';'mitochondrial carnitine acyltransferase (C20:0)';'cytoplasmatic carnitine acyltransferase (C22:0)';'fatty acylcarnitine transport (C22:0)';'mitochondrial carnitine acyltransferase (C22:0)';'cytoplasmatic carnitine acyltransferase (C24:0)';'fatty acylcarnitine transport (C24:0)';'mitochondrial carnitine acyltransferase (C24:0)';'cytoplasmatic carnitine acyltransferase (C26:0)';'fatty acylcarnitine transport (C26:0)';'mitochondrial carnitine acyltransferase (C26:0)';'cytoplasmatic carnitine acyltransferase (C18:1)';'fatty acylcarnitine transport (C18:1)';'mitochondrial carnitine acyltransferase (C18:1)';'cytoplasmatic carnitine acyltransferase (C16:1)';'fatty acylcarnitine transport (C16:1)';'mitochondrial carnitine acyltransferase (C16:1)';'cytoplasmatic carnitine acyltransferase (C18:0)';'fatty acylcarnitine transport (C18:0)';'mitochondrial carnitine acyltransferase (C18:0)';'cytoplasmatic carnitine acyltransferase (C12:0)';'fatty acylcarnitine transport (C12:0)';'mitochondrial carnitine acyltransferase (C12:0)';'cytoplasmatic carnitine acyltransferase (C14:0)';'fatty acylcarnitine transport (C14:0)';'mitochondrial carnitine acyltransferase (C14:0)';'malic enzyme (NADP) cytoplasmic'};
rxnsToAdd.grRules={'RHTO_04605';'RHTO_04602';'RHTO_04602';'RHTO_04605';'RHTO_04605';'RHTO_04602';'RHTO_04605';'';'RHTO_01903';'RHTO_01354';'RHTO_01064 and RHTO_01341 and RHTO_00095';'RHTO_01903';'RHTO_01354';'RHTO_01064 and RHTO_01341 and RHTO_00095';'RHTO_01903';'RHTO_01354';'RHTO_01064 and RHTO_01341 and RHTO_00095';'RHTO_01903';'RHTO_01354';'RHTO_01064 and RHTO_01341 and RHTO_00095';'RHTO_01903';'RHTO_01354';'RHTO_01064 and RHTO_01341 and RHTO_00095';'RHTO_01903';'RHTO_01354';'RHTO_01064 and RHTO_01341 and RHTO_00095';'RHTO_01903';'RHTO_01354';'RHTO_01064 and RHTO_01341 and RHTO_00095';'RHTO_01903';'RHTO_01354';'RHTO_01064 and RHTO_01341 and RHTO_00095';'RHTO_01903';'RHTO_01354';'RHTO_01064 and RHTO_01341 and RHTO_00095';'RHTO_01903';'RHTO_01354';'RHTO_01064 and RHTO_01341 and RHTO_00095';'RHTO_03795'};
model=addRxns(model,rxnsToAdd,3,'',false,true); clear rxnsToAdd

%% Mitochondrial beta-oxidation
metsToAdd.metNames={'(R)-3-hydroxy-cis-dodec-5-enoyl-CoA';'(R)-3-hydroxy-cis-hexadec-7-enoyl-CoA';'(R)-3-hydroxy-cis-hexadec-9-enoyl-CoA';'(R)-3-hydroxy-cis-octadec-9-enoyl-CoA';'(R)-3-hydroxy-cis-tetradec-5-enoyl-CoA';'(R)-3-hydroxy-cis-tetradec-7-enoyl-CoA';'(R)-3-hydroxybutanoyl-CoA';'(R)-3-hydroxydecanoyl-CoA';'(R)-3-hydroxydocosanoyl-CoA';'(R)-3-hydroxyhexanoyl-CoA';'(R)-3-hydroxyicosanoyl-CoA';'(R)-3-hydroxylauroyl-CoA';'(R)-3-hydroxyoctanoyl-CoA';'(R)-3-hydroxytetracosanoyl-CoA';'(S)-3-hydroxyhexacosanoyl-CoA';'(S)-3-hydroxypalmitoyl-CoA';'(S)-3-hydroxytetradecanoyl-CoA';'3-hydroxyoctadecanoyl-CoA';'3-oxo-cis-dodec-5-enoyl-CoA';'3-oxo-cis-hexadec-7-enoyl-CoA';'3-oxo-cis-hexadec-9-enoyl-CoA';'3-oxo-cis-octadec-9-enoyl-CoA';'3-oxo-cis-tetradec-5-enoyl-CoA';'3-oxo-cis-tetradec-7-enoyl-CoA';'3-oxodecanoyl-CoA';'3-oxodocosanoyl-CoA';'3-oxohexacosanoyl-CoA';'3-oxohexanoyl-CoA';'3-oxoicosanoyl-CoA';'3-oxolauroyl-CoA';'3-oxooctadecanoyl-CoA';'3-oxooctanoyl-CoA';'3-oxopalmitoyl-CoA';'3-oxotetracosanoyl-CoA';'3-oxotetradecanoyl-CoA';'butanoyl-CoA';'cis-dec-3-enoyl-CoA';'cis-dodec-3-enoyl-CoA';'cis-dodec-5-enoyl-CoA';'cis-hexadec-7-enoyl-CoA';'cis-tetradec-5-enoyl-CoA';'cis-tetradec-7-enoyl-CoA';'decanoyl-CoA';'hexadec-2-enoyl-CoA';'hexanoyl-CoA';'octanoyl-CoA';'trans-2,cis-5-dodecadienoyl-CoA';'trans-2,cis-5-tetradecadienoyl-CoA';'trans-2,cis-7-hexadecadienoyl-CoA';'trans-2,cis-7-tetradecadienoyl-CoA';'trans-2,cis-9-hexadecadienoyl-CoA';'trans-2,cis-9-octadecadienoyl-CoA';'trans-but-2-enoyl-CoA';'trans-dec-2-enoyl-CoA';'trans-docos-2-enoyl-CoA';'trans-dodec-2-enoyl-CoA';'trans-hex-2-enoyl-CoA';'trans-hexacos-2-enoyl-CoA';'trans-icos-2-enoyl-CoA';'trans-oct-2-enoyl-CoA';'trans-octadec-2-enoyl-CoA';'trans-tetracos-2-enoyl-CoA';'trans-tetradec-2-enoyl-CoA'};
metsToAdd.compartments={'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m';'m'};
metsToAdd.mets=generateNewIds(model,'mets','st_',length(metsToAdd.metNames));
%metsToAdd.metFormulas
model=addMets(model,metsToAdd); clear metsToAdd;

rxnsToAdd.equations={'NAD[m] + (S)-3-hydroxytetradecanoyl-CoA[m] <=> H+[m] + NADH[m] + 3-oxotetradecanoyl-CoA[m]';'coenzyme A[m] + 3-oxooctadecanoyl-CoA[m] <=> acetyl-CoA[m] + palmitoyl-CoA[m]';'coenzyme A[m] + 3-oxohexacosanoyl-CoA[m] <=> acetyl-CoA[m] + tetracosanoyl-CoA[m]';'coenzyme A[m] + 3-oxopalmitoyl-CoA[m] <=> acetyl-CoA[m] + myristoyl-CoA[m]';'coenzyme A[m] + 3-oxotetradecanoyl-CoA[m] <=> acetyl-CoA[m] + lauroyl-CoA[m]';'coenzyme A[m] + 3-oxolauroyl-CoA[m] <=> acetyl-CoA[m] + decanoyl-CoA[m]';'FAD[m] + decanoyl-CoA[m] <=> FADH2[m] + trans-dec-2-enoyl-CoA[m]';'FAD[m] + lauroyl-CoA[m] <=> FADH2[m] + trans-dodec-2-enoyl-CoA[m]';'FAD[m] + hexacosanoyl-CoA[m] <=> FADH2[m] + trans-hexacos-2-enoyl-CoA[m]';'FAD[m] + palmitoyl-CoA[m] <=> FADH2[m] + hexadec-2-enoyl-CoA[m]';'FAD[m] + stearoyl-CoA[m] <=> FADH2[m] + trans-octadec-2-enoyl-CoA[m]';'FAD[m] + myristoyl-CoA[m] <=> FADH2[m] + trans-tetradec-2-enoyl-CoA[m]';'FAD[m] + butanoyl-CoA[m] <=> FADH2[m] + trans-but-2-enoyl-CoA[m]';'FAD[m] + hexanoyl-CoA[m] <=> FADH2[m] + trans-hex-2-enoyl-CoA[m]';'FAD[m] + octanoyl-CoA[m] <=> FADH2[m] + trans-oct-2-enoyl-CoA[m]';'FAD[m] + icosanoyl-CoA[m] <=> FADH2[m] + trans-icos-2-enoyl-CoA[m]';'FAD[m] + docosanoyl-CoA[m] <=> FADH2[m] + trans-docos-2-enoyl-CoA[m]';'FAD[m] + tetracosanoyl-CoA[m] <=> FADH2[m] + trans-tetracos-2-enoyl-CoA[m]';'FAD[m] + palmitoleoyl-CoA(4-)[m] <=> FADH2[m] + trans-2,cis-9-hexadecadienoyl-CoA[m]';'FAD[m] + cis-tetradec-7-enoyl-CoA[m] <=> FADH2[m] + trans-2,cis-7-tetradecadienoyl-CoA[m]';'FAD[m] + cis-dodec-5-enoyl-CoA[m] <=> FADH2[m] + trans-2,cis-5-dodecadienoyl-CoA[m]';'FAD[m] + oleoyl-CoA[m] <=> FADH2[m] + trans-2,cis-9-octadecadienoyl-CoA[m]';'FAD[m] + cis-hexadec-7-enoyl-CoA[m] <=> FADH2[m] + trans-2,cis-7-hexadecadienoyl-CoA[m]';'FAD[m] + cis-tetradec-5-enoyl-CoA[m] <=> FADH2[m] + trans-2,cis-5-tetradecadienoyl-CoA[m]';'H2O[m] + trans-dec-2-enoyl-CoA[m] <=> (R)-3-hydroxydecanoyl-CoA[m]';'H2O[m] + trans-dodec-2-enoyl-CoA[m] <=> (R)-3-hydroxylauroyl-CoA[m]';'H2O[m] + trans-tetradec-2-enoyl-CoA[m] <=> (S)-3-hydroxytetradecanoyl-CoA[m]';'H2O[m] + hexadec-2-enoyl-CoA[m] <=> (S)-3-hydroxypalmitoyl-CoA[m]';'H2O[m] + trans-octadec-2-enoyl-CoA[m] <=> 3-hydroxyoctadecanoyl-CoA[m]';'H2O[m] + trans-hexacos-2-enoyl-CoA[m] <=> (S)-3-hydroxyhexacosanoyl-CoA[m]';'H2O[m] + trans-but-2-enoyl-CoA[m] <=> (R)-3-hydroxybutanoyl-CoA[m]';'H2O[m] + trans-hex-2-enoyl-CoA[m] <=> (R)-3-hydroxyhexanoyl-CoA[m]';'H2O[m] + trans-oct-2-enoyl-CoA[m] <=> (R)-3-hydroxyoctanoyl-CoA[m]';'H2O[m] + trans-icos-2-enoyl-CoA[m] <=> (R)-3-hydroxyicosanoyl-CoA[m]';'H2O[m] + trans-docos-2-enoyl-CoA[m] <=> (R)-3-hydroxydocosanoyl-CoA[m]';'H2O[m] + trans-tetracos-2-enoyl-CoA[m] <=> (R)-3-hydroxytetracosanoyl-CoA[m]';'H2O[m] + trans-2,cis-9-hexadecadienoyl-CoA[m] <=> (R)-3-hydroxy-cis-hexadec-9-enoyl-CoA[m]';'H2O[m] + trans-2,cis-7-tetradecadienoyl-CoA[m] <=> (R)-3-hydroxy-cis-tetradec-7-enoyl-CoA[m]';'H2O[m] + trans-2,cis-5-dodecadienoyl-CoA[m] <=> (R)-3-hydroxy-cis-dodec-5-enoyl-CoA[m]';'H2O[m] + trans-2,cis-9-octadecadienoyl-CoA[m] <=> (R)-3-hydroxy-cis-octadec-9-enoyl-CoA[m]';'H2O[m] + trans-2,cis-7-hexadecadienoyl-CoA[m] <=> (R)-3-hydroxy-cis-hexadec-7-enoyl-CoA[m]';'H2O[m] + trans-2,cis-5-tetradecadienoyl-CoA[m] <=> (R)-3-hydroxy-cis-tetradec-5-enoyl-CoA[m]';'NAD[m] + (R)-3-hydroxydecanoyl-CoA[m] <=> H+[m] + NADH[m] + 3-oxodecanoyl-CoA[m]';'NAD[m] + (R)-3-hydroxylauroyl-CoA[m] <=> H+[m] + NADH[m] + 3-oxolauroyl-CoA[m]';'NAD[m] + (S)-3-hydroxypalmitoyl-CoA[m] <=> H+[m] + NADH[m] + 3-oxopalmitoyl-CoA[m]';'NAD[m] + 3-hydroxyoctadecanoyl-CoA[m] <=> H+[m] + NADH[m] + 3-oxooctadecanoyl-CoA[m]';'NAD[m] + (S)-3-hydroxyhexacosanoyl-CoA[m] <=> H+[m] + NADH[m] + 3-oxohexacosanoyl-CoA[m]';'NAD[m] + (R)-3-hydroxybutanoyl-CoA[m] <=> acetoacetyl-CoA[m] + H+[m] + NADH[m]';'NAD[m] + (R)-3-hydroxyhexanoyl-CoA[m] <=> H+[m] + NADH[m] + 3-oxohexanoyl-CoA[m]';'NAD[m] + (R)-3-hydroxyoctanoyl-CoA[m] <=> H+[m] + NADH[m] + 3-oxooctanoyl-CoA[m]';'NAD[m] + (R)-3-hydroxyicosanoyl-CoA[m] <=> H+[m] + NADH[m] + 3-oxoicosanoyl-CoA[m]';'NAD[m] + (R)-3-hydroxydocosanoyl-CoA[m] <=> H+[m] + NADH[m] + 3-oxodocosanoyl-CoA[m]';'NAD[m] + (R)-3-hydroxytetracosanoyl-CoA[m] <=> H+[m] + NADH[m] + 3-oxotetracosanoyl-CoA[m]';'NAD[m] + (R)-3-hydroxy-cis-hexadec-9-enoyl-CoA[m] <=> H+[m] + NADH[m] + 3-oxo-cis-hexadec-9-enoyl-CoA[m]';'NAD[m] + (R)-3-hydroxy-cis-tetradec-7-enoyl-CoA[m] <=> H+[m] + NADH[m] + 3-oxo-cis-tetradec-7-enoyl-CoA[m]';'NAD[m] + (R)-3-hydroxy-cis-dodec-5-enoyl-CoA[m] <=> H+[m] + NADH[m] + 3-oxo-cis-dodec-5-enoyl-CoA[m]';'NAD[m] + (R)-3-hydroxy-cis-octadec-9-enoyl-CoA[m] <=> H+[m] + NADH[m] + 3-oxo-cis-octadec-9-enoyl-CoA[m]';'NAD[m] + (R)-3-hydroxy-cis-hexadec-7-enoyl-CoA[m] <=> H+[m] + NADH[m] + 3-oxo-cis-hexadec-7-enoyl-CoA[m]';'NAD[m] + (R)-3-hydroxy-cis-tetradec-5-enoyl-CoA[m] <=> H+[m] + NADH[m] + 3-oxo-cis-tetradec-5-enoyl-CoA[m]';'acetoacetyl-CoA[m] + coenzyme A[m] <=> 2 acetyl-CoA[m]';'coenzyme A[m] + 3-oxohexanoyl-CoA[m] <=> acetyl-CoA[m] + butanoyl-CoA[m]';'coenzyme A[m] + 3-oxooctanoyl-CoA[m] <=> acetyl-CoA[m] + hexanoyl-CoA[m]';'coenzyme A[m] + 3-oxoicosanoyl-CoA[m] <=> acetyl-CoA[m] + stearoyl-CoA[m]';'coenzyme A[m] + 3-oxodocosanoyl-CoA[m] <=> acetyl-CoA[m] + icosanoyl-CoA[m]';'coenzyme A[m] + 3-oxotetracosanoyl-CoA[m] <=> acetyl-CoA[m] + docosanoyl-CoA[m]';'coenzyme A[m] + 3-oxo-cis-hexadec-9-enoyl-CoA[m] <=> acetyl-CoA[m] + cis-tetradec-7-enoyl-CoA[m]';'coenzyme A[m] + 3-oxo-cis-tetradec-7-enoyl-CoA[m] <=> acetyl-CoA[m] + cis-dodec-5-enoyl-CoA[m]';'coenzyme A[m] + 3-oxo-cis-dodec-5-enoyl-CoA[m] <=> acetyl-CoA[m] + cis-dec-3-enoyl-CoA[m]';'coenzyme A[m] + 3-oxo-cis-octadec-9-enoyl-CoA[m] <=> acetyl-CoA[m] + cis-hexadec-7-enoyl-CoA[m]';'coenzyme A[m] + 3-oxo-cis-hexadec-7-enoyl-CoA[m] <=> acetyl-CoA[m] + cis-tetradec-5-enoyl-CoA[m]';'coenzyme A[m] + 3-oxo-cis-tetradec-5-enoyl-CoA[m] <=> acetyl-CoA[m] + cis-dodec-3-enoyl-CoA[m]'};
rxnsToAdd.rxns={'mo_0057';'mo_0100';'mo_0101';'mo_0102';'mo_0105';'mo_0107';'mo_0120';'mo_0121';'mo_0122';'mo_0123';'mo_0124';'mo_0125';'mo_2236';'mo_2237';'mo_2238';'mo_2239';'mo_2240';'mo_2241';'mo_2242';'mo_2243';'mo_2244';'mo_2245';'mo_2246';'mo_2247';'mo_2248';'mo_2249';'mo_2250';'mo_2251';'mo_2252';'mo_2253';'mo_2254';'mo_2255';'mo_2256';'mo_2257';'mo_2258';'mo_2259';'mo_2260';'mo_2261';'mo_2262';'mo_2263';'mo_2264';'mo_2265';'mo_2266';'mo_2267';'mo_2268';'mo_2269';'mo_2270';'mo_2271';'mo_2272';'mo_2273';'mo_2274';'mo_2275';'mo_2276';'mo_2277';'mo_2278';'mo_2279';'mo_2280';'mo_2281';'mo_2282';'mo_2283';'mo_2284';'mo_2285';'mo_2286';'mo_2287';'mo_2288';'mo_2289';'mo_2290';'mo_2291';'mo_2292';'mo_2293';'mo_2294'};
rxnsToAdd.rxnNames={'3-hydroxyacyl-CoA dehydrogenase (3-oxotetradecanoyl-CoA)';'acetyl-CoA C-acyltransferase (palmitoyl-CoA)';'acetyl-CoA C-acyltransferase (tetracosanoyl-CoA)';'acetyl-CoA C-acyltransferase (myristoyl-CoA)';'acetyl-CoA C-acyltransferase (lauroyl-CoA)';'acetyl-CoA C-acyltransferase (decanoyl-CoA)';'acyl-CoA oxidase (decanoyl-CoA)';'acyl-CoA oxidase (dodecanoyl-CoA)';'acyl-CoA oxidase (hexacosanoyl-CoA)';'acyl-CoA oxidase (hexadecanoyl-CoA)';'acyl-CoA oxidase (octadecanoyl-CoA)';'acyl-CoA oxidase (tetradecanoyl-CoA)';'acyl-CoA oxidase (butanoyl-CoA)';'acyl-CoA oxidase (hexanoyl-CoA)';'acyl-CoA oxidase (octanoyl-CoA)';'acyl-CoA oxidase (icosanoyl-CoA)';'acyl-CoA oxidase (docosanoyl-CoA)';'acyl-CoA oxidase (tetracosanoyl-CoA)';'acyl-CoA oxidase (palmitoleoyl-CoA)';'acyl-CoA oxidase (cis-tetradec-7-enoyl-CoA)';'acyl-CoA oxidase (cis-dodec-5-enoyl-CoA)';'acyl-CoA oxidase (oleoyl-CoA)';'acyl-CoA oxidase (cis-hexadec-7-enoyl-CoA)';'acyl-CoA oxidase (cis-tetradec-5-enoyl-CoA)';'2-enoyl-CoA hydratase (3-hydroxydecanoyl-CoA)';'2-enoyl-CoA hydratase (3-hydroxydodecanoyl-CoA)';'2-enoyl-CoA hydratase (3-hydroxytetradecanoyl-CoA)';'2-enoyl-CoA hydratase (3-hydroxyhexadecanoyl-CoA)';'2-enoyl-CoA hydratase (3-hydroxyoctadecanoyl-CoA)';'2-enoyl-CoA hydratase (3-hydroxyhexacosanoyl-CoA)';'2-enoyl-CoA hydratase (3-hydroxybutanoyl-CoA)';'2-enoyl-CoA hydratase (3-hydroxyhexanoyl-CoA)';'2-enoyl-CoA hydratase (3-hydroxyoctanoyl-CoA)';'2-enoyl-CoA hydratase (3-hydroxyicosanoyl-CoA)';'2-enoyl-CoA hydratase (3-hydroxydocosanoyl-CoA)';'2-enoyl-CoA hydratase (3-hydroxytetracosanoyl-CoA)';'2-enoyl-CoA hydratase (3-hydroxy-cis-hexadec-9-enoyl-CoA)';'2-enoyl-CoA hydratase (3-hydroxy-cis-tetradec-7-enoyl-CoA)';'2-enoyl-CoA hydratase (3-hydroxy-cis-dodec-5-enoyl-CoA)';'2-enoyl-CoA hydratase (3-hydroxy-cis-octadec-9-enoyl-CoA)';'2-enoyl-CoA hydratase (3-hydroxy-cis-hexadec-7-enoyl-CoA)';'2-enoyl-CoA hydratase (3-hydroxy-cis-tetradec-5-enoyl-CoA)';'3-hydroxyacyl-CoA dehydrogenase (3-oxodecanoyl-CoA)';'3-hydroxyacyl-CoA dehydrogenase (3-oxododecanoyl-CoA)';'3-hydroxyacyl-CoA dehydrogenase (3-oxohexadecanoyl-CoA)';'3-hydroxyacyl-CoA dehydrogenase (3-oxooctadecanoyl-CoA)';'3-hydroxyacyl-CoA dehydrogenase (3-oxohexacosanoyl-CoA)';'3-hydroxyacyl-CoA dehydrogenase (3-oxobutanoyl-CoA)';'3-hydroxyacyl-CoA dehydrogenase (3-oxohexanoyl-CoA)';'3-hydroxyacyl-CoA dehydrogenase (3-oxooctanoyl-CoA)';'3-hydroxyacyl-CoA dehydrogenase (3-oxoicosanoyl-CoA)';'3-hydroxyacyl-CoA dehydrogenase (3-oxodocosanoyl-CoA)';'3-hydroxyacyl-CoA dehydrogenase (3-oxotetracosanoyl-CoA)';'3-hydroxyacyl-CoA dehydrogenase (3-oxo-cis-hexadec-9-enoyl-CoA)';'3-hydroxyacyl-CoA dehydrogenase (3-oxo-cis-tetradec-7-enoyl-CoA)';'3-hydroxyacyl-CoA dehydrogenase (3-oxo-cis-dodec-5-enoyl-CoA)';'3-hydroxyacyl-CoA dehydrogenase (3-oxo-cis-octadec-9-enoyl-CoA)';'3-hydroxyacyl-CoA dehydrogenase (3-oxo-cis-hexadec-7-enoyl-CoA)';'3-hydroxyacyl-CoA dehydrogenase (3-oxo-cis-tetradec-5-enoyl-CoA)';'acetyl-CoA C-acyltransferase (acetyl-CoA)';'acetyl-CoA C-acyltransferase (butanoyl-CoA)';'acetyl-CoA C-acyltransferase (hexanoyl-CoA)';'acetyl-CoA C-acyltransferase (stearoyl-CoA)';'acetyl-CoA C-acyltransferase (icosanoyl-CoA)';'acetyl-CoA C-acyltransferase (docosanoyl-CoA)';'acetyl-CoA C-acyltransferase (cis-tetradec-7-enoyl-CoA)';'acetyl-CoA C-acyltransferase (cis-dodec-5-enoyl-CoA)';'acetyl-CoA C-acyltransferase (cis-dec-3-enoyl-CoA)';'acetyl-CoA C-acyltransferase (cis-hexadec-7-enoyl-CoA)';'acetyl-CoA C-acyltransferase (cis-tetradec-5-enoyl-CoA)';'acetyl-CoA C-acyltransferase (cis-dodec-3-enoyl-CoA)'};
rxnsToAdd.grRules={'RHTO_05537';'RHTO_00476';'RHTO_00477';'RHTO_00478';'RHTO_00479';'RHTO_00480';'RHTO_04971 or RHTO_06738 or RHTO_05797 or RHTO_03594 or RHTO_01625 or RHTO_00397 or RHTO_05970 or RHTO_04579';'RHTO_04971 or RHTO_06738 or RHTO_05797 or RHTO_03594 or RHTO_01625 or RHTO_00397 or RHTO_05970 or RHTO_04580';'RHTO_04971 or RHTO_06738 or RHTO_05797 or RHTO_03594 or RHTO_01625 or RHTO_00397 or RHTO_05970 or RHTO_04581';'RHTO_04971 or RHTO_06738 or RHTO_05797 or RHTO_03594 or RHTO_01625 or RHTO_00397 or RHTO_05970 or RHTO_04582';'RHTO_04971 or RHTO_06738 or RHTO_05797 or RHTO_03594 or RHTO_01625 or RHTO_00397 or RHTO_05970 or RHTO_04583';'RHTO_04971 or RHTO_06738 or RHTO_05797 or RHTO_03594 or RHTO_01625 or RHTO_00397 or RHTO_05970 or RHTO_04584';'RHTO_04971 or RHTO_06738 or RHTO_05797 or RHTO_03594 or RHTO_01625 or RHTO_00397 or RHTO_05970 or RHTO_04585';'RHTO_04971 or RHTO_06738 or RHTO_05797 or RHTO_03594 or RHTO_01625 or RHTO_00397 or RHTO_05970 or RHTO_04586';'RHTO_04971 or RHTO_06738 or RHTO_05797 or RHTO_03594 or RHTO_01625 or RHTO_00397 or RHTO_05970 or RHTO_04587';'RHTO_04971 or RHTO_06738 or RHTO_05797 or RHTO_03594 or RHTO_01625 or RHTO_00397 or RHTO_05970 or RHTO_04588';'RHTO_04971 or RHTO_06738 or RHTO_05797 or RHTO_03594 or RHTO_01625 or RHTO_00397 or RHTO_05970 or RHTO_04589';'RHTO_04971 or RHTO_06738 or RHTO_05797 or RHTO_03594 or RHTO_01625 or RHTO_00397 or RHTO_05970 or RHTO_04590';'RHTO_04971 or RHTO_06738 or RHTO_05797 or RHTO_03594 or RHTO_01625 or RHTO_00397 or RHTO_05970 or RHTO_04591';'RHTO_04971 or RHTO_06738 or RHTO_05797 or RHTO_03594 or RHTO_01625 or RHTO_00397 or RHTO_05970 or RHTO_04592';'RHTO_04971 or RHTO_06738 or RHTO_05797 or RHTO_03594 or RHTO_01625 or RHTO_00397 or RHTO_05970 or RHTO_04593';'RHTO_04971 or RHTO_06738 or RHTO_05797 or RHTO_03594 or RHTO_01625 or RHTO_00397 or RHTO_05970 or RHTO_04594';'RHTO_04971 or RHTO_06738 or RHTO_05797 or RHTO_03594 or RHTO_01625 or RHTO_00397 or RHTO_05970 or RHTO_04595';'RHTO_04971 or RHTO_06738 or RHTO_05797 or RHTO_03594 or RHTO_01625 or RHTO_00397 or RHTO_05970 or RHTO_04596';'RHTO_05407 or RHTO_04298 or RHTO_02517';'RHTO_05407 or RHTO_04298 or RHTO_02518';'RHTO_05407 or RHTO_04298 or RHTO_02519';'RHTO_05407 or RHTO_04298 or RHTO_02520';'RHTO_05407 or RHTO_04298 or RHTO_02521';'RHTO_05407 or RHTO_04298 or RHTO_02522';'RHTO_05407 or RHTO_04298 or RHTO_02523';'RHTO_05407 or RHTO_04298 or RHTO_02524';'RHTO_05407 or RHTO_04298 or RHTO_02525';'RHTO_05407 or RHTO_04298 or RHTO_02526';'RHTO_05407 or RHTO_04298 or RHTO_02527';'RHTO_05407 or RHTO_04298 or RHTO_02528';'RHTO_05407 or RHTO_04298 or RHTO_02529';'RHTO_05407 or RHTO_04298 or RHTO_02530';'RHTO_05407 or RHTO_04298 or RHTO_02531';'RHTO_05407 or RHTO_04298 or RHTO_02532';'RHTO_05407 or RHTO_04298 or RHTO_02533';'RHTO_05407 or RHTO_04298 or RHTO_02534';'RHTO_05520';'RHTO_05521';'RHTO_05522';'RHTO_05523';'RHTO_05524';'RHTO_05525';'RHTO_05526';'RHTO_05527';'RHTO_05528';'RHTO_05529';'RHTO_05530';'RHTO_05531';'RHTO_05532';'RHTO_05533';'RHTO_05534';'RHTO_05535';'RHTO_05536';'RHTO_00481';'RHTO_00482';'RHTO_00483';'RHTO_00484';'RHTO_00485';'RHTO_00486';'RHTO_00487';'RHTO_00488';'RHTO_00489';'RHTO_00490';'RHTO_00491';'RHTO_00492'};
model=addRxns(model,rxnsToAdd,3,'',false,true); clear rxnsToAdd

%% Some lipid related reactions
% More extensive curation is done after gap-filling
metsToAdd.metNames={'linoleoyl-CoA';'linoleoyl-CoA';'linoleoyl-CoA';'octadec-9-ynoyl-CoA'};
metsToAdd.compartments={'c';'erm';'p';'p'};
metsToAdd.mets=generateNewIds(model,'mets','st_',length(metsToAdd.metNames));
%metsToAdd.metFormulas
model=addMets(model,metsToAdd); clear metsToAdd;

rxnsToAdd.equations={'8 coenzyme A[p] + 8 H2O[p] + 8 NAD[p] + 6 oxygen[p] + octadec-9-ynoyl-CoA[p] => 9 acetyl-CoA[p] + 8 H+[p] + 8 NADH[p] + 6 hydrogen peroxide[p]';
    'ATP[c] + citrate[c] + coenzyme A[c] => acetyl-CoA[c] + ADP[c] + oxaloacetate[c] + phosphate[c]';
    'H+[erm] + oxygen[erm] + NADH[erm] + oleoyl-CoA[erm] => 2 H2O[erm] + NAD[erm] + linoleoyl-CoA[erm]';
    'ATP[c] + H2O[c] + linoleoyl-CoA[c] => ADP[c] + H+[c] + phosphate[c] + linoleoyl-CoA[p]'};
rxnsToAdd.rxns=generateNewIds(model,'rxns','t_',length(rxnsToAdd.equations));;
rxnsToAdd.rxnNames={'fatty acid oxidation (C18:2)';
    'ATP-citrate lyase';
    'oleoyl-CoA desaturase (n-C18:1CoA - n-C18:2CoA), ER membrane';
    'fatty acyl-CoA transport via ABC system (C18:2)'};
rxnsToAdd.grRules={'RHTO_03890';'RHTO_03915';'RHTO_03911';'RHTO_06195 and RHTO_04105'};
model=addRxns(model,rxnsToAdd,3,'',false,true); clear rxnsToAdd

model=deleteUnusedGenes(model);

save('../../scrap/model_r3.mat','model');
cd('..'); newCommit(model); cd('reconstruction')