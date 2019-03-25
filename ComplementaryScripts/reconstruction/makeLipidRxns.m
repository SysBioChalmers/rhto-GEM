function [newEqns, newNames] = makeLipidRxns(templEqn, templName, reduced, comps, chain1, chain2, chain3, chain4)
% makeLipidRxns
%   Given a template reaction equation and name, new equations and names
%   generated using the specified fatty acid chains on the indicated
%   positions, while modifying the compartment.
%
%   templEqn        template reaction equation (see example below)
%   templName       template reaction name (see example below)
%   reduced         logical whether the number of reactions should be
%                   reduced by excluding duplicate reactions that have
%                   similar total chain profile (see example below)
%   comps           string or cell array of compartment abbreviations, such
%                   as 'ce' or {'mm','erm'}
%   chain1          cell array with strings of fatty acid chains for the
%                   first position, e.g. {'16:0','18:1'}
%   chain2          cell array with strings of fatty acid chains for the
%                   second position, if relevant (opt)
%   chain2          cell array with strings of fatty acid chains for the
%                   third position, if relevant (opt)
%   chain2          cell array with strings of fatty acid chains for the
%                   fourth position, if relevant (opt)
%
%   newEqns         cell array of strings with new reaction equations
%   newNames        cell array of strings with new reaction names
%
%   Example of use: assume the following parameters
%   templEqn    = 'CHAIN1_COA[COMP] + gly-3-P[COMP] => CoA[erm] +
%                  1-acyl-gly-3-P (CHAIN1_NUMB)[COMP]'
%   templName   = 'gly-3-p acyltransferase (CHAIN1), COMP'
%   comps       = {'erm','lp'}
%   chain1      = {'16:0','18:3'}
%
%   [newEqns newNames] =
%   makeLipidRxns(templEqn,templName,false,comps,chain)
%
%   will generate four reaction names and equations, two for each
%   COMP compartment and two for each CHAIN fatty acid chain. CHAIN1 to
%   CHAIN4 refer to the fatty acid chain positions (were relevant, in the
%   example there is only one fatty acid chain). _COA specifies that the
%   fatty-acyl-CoA should be used (e.g. palmitoyl-CoA), while _NUMB
%   specifies that the fatty acid should be represented in number style
%   (e.g. 16:0). A third alternative (not shown here), is _ATE, where the
%   fatty acid is presented as e.g. palmitate.
%
%   If 'reduced' is set to true, a reduced set of reactions are given,
%   filtered for duplicate position-aspecific fatty acids. Example:
%   Reaction 1: chain 1 = 16:0, chain 2: 16:1, chain 3: 18:0
%   Reaction 2: chain 2 = 16:1, chain 2: 18:0, chain 3: 16:0
%   Reaction 3: chain 2 = 16:0, chain 2: 18:0, chain 3: 16:0
%   Both reaction 1 and 2 involve 1x 16:0, 1x 16:1 and 1x 18:0. Reaction 2
%   will therefore be filtered out. This reduces the number of possible
%   combinations, but retains the overall profile.
%
%   Compartment and fatty acid names and abbreviations might have to be
%   adjusted for your model, this is done in the code below, under fillEqns
%   and fillComps.
%
%   Eduard Kerkhoven, 2019-03-08


[newEqns, newNames, chainTrack] = fillEqns(templEqn, templName, 1, chain1);

if nargin > 5 && ~isempty(chain2{1})
    [newEqns, newNames, chainTrack] = fillEqns(newEqns, newNames, 2, chain2, chainTrack);
end
if nargin > 6 && ~isempty(chain3{1})
    [newEqns, newNames, chainTrack] = fillEqns(newEqns, newNames, 3, chain3, chainTrack);
end
if nargin > 7 && ~isempty(chain4{1})
    [newEqns, newNames, chainTrack] = fillEqns(newEqns, newNames, 4, chain4, chainTrack);
end

if reduced == true
    [~, keep, ~] = unique(chainTrack,'rows','stable');
    newEqns     = newEqns(keep);
    newNames    = newNames(keep);
end

if ~isempty(comps)
    [newEqns, newNames] = fillComps(newEqns, newNames, comps);
end
end

function [newEqns, newNames, chainTrack] = fillEqns(templEqns, templNames, position, chains, chainTrack)

chainDB.numb     = {'16:0', '18:0', '18:1', '18:2', '18:3'};
chainDB.ate      = {'palmitate', 'stearate', 'oleate', 'linoleate', ...
    'linolenate'};
chainDB.coa      = {'palmitoyl-CoA', 'stearoyl-CoA', 'oleoyl-CoA', ...
    'linoleoyl-CoA', 'linolenoyl-CoA'};

if nargin<5
    chainTrack = zeros(length(chains),length(chainDB.numb));
else
    chainTrack = repmat(chainTrack, length(chains),1);
end

POS = ['CHAIN' num2str(position)];

chainNumbs      = chains;
[~, chainLoc]   = ismember(chains,chainDB.numb);
chainAtes       = chainDB.ate(chainLoc);
chainCoas       = chainDB.coa(chainLoc);

if ~iscell(templEqns)
    newEqns         = cell(length(chains),1);
    newEqns(:)      = {templEqns};
    newNames        = cell(length(chains),1);
    newNames(:)     = {templNames};
else
    newEqns         = repmat(templEqns,length(chains),1);
    newNames        = repmat(templNames,length(chains),1);
end

idx = 1:1:length(newEqns);
idx = reshape(idx,[length(idx)/length(chains),length(chains)]);

for j = 1:length(chains)
    newEqns(idx(:,j)) = regexprep(newEqns(idx(:,j)),[POS '_NUMB'],chainNumbs{j});
    newEqns(idx(:,j)) = regexprep(newEqns(idx(:,j)),[POS '_ATE'],chainAtes{j});
    newEqns(idx(:,j)) = regexprep(newEqns(idx(:,j)),[POS '_COA'],chainCoas{j});
    newNames(idx(:,j))= regexprep(newNames(idx(:,j)),[POS '_COA'],chainCoas{j});
    newNames(idx(:,j))= regexprep(newNames(idx(:,j)),[POS '_ATE'],chainCoas{j});
    newNames(idx(:,j))= regexprep(newNames(idx(:,j)),POS,chainNumbs{j});
    chainTrack(idx(:,j),chainLoc(j)) = chainTrack(idx(:,j),chainLoc(j)) + 1;
end
end

function [newEqns, newNames] = fillComps(templEqns, templNames, comps)

compsDB.id      = {'c','ce','e','er','erm','g','gm','lp','m','mm','n',...
    'p','v','vm'};
compsDB.name      = {'cytoplasm','cell envelope','extracellular',...
    'ER reticulum','ER membrane','Golgi',...
    'Golgi membrane','lipid particle','mitochondrion',...
    'mitochondrial membrane','nucleus','peroxisome','vacuole',...
    'vacuolar membrane'};

if ~iscell(templEqns)
    templEqns   = {templEqns};
end

[~, compLoc]    = ismember(comps,compsDB.id);
compNames       = compsDB.name(compLoc);

if ~iscell(comps)
    comps       = {comps};
    newEqns     = cell(length(templEqns),1);
    newEqns(:)  = templEqns;
    newNames    = cell(length(templNames),1);
    newNames(:) = templNames;    
else
    newEqns     = repmat(templEqns,length(comps),1);
    newNames    = repmat(templNames,length(comps),1);    
end

idx = 1:1:length(newEqns);
idx = reshape(idx,[length(idx)/length(comps),length(comps)]);

for i=1:length(comps)
    newEqns(idx(:,i))  = regexprep(newEqns(idx(:,i)),'\[COMP\]',['\[' comps{i} '\]']);
    newNames(idx(:,i)) = regexprep(newNames(idx(:,i)),'COMP',['' compNames{i} '']);
end
end
