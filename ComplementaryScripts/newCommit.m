function newCommit(model)
% newCommit
%   Prepares a new commit of the Rhto model in a development branch, for
%   direct submission to GitHub. In contrast to new releases, new commits
%   in development branches do not contain .mat and .xlsx files. This
%   function should be run from the ComplementaryScripts directory.
%
%   model   RAVEN model structure to be used in commit. If no model
%           structure is provided, the stored XML file will be imported
%           prior to exporting other files
%
%   Usage: newCommit(model)
%
% Eduard Kerkhoven, 2018-05-15

%Check if in master:
currentBranch = git('rev-parse --abbrev-ref HEAD');
if strcmp(currentBranch,'master')
    error('ERROR: in master branch. For new releases, use newRelease.m')
end

if ~exist('model')
    %Load model:
    model = importModel('../ModelFiles/xml/rhto.xml');
end

%Sort model
[~,i]=sort(model.rxns);
model=permuteModel(model,i,'rxns');
[~,i]=sort(model.mets);
model=permuteModel(model,i,'mets');
%Save model
exportForGit(model,'rhto','../',{'txt', 'xml', 'yml'});
end