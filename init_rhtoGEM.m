% Define path to folders
root  = pwd();
scripts = [root '/ComplementaryScripts']; 
data    = [root '/ComplementaryData'];

% Check whether RAVEN is installed
if ~exist('ravenCobraWrapper.m')
    warning(['RAVEN does not seem to be installed. For model curation, '...
        'you are advised to install RAVEN Toolbox, which is available '...
        'from http://github.com/SysBioChalmers/RAVEN. For model '...
        'simulations, rhto-GEM is compatible with e.g. COBRA Toolbox.'])
end

% Clear
