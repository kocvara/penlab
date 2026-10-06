% PenLab "installer" - to be called from PenLab main directory
%  - adds paths to Matlab Path list
%  - compiles mex files if they are missing for this platform (see build_mex)

% This file is a part of PENLAB package distributed under GPLv3 license
% Copyright (c) 2013 by  J. Fiala, M. Kocvara, M. Stingl
% Last Modified: 27 Nov 2013

dirs = {'doc_generator', 'examples', 'interfaces/BMI', ...
        'interfaces/NLP_AMPL', 'interfaces/PMI', 'interfaces/LMI', ...
        'utilities', 'source'};
currentdir = pwd;

disp('This scripts adds the following PenLab directories to Matlab path:');
len=length(dirs);
for i=1:len
  disp([' * ', fullfile(currentdir, dirs{i})]);
end
disp('To make the change permanent, call savepath.');

for i=1:len
  addpath(fullfile(currentdir, dirs{i}));
end

if (exist('mextrcolumn','file')~=3)
  disp('Mex files are not compiled for this platform, compiling...');
  build_mex;
end
