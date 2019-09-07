% PenLab "installer" - to be called from PenLab main directory
%  - adds paths to Matlab Path list
%  - later -> compile mex

% This file is a part of PENLAB package distributed under GPLv3 license
% Copyright (c) 2013-2019 by  J. Fiala, M. Kocvara, M. Stingl
% Last Modified: 6 Sep 2019

dirs = {'doc_generator', 'examples', 'interfaces/BMI', ...
    'interfaces/NLP_AMPL', 'interfaces/PMI', 'interfaces/LMI', ...
    'utilities', 'source'};
if ismac
    slash = '/';
else
    slash = '\';
end
currentdir = pwd;

lastchar = currentdir(length(currentdir));
if (lastchar~='/' && lastchar~='\')
    currentdir = [currentdir, slash];
end

disp('This scripts adds the following PenLab directories to Matlab path:');
len=length(dirs);
for i=1:len
    disp([' * ', currentdir, dirs{i}]);
end
disp('To make the change permanent, call savepath.');

for i=1:len
    addpath([currentdir, dirs{i}]);
end

