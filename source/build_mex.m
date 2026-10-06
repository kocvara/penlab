% Compile PenLab mex files for the current platform (needs a C compiler,
% see 'mex -setup C'). The binaries are placed next to the sources.
function build_mex()

  srcdir = fileparts(mfilename('fullpath'));
  files = {'mexsumsparse.c', 'mextrcolumn.c', 'mextrdsdsmat.c'};
  for i=1:length(files)
    fprintf('Compiling %s\n', files{i});
    mex('-O', '-largeArrayDims', '-outdir', srcdir, fullfile(srcdir, files{i}));
  end

end
