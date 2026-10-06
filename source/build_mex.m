% Compile PenLab mex files for the current platform (needs a C compiler,
% see 'mex -setup C'). The binaries are placed next to the sources.
%
% build_mex()        ... the three mex files used by the solver
% build_mex(asldir)  ... in addition, the AMPL interface 'amplf', linked
%                        against the AMPL Solver Library compiled in asldir
%                        (with amplsolver.a and arith.h), e.g., from
%                        https://netlib.org/ampl/solvers.tgz:
%                          ./configurehere
%                          make CFLAGS="-O2 -fPIC -DASL_NO_FPINITMT"
%                        (on macOS add -mmacosx-version-min=12.0)
function build_mex(asldir)

  srcdir = fileparts(mfilename('fullpath'));
  files = {'mexsumsparse.c', 'mextrcolumn.c', 'mextrdsdsmat.c'};
  for i=1:length(files)
    fprintf('Compiling %s\n', files{i});
    mex('-O', '-largeArrayDims', '-outdir', srcdir, fullfile(srcdir, files{i}));
  end

  if (nargin>=1)
    fprintf('Compiling amplf against ASL in %s\n', asldir);
    mexdir = fullfile(srcdir, '..', 'utilities', 'mex');
    mex('-O', '-largeArrayDims', '-DMEMORY_MATLAB', '-DASL_NO_FPINITMT', ...
        ['-I' asldir], fullfile(mexdir, 'my_amplfunc4.c'), ...
        fullfile(asldir, 'amplsolver.a'), '-output', 'amplf', '-outdir', srcdir);
  end

end
