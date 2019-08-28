  penm=ex3_define_by_anonymous();    % to define structure
  prob=penlab(penm);    % to convert the structure and initialize the problem
  prob.solve();         % to start the solver
  prob.x   