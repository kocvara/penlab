clear bmidata;
bmidata.name='Example 3 from BMI manual';
bmidata.Nx=3;
bmidata.Na=1;
bmidata.Ng=4;
bmidata.B=[-1 0 0;1 0 0;0 -1 0;0 1 0]; % ...... matrix defining the linear inequality constraints Bx<=d
%            dimension Ng x Nx
bmidata.d=[.5; 2; 3; 7];% ...... rhs for linear constraints
bmidata.c=[0; 0; 1];% ...... dim (Nx,1), coefficients of the linear objective function
bmidata.NaDims=[3];% . vector of sizes of matrix constraints (diagonal blocks)
bmidata.A{1,1} = sparse([-10 -.5 -2;-.5 4.5 0; -2 0 0]);
bmidata.A{1,2} = -sparse([9 .5 0; .5 0 -3; 0 -3 -1]);
bmidata.A{1,3} = -sparse([-1.8 -.1 -.4; -.1 1.2 -1; -.4 -1 0]);
bmidata.A{1,4} = speye(3,3);
%...... cell array (matrix) of A{k,l} for k=1,...,Na matrix constraint
%            for l=1 ~ absolute term, l=2..Nx+1 coeficient matrices
%            (some of them might be empty)

bmidata.K{1,1,2} = -sparse([0 0 2;0 -5.5 3; 2 3 0]);
bmidata.indK{1} = [1 2; 0 0]; %array of ij indices for existing K_ij matrices
% one row for each variable (but the last one)
% first column gives number of nonzero indices j, following columns are
% these indices

bmidata.Adep{1}=[1 2 3]; %... rather create here vvv !!!
