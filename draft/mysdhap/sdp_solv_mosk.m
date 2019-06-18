

function X=sdp_solv_mosk(W)

%indices
n=size(W,1);
variables_number=n*(n+1)/2;




%%%%%%% objective function %%%%%%

K_tr=[]; L_tr=[]; % uper triangular (with diag) index of a n*n matrix, fristly, first column, then, 2nd column
for i=1:n
    K_tr=[K_tr,i:n];
    L_tr=[L_tr,i*ones(1,n-i+1)];
end
J=ones(1,variables_number);

W_vec=zeros(1,variables_number);
for i=1:variables_number
    W_vec(i)= W(K_tr(i),L_tr(i));
end

prob.bardim     = [n];
prob.barc.subj  = J;
prob.barc.subk  = K_tr;
prob.barc.subl  = L_tr;
prob.barc.val   = W_vec;



%%%%%%% constraint function %%%%%%


% first constraint X_ij>-.5
variables_number_nondiag=n*(n-1)/2;
I=1:variables_number_nondiag; % each constrains has only one nonzero value in its matrix
J=ones(1,variables_number_nondiag); % only one sdp matrix
K=[]; L=[]; % uper triangular (non-diag) index of a n*n matrix, fristly, first column, then, 2nd column
for i=1:n
    K=[K,i+1:n];
    L=[L,i*ones(1,n-i)];
end
val=.5*ones(1,variables_number_nondiag);

blc = zeros(1,variables_number_nondiag)-0.5;
buc = zeros(1,variables_number_nondiag)+inf;


% % second constraint X_ii=1
I=[I, variables_number_nondiag+1:variables_number_nondiag+n];
J=[J, ones(1,n)]; % there is only one sdp matrix
K=[K, 1:n];
L=[L, 1:n];
val=[val, ones(1,n)];

blc = [blc, ones(1,n)]; % equality contraints, so upper and lower are the same
buc = [buc, ones(1,n)];



%%% constraints

prob.bara.subi = I;
prob.bara.subj = J;
prob.bara.subk = K;
prob.bara.subl = L;
prob.bara.val  = val;
prob.blc = blc;
prob.buc = buc;


prob.a=sparse([], [], [], variables_number_nondiag+n, 0);
%[r, res] = mosekopt('symbcon');
%[r,res] = mosekopt('minimize info',prob);
[r, res] = mosekopt('minimize echo(0)', prob);
X = zeros(n);
K_tr=[]; L_tr=[]; % uper triangular (with diag) index of a n*n matrix, fristly, first column, then, 2nd column
for i=1:n
    K_tr=[K_tr,i:n];
    L_tr=[L_tr,i*ones(1,n-i+1)];
end
x_vec=res.sol.itr.barx;
for i=1:variables_number
    X(K_tr(i),L_tr(i))=x_vec(i);
end

X = X + tril(X,-1)';






end
