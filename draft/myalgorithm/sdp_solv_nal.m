
function X_out=sdp_solv_nal(W)

%W=[1,2,3;4,5,6;7,8,9];
%W=W+W';
N=size(W,1);



blk{1,1} = 's';
blk{1,2} = N;
C=cell(1);
L=cell(1);
At=cell(1);

C{1}=sparse(W);
L{1}=-.5*ones(N);

A_all=[];

for i=1:N
    A_i=zeros(N,N); 
    A_i(i,i)=1;
    A_i_svec=svec(blk(1,:),A_i); 
    A_all=[A_all, A_i_svec];
end


%A_2=zeros(N); A_2(2,2)=1; At_2=vec(A_2);
%At{1}=[At_1';At_2'];
% At{1}=[At_1';At_2'];
% b1=[1; 1]; %ones(N,1);

At{1}= A_all; %svec(blk(1,:),rand(3));

b=ones(1,N);

OPTIONS.printlevel=0;

X_init=cell(1);X_init{1}=ones(N);%X_0;

% X_0=-.5*ones(N);X_0(1:N+1:N^2)=ones(1,N);
% X_init=cell(1);X_init{1}=X_0;

%[obj,X,s,y,Z1,Z2,y2,v,info,runhist]=sdpnalplus(blk,At,C,b,L,[],[],[],[],OPTIONS);
[obj,X,s,y,Z1,Z2,y2,v,info,runhist]=sdpnalplus(blk,At,C,b,L,[],[],[],[],OPTIONS,X_init);
 X_out=X{1};        
end     
         