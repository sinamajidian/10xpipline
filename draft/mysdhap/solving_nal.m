% solving
% clearvars
name_mat='data/simulation1a/1a_1.mat';%'data/65a_1_20/20.mat';%'data/simulation1a/1a_1.mat';%'data/old/sim2_b_R1';%'data/65a_1_20/20.mat' %'data/simulation1a/1a_0.mat'; %'data/simulation1a/1_a_1.mat'; %'sim2_b_R0.mat';
name_hap='data/simulation1a/nal_nim1_.hap'; %'data/65a_1_20/sdhap_matlab_nal.hap';
load(name_mat) 

%R1=R;
%R=R1(1:200,1:200)

k=3;
%R=R(1:40,:);
N=size(R,1);
l=size(R,2);
W=zeros(N,N);
diag_const=0;
W(1,1)=diag_const;
for i=2:N
    line_i=R(i,:);
    W(i,i)=diag_const;
    for j=1:(i-1)
        line_j=R(j,:);
        SNP_shared=sum( (line_i~=0) & (line_j~=0));
        line_j(line_j==0)=33; % 33 is arbitary in order to check only in nonzero values
        allele_shared=sum( line_i==line_j);
        if SNP_shared>0
            W(i,j)=(2*allele_shared-SNP_shared)/SNP_shared;
        else
             W(i,j)=0;
        end
        
        W(j,i)=W(i,j);
    end
end

size(W)
%  
% 
% W=[1:3;4:6;7:9];
% N=3;
% V=zeros(N,2)+1;
% i=1;
% 
% 
% for i=1:N
% V_i=sum(W(i,:)'.*V);
% V(i,:)=V_i/norm(V_i);
% end






X=sdp_solv_nal(-W);


%X=sparse(X);
[Q, sig]=eig(X);
[val_eig, idx]=sort(diag(sig), 'descend'); % ascend[1,2,3]  descend [2,3,1]
three_ind_largest=idx(1:3); % sometimes  the third is zero
val_eig(1:7)'

V=Q(:,three_ind_largest)*sqrt(sig(three_ind_largest,three_ind_largest));


V_arch=V;



object_all=[]
indx_all=[];

for ii=1:100*N
    Z=normrnd(0,1,[k,k]); %Z=normalize(Z);
    VZ=V*Z;
    [val, index]=max(VZ');
    
    X_estimated=ones(N,N);
    for i=1:N
        for j=1:N
            if index(i)~=index(j)
                X_estimated(i,j)=-1;
            end
        end
    end
    object_all=[object_all; trace(W*X_estimated)  ];
    indx_all=[indx_all;index];
end
[vall,i_best]=max(object_all);
index_best=indx_all(i_best,:);

size(index_best)

R=full(R);
H=zeros(3,l);
for i_k=1:k
   R(index_best==i_k,:);
   H(i_k,:)=sum(R(index_best==i_k,:))>0;
end

H1=H'+1;



%51:70
%H2=[[23:26]', H'+1]

%save('sim2_b_R.mat','-v7.3')


%%% greedy refinement
element_num=size(H,1)*size(H,2);

H_new=2*H-1;
mec_all=[mec_calculator(R,H_new)];
for kk=1:8 % I can save all of them and extract the best mec, no it can become worde
for ii=1:element_num
    
    H_check=H_new;
    H_check(ii)=-H_new(ii);
    if mec_calculator(R,H_check)< mec_calculator(R,H_new)
        H_new=H_check;
    end   
end

mec_all=[mec_all, mec_calculator(R,H_new)];
end
mec_all
H_final=H_new;

(H_final'+1)/2+1;
% H_sd=2*(H_sd-1)-1;
% sum(abs(H_sd(1,:)-H_final(3,:)))
indces_block=hap_index'-1;  % The output file will be like sdhap. index starts from zero
fileID_hap = fopen(name_hap,'w');
fprintf(fileID_hap,'Block 1\t Length of haplotype block %d\t Number of read %d\t Total MEC s \n',length(indces_block),N);
H_with_ind=[indces_block, (H_final'+1)/2+1]
fprintf(fileID_hap,'%d\t%d\t%d\t%d\n',H_with_ind');
% 



% 
% 
% %python2 $hapcompare ../../genome_generated/simulated_varianthaplos.txt sdhap_matlab.hap  -t -v 

