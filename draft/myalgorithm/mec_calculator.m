
function mec=mec_calculator(R_matrix,H_candidate)
% R {+1,-1,0}   h {1,-1}
N=size(R_matrix,1);
l=size(H_candidate,1);

mec=0;
hd=zeros(1,3);
for i=1:N
    for j=1:l
    hd(j)=sum(abs(R_matrix(i,:)-H_candidate(j,:))==2);
  
    end
   mec=mec+min(hd); 
end


end



