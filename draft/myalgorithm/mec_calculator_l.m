
function mec=mec_calculator_l(R_matrix,H_candidate,l)

% R {+1,-1,0}   h {1,-1}
%R_matrix_l=R_matrix(:,l);
[row_list,~,~]=find(R_matrix(:,l)); % row list that are non zero in l-th column

mec=0;
for i=1:length(row_list)
    R_i=R_matrix(row_list(i),:);
    [~, R_i_ind,R_i_val]=find(R_i); %[I,J,value]=find() % this is a vector so I=ones
    H_i=H_candidate(:,R_i_ind);
    diff=repmat(R_i_val,3,1)-H_i;
    sm=sum(abs(diff),2)/2;
    mec=mec+min(sm);
end

end



