



adress_folder='data/tst/'%/5m_.001_1_cov10/';

all_files = dir(adress_folder);
all_dir = all_files([all_files(:).isdir]);
num_dir = numel(all_dir)-2;


%fragment_file='data/5m_.001_1_cov10/0/frag0_1.txt';
fragment_file_list=[];
for i=1:num_dir
    all_files_i = dir(strcat(adress_folder,num2str(i)));
    file_num=length({all_files_i.name})-2;
    frags_num=file_num/3; % one conneceted_i_j.txt, fragi_j_sd.txt, fragi_j.txt
for j=1:frags_num
    
     fragment_file=strcat(adress_folder,num2str(i),'/frag',num2str(i),'_',num2str(j),'.txt');
fragment_file_list=[fragment_file_list;fragment_file];

end
end


 for i=1:size(fragment_file_list,1)
     fragment_file=fragment_file_list(i,:)
% %H_final=hap10(fragment_file);
     
 end
    