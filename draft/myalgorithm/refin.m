


function H_final=refin(R,H)

[K, L]=size(H); % k rows , k ploidy level

H_best_upnow=H;
KK=4;
mec_all=zeros(1,KK);

mec_all(1)=mec_calculator(R,H_best_upnow);
kk=1;
continu=1;
while (kk<KK && continu)  % I can save all of them and extract the best mec, no it can become worde
    kk=kk+1;
    
    for k=1:K % row  of H
        for l=1:L % col of H     
            H_check=H_best_upnow;
            H_check(k,l) = -H_best_upnow(k,l);
            mec_check_l=mec_calculator_l(R,H_check,l);
            mec_best_l=mec_calculator_l(R,H_best_upnow,l);
            if mec_check_l < mec_best_l
                H_best_upnow=H_check;
            end
        end
    end


mec_all(kk) = mec_calculator(R,H_best_upnow);

if ( (mec_all(kk)==mec_all(kk-1)) || mec_all(kk)==0 )
    continu=0;
end

mec_all_end=mec_all(1:kk);
mec_all_end
H_final=H_best_upnow;
end
