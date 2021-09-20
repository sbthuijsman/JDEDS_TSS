function trans = quick_prune_exact_trans(trans,transremove)
tsize = size(transremove,2);
tindices = zeros(1,tsize);
    for ti=1:tsize
        tindices(ti)=find(sum(trans==transremove(:,ti))==3);
    end
trans(:,tindices)=[];
end