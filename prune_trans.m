function trans = prune_trans(trans,present_states)
i=1;
while i<=length(trans)
    if present_states(trans(1,i))==0 || present_states(trans(3,i))==0
        trans(:,i)=[];
        continue
    end
    i=i+1;
end