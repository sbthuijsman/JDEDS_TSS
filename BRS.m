function Xplus = BRS(X,Sigma,trans,Xinit)
Xplus_1 = (Xinit.*X)~=0; Xplus_2 = zeros(1,length(X));
Sigma=Sigma~=0;
while(true)
    Xplus_2 = (Xplus_1+Xplus_2)~=0;
    for i=1:length(trans)
        t = trans(:,i);
        if X(t(1))==1 && Sigma(t(2))==1 && X(t(3))==1
        set_n_trans_trav(get_n_trans_trav+1); %number of traversed edges
            if Xplus_2(t(1))~=1 && Xplus_1(t(3))==1
                Xplus_2(t(1)) = 1;
                Xplus_2 = (Xplus_1+Xplus_2)~=0; %union with 1 and 0 elements
            end
        end
    end

    if isequal(Xplus_2,Xplus_1)
        Xplus = Xplus_2;
        break;
    else
        Xplus_1=Xplus_2;
    end
end

end