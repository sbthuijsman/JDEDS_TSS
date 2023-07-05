%% DLSS (Ouedrago synth)
function [Y, G] = DLSS(X, adj, adjrev, adj_urev, X0, Xm)
    Xk1=X;
    while(true)
        Nk = FRS(Xk1,adjrev,Xk1.*Xm); %BRS
        Bk = FRS(Xk1,adj_urev,(X.*(Xk1-Nk))>0); %BRS

        Xk2 = X.*((Xk1-Bk))>0;
        if isequal(Xk2,Xk1)
            G = Xk2;
            break;
        else
            Xk1 = Xk2;
        end
    end

    Y = FRS(G,adj,X0.*G);

end