function [X, Sigma_c, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, transminindices, transplus, X0plus, X0min, Xmplus, Xmmin] = ITSS(X, Sigma_c, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, transminindices, transplus, X0plus, X0min, Xmplus, Xmmin, skipCalc, doGroup)

if doGroup==true
    error();
end

for xd=1:length(X0plus)
    if X0plus(xd)==1
        [Y, G, skipped, ~] = TSSAIS(X, Sigma_c, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, xd, skipCalc, doGroup);
        if skipped==false
            X0(xd)=1;
            X0plus(xd)=0;
        end
    end
end
for xd=1:length(X0min)
    if X0min(xd)==1
        [Y, G, skipped, ~] = TSSRIS(X, Sigma_c, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, xd, skipCalc, doGroup);
        if skipped==false
            X0(xd)=0;
            X0min(xd)=0;
        end
    end
end

for xd=1:length(Xmplus)
    if Xmplus(xd)==1
        [Y, G, skipped, ~] = TSSAMS(X, Sigma_c, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, xd, skipCalc, doGroup);
        if skipped==false
            Xm(xd)=1;
            Xmplus(xd)=0;
        end
    end
end
for xd=1:length(Xmmin)
    if Xmmin(xd)==1
        [Y, G, skipped, ~] = TSSRMS(X, Sigma_c, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, xd, skipCalc, doGroup);
        if skipped==false
            Xm(xd)=0;
            Xmmin(xd)=0;
        end
    end
end

transplus_skipped = [];
for ti = 1:size(transplus,2)
    transdelta = transplus(:,ti);
    [Y, G, skipped, ~, ~] = TSSAT(X, Sigma_c, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, transdelta, skipCalc, doGroup);
    if skipped==false
        transX = [transX transdelta]; %comment later, keep up to date for checks
        xor = transplus(1,ti);
        sig = transplus(2,ti);
        xtar = transplus(3,ti);
        adj{xor}(end+1)=xtar;
        adjrev{xtar}(end+1)=xor;
        if Sigma_u(sig)
            adj_urev{xtar}(end+1)=xor;
        end
    else
        transplus_skipped = [transplus_skipped transdelta];
    end
end
transplus = transplus_skipped;

transmin_skipped = [];
for ti = transminindices
    [Y, G, skipped, ~] = TSSRT(X, Sigma_c, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, ti, skipCalc, doGroup);
    if skipped==false
        %remove only one each time
        xor = transX(1,ti);
        sig = transX(2,ti);
        xtar = transX(3,ti);
        adj{xor}(find(adj{xor}==xtar,1))=[]; %removes first element
        adjrev{xtar}(find(adjrev{xtar}==xor,1))=[];
        if Sigma_u(sig)
            adj_urev{xtar}(find(adj_urev{xtar}==xor,1))=[];
        end
    else
        transmin_skipped = [transmin_skipped ti];
    end
end
transminindices = transmin_skipped;

end