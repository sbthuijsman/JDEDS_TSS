function [X, Sigma_c, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, transminindices, transplus, X0plus, X0min, Xmplus, Xmmin] = GTSS(X, Sigma_c, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, transminindices, transplus, X0plus, X0min, Xmplus, Xmmin, skipCalc, doGroup)

if doGroup==false
    error();
end

if sum(X0plus)>0
    [Y, G, skipped, deltaCalc] = TSSAIS(X, Sigma_c, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, X0plus, skipCalc, doGroup);
    if skipped==false
        X0 = (X0+deltaCalc)>0;
        X0plus = (X0plus-deltaCalc)>0;
        return
    end
end

if sum(X0min)>0
    [Y, G, skipped, deltaCalc] = TSSRIS(X, Sigma_c, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, X0min, skipCalc, doGroup);
    if skipped==false
        X0 = (X0-deltaCalc)>0;
        X0min = (X0min-deltaCalc)>0;
        return
    end
end

if sum(Xmplus)>0
    [Y, G, skipped, deltaCalc] = TSSAMS(X, Sigma_c, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, Xmplus, skipCalc, doGroup);
    if skipped==false
        Xm = (Xm+deltaCalc)>0;
        Xmplus = (Xmplus-deltaCalc)>0;
        return
    end
end


if sum(Xmmin)>0
    [Y, G, skipped, deltaCalc] = TSSRMS(X, Sigma_c, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, Xmmin, skipCalc, doGroup);
    if skipped==false
        Xm = (Xm-deltaCalc)>0;
        Xmmin = (Xmmin-deltaCalc)>0;
        return
    end
end


if size(transplus,2)>0
    [Y, G, skipped, deltaCalc, deltaCalcInd] = TSSAT(X, Sigma_c, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, transplus, skipCalc, doGroup);
    if skipped==false

        for i=1:size(deltaCalc,2)
            xor=deltaCalc(1,i);
            sig=deltaCalc(2,i);
            xtar=deltaCalc(3,i); 
        
            adj{xor}(end+1)=xtar;
            adjrev{xtar}(end+1)=xor;
            if Sigma_u(sig)
                adj_urev{xtar}(end+1)=xor;
            end
        end

        %remove deltaCalc from transplus, has been processed
        transplus(:,deltaCalcInd)=[];

        
        return
    end
end

if size(transminindices,2)>0
    [Y, G, skipped, deltaCalc] = TSSRT(X, Sigma_c, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, transminindices, skipCalc, doGroup);
    if skipped==false

            for i=deltaCalc
                %remove only one each time
                xor=transX(1,i);
                sig=transX(2,i);
                xtar=transX(3,i); 
                adj{xor}(find(adj{xor}==xtar,1))=[]; %removes first element
                adjrev{xtar}(find(adjrev{xtar}==xor,1))=[];
                if Sigma_u(sig)
                    adj_urev{xtar}(find(adj_urev{xtar}==xor,1))=[];
                end

            end 

        %remove deltaCalc from transmin, has been processed
        transminindices = setdiff(transminindices,deltaCalc);

        return
    end
end

end