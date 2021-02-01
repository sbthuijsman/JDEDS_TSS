function [X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transmin, transplus, X0plus, X0min, Xmplus, Xmmin] = GTSS(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transmin, transplus, X0plus, X0min, Xmplus, Xmmin, skipCalc, doGroup)

if doGroup==false
    error();
end

if sum(X0plus)>0
    [Y, G, skipped, deltaCalc] = TSSAIS(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, X0plus, skipCalc, doGroup);
    if skipped==false
        X0 = (X0+deltaCalc)>0;
        X0plus = (X0plus-deltaCalc)>0;
        return
    end
end

if sum(X0min)>0
    [Y, G, skipped, deltaCalc] = TSSRIS(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, X0min, skipCalc, doGroup);
    if skipped==false
        X0 = (X0-deltaCalc)>0;
        X0min = (X0min-deltaCalc)>0;
        return
    end
end

if sum(Xmplus)>0
    [Y, G, skipped, deltaCalc] = TSSAMS(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, Xmplus, skipCalc, doGroup);
    if skipped==false
        Xm = (Xm+deltaCalc)>0;
        Xmplus = (Xmplus-deltaCalc)>0;
        return
    end
end


if sum(Xmmin)>0
    [Y, G, skipped, deltaCalc] = TSSRMS(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, Xmmin, skipCalc, doGroup);
    if skipped==false
        Xm = (Xm-deltaCalc)>0;
        Xmmin = (Xmmin-deltaCalc)>0;
        return
    end
end

if size(transplus,2)>0
    [Y, G, skipped, deltaCalc] = TSSAT(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transplus, skipCalc, doGroup);
    if skipped==false
        transX = [transX deltaCalc];
        %remove deltaCalc from transplus, has been processed
        for ti = 1:size(deltaCalc,2)
            xor = deltaCalc(1,ti);
            sigma = deltaCalc(2,ti);
            xtar = deltaCalc(3,ti);

            edgfrmxor = find(transplus(1,:)==xor);

            for j=1:length(edgfrmxor)
                t = transplus(:,edgfrmxor(j));
                if t(2)==sigma && t(3)==xtar
                    transplus(:,edgfrmxor(j))=[];
                    break; %this assumes this exact transition only exists once in transX
                end
            end
        end
        
        return
    end
end

if size(transmin,2)>0
    [Y, G, skipped, deltaCalc] = TSSRT(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transmin, skipCalc, doGroup);
    if skipped==false
        %remove deltaCalc from transX, has been processed
        for ti = 1:size(deltaCalc,2)
            xor = deltaCalc(1,ti);
            sigma = deltaCalc(2,ti);
            xtar = deltaCalc(3,ti);

            edgfrmxor = find(transX(1,:)==xor);

            for j=1:length(edgfrmxor)
                t = transX(:,edgfrmxor(j));
                if t(2)==sigma && t(3)==xtar
                    transX(:,edgfrmxor(j))=[];
                    break; %this assumes this exact transition only exists once in transX
                end
            end
        end

        %remove deltaCalc from transmin, has been processed
        for ti = 1:size(deltaCalc,2)
            xor = deltaCalc(1,ti);
            sigma = deltaCalc(2,ti);
            xtar = deltaCalc(3,ti);

            edgfrmxor = find(transmin(1,:)==xor);

            for j=1:length(edgfrmxor)
                t = transmin(:,edgfrmxor(j));
                if t(2)==sigma && t(3)==xtar
                    transmin(:,edgfrmxor(j))=[];
                    break; %this assumes this exact transition only exists once in transX
                end
            end
        end

        return
    end
end

end