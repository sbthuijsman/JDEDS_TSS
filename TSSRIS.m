function [Yprime, Gprime, skipped, deltaCalc] = TSSRIS(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, xdelta, skipCalc, doGroup)
    deltaCalc=[];
    skipped=false;
    
    if skipCalc==false && doGroup==false
        %we are performing ITSS

        if Y(xdelta)==1
            arg = X0;
            arg(xdelta)=0;
            Yprime = FRS(G,(Sigma_c+Sigma_u)>0,transX,arg);
            Gprime = G;
        else
            Yprime = Y; Gprime=G;
        end

    elseif skipCalc==true && doGroup==false
        %skip FRS, BRS and SS calls

        if Y(xdelta)==1
            skipped = true; Yprime = Y; Gprime=G;

        else
            Yprime = Y; Gprime=G;
        end
        
    elseif skipCalc==false && doGroup==true
        %perform grouped FRS, BRS, or SS call
        
        Gr1=zeros(size(xdelta));
        for xd=1:length(xdelta)
            if xdelta(xd)==1
                if Y(xd)==1
                    Gr1(xd)=1;
                else
                    0;
                end
            end
        end
        
        if sum(Gr1)>0
            arg = ((X0-Gr1)>0);
            Yprime = FRS(Y,(Sigma_c+Sigma_u)>0,transX,arg);
            Gprime = G;
            deltaCalc = Gr1;
        else
            skipped = true; Yprime = Y; Gprime=G;
        end
        
    end
end