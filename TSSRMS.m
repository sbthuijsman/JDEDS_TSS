function [Yprime, Gprime, skipped, deltaCalc] = TSSRMS(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, xdelta, skipCalc, doGroup)
    deltaCalc=[];
    skipped=false;
    
    if skipCalc==false && doGroup==false
        %we are performing ITSS
        
        if (G(xdelta)-Y(xdelta))==1
            arg1 = ((Y+X0)>0);
            arg2 = ((Y+Xm)>0);
            arg2(xdelta)=0;
            [Yprime,Gprime] = DLSS(G,Sigma_c,Sigma_u,transX,arg1,arg2);
        elseif Y(xdelta)==1
            arg = Xm;
            arg(xdelta)=0;
            [Yprime,Gprime] = DLSS(G,Sigma_c,Sigma_u,transX,X0,arg);
        else
            Yprime = Y; Gprime=G;
        end

    elseif skipCalc==true && doGroup==false
        %skip FRS, BRS and SS calls
        
        if (G(xdelta)-Y(xdelta))==1
            skipped = true; Yprime = Y; Gprime=G;

        elseif Y(xdelta)==1
            skipped = true; Yprime = Y; Gprime=G;

        else
            Yprime = Y; Gprime=G;
        end

    elseif skipCalc==false && doGroup==true
        %perform grouped FRS, BRS, or SS call
        
        Gr1=zeros(size(xdelta));
        Gr2=zeros(size(xdelta));
        for xd=1:length(xdelta)
            if xdelta(xd)==1
                if (G(xd)-Y(xd))==1
                    Gr1(xd)=1;
                elseif Y(xd)==1
                    Gr2(xd)=1;
                else
                    0;
                end
            end
        end
        
        if sum(Gr1)>0
            arg1 = ((Y+X0)>0);
            arg2 = ((Y+Xm+Gr1)>0);
            arg2 = ((arg2-Gr1)>0);
            [Yprime,Gprime] = DLSS(G,Sigma_c,Sigma_u,transX,arg1,arg2);
            deltaCalc=Gr1;
        elseif sum(Gr2)>0
            arg = ((Xm-Gr1)>0);
            [Yprime,Gprime] = DLSS(G,Sigma_c,Sigma_u,transX,X0,arg);
            deltaCalc=Gr2;
        else
            skipped =true; Yprime = Y; Gprime=G;
        end
        
    end
end