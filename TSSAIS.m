function [Yprime, Gprime, skipped, deltaCalc] = TSSAIS(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, xdelta, skipCalc, doGroup)
    deltaCalc=[];
    skipped=false;
    
    if skipCalc==false && doGroup==false
        %we are performing ITSS
        if length(xdelta)~=1
            error = "unexpected delta length"
            error()
        end

        GsetminY=((G-Y)>0);

        if GsetminY(xdelta)==1
            arg = Y;
            arg(xdelta)=1;
            Yprime = FRS(G,(Sigma_c+Sigma_u)>0,transX,arg);
            Gprime = G;
        else
            Yprime = Y; Gprime=G;
        end

    elseif skipCalc==true && doGroup==false
        %skip FRS, BRS and SS calls
        if length(xdelta)~=1
            error = "unexpected delta length"
            error()
        end
        
        GsetminY=((G-Y)>0);

        if GsetminY(xdelta)==1
            skipped = true; Yprime = Y; Gprime=G;
            
        else
             Yprime = Y; Gprime=G;
        end
        
    elseif skipCalc==false && doGroup==true
        %perform grouped FRS, BRS, or SS call
        
        GsetminY=((G-Y)>0);
        Gr1=zeros(size(xdelta));
        for xd=1:length(xdelta)
            if xdelta(xd)==1
                if GsetminY(xd)==1
                    Gr1(xd)=1;
                else
                    0;
                end
            end
        end
        
        if sum(Gr1)>0
            arg = ((Y+Gr1)>0);
            Yprime = FRS(G,(Sigma_c+Sigma_u)>0,transX,arg);
            Gprime = G;
            deltaCalc = Gr1;
        else
            skipped = true; Yprime = Y; Gprime=G;
        end
        
    end
end