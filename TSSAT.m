function [Yprime, Gprime, skipped, deltaCalc, deltaCalcInd] = TSSAT(X, ~, Sigma_u, ~, adj, adjrev, adj_urev, X0, Xm, Y, G, transdelta, skipCalc, doGroup)
    deltaCalc=[];
    deltaCalcInd=-1;
    skipped=false;
        
    if skipCalc==false && doGroup==false
        %we are performing ITSS

        xor=transdelta(1);
        sigma=transdelta(2);
        xtar=transdelta(3);
        
        if Y(xor)==1 && (G(xtar)-Y(xtar))==1
            adj{xor}(end+1)=xtar;

            Yprime = FRS(G,adj,Y);
            Gprime = G;        

        elseif G(xor)==1 && (X(xtar)-G(xtar))==1 && Sigma_u(sigma)==1

            if ~isempty(adj{xor})
                for xtar = adj{xor}
                    adjrev{xtar} = adjrev{xtar}(adjrev{xtar}~=xor); %remove all xor from adjrev{xtar}
                    adj_urev{xtar} = adj_urev{xtar}(adj_urev{xtar}~=xor); %remove all xor from adjrev{xtar}
                end
            end
            adj{xor}=[];

                
            arg1 = X0;
            arg1(xor)=0;
            arg2 = Xm;
            arg2(xor)=0;
            [Yprime,Gprime] = DLSS(G, adj, adjrev, adj_urev, arg1, arg2);
            
        elseif (X(xor)-G(xor))==1 && xor~=xtar
            adj{xor}(end+1)=xtar;
            adjrev{xtar}(end+1)=xor;
            if Sigma_u(sigma)
                adj_urev{xtar}(end+1)=xor;
            end

            arg1 = ((Y+X0)>0);
            arg2 = ((Y+Xm)>0);
            [Yprime,Gprime] = DLSS(X, adj, adjrev, adj_urev, arg1, arg2);
            
        else
            Yprime = Y; Gprime=G;
        end

    elseif skipCalc==true && doGroup==false
        %skip FRS, BRS and SS calls
        
        xor=transdelta(1);
        sigma=transdelta(2);
        xtar=transdelta(3);
                
        if Y(xor)==1 && (G(xtar)-Y(xtar))==1
            skipped = true; Yprime = Y; Gprime=G;
            
        elseif G(xor)==1 && (X(xtar)-G(xtar))==1 && Sigma_u(sigma)==1
            skipped = true; Yprime = Y; Gprime=G;
            
        elseif (X(xor)-G(xor))==1 && xor~=xtar
            skipped = true; Yprime = Y; Gprime=G;
            
        else
            Yprime = Y; Gprime=G;
        end
        
    elseif skipCalc==false && doGroup==true
        %perform grouped FRS, BRS, or SS call
                
        Gr1=[];
        Gr2=[];
        Gr3=[];
        Gr1ind=[];
        Gr2ind=[];
        Gr3ind=[];
        for trans=1:size(transdelta,2)
            td=(transdelta(:,trans));
            xor=td(1);
            sigma=td(2);
            xtar=td(3);
            
            if Y(xor)==1 && (G(xtar)-Y(xtar))==1
                Gr1=[Gr1 td];
                Gr1ind=[Gr1ind trans];
            elseif G(xor)==1 && (X(xtar)-G(xtar))==1 && Sigma_u(sigma)==1
                Gr2=[Gr2 td];
                Gr2ind=[Gr2ind trans];
            elseif (X(xor)-G(xor))==1 && xor~=xtar
                Gr3=[Gr3 td];
                Gr3ind=[Gr3ind trans];
            else
                0;
            end
        end
        
        if sum(Gr1)>0
            for i=1:size(Gr1,2)
                adj{Gr1(1,i)}(end+1)=Gr1(3,i);
            end
            Yprime = FRS(G,adj,Y);
            Gprime = G;
            deltaCalc=Gr1;
            deltaCalcInd=Gr1ind;
            
        elseif sum(Gr2)>0
            arg1 = X0;
            arg2 = Xm;
            
            %for all trans in group
            for ti = 1:size(Gr2,2)
                xor = Gr2(1,ti);
                
                arg1(xor)=0;
                arg2(xor)=0;
                if ~isempty(adj{xor})
                    for xtar = adj{xor}
                        adjrev{xtar} = adjrev{xtar}(adjrev{xtar}~=xor); %remove all xor from adjrev{xtar}
                        adj_urev{xtar} = adj_urev{xtar}(adj_urev{xtar}~=xor); %remove all xor from adj_urev{xtar}
                    end
                end
                adj{xor}=[];

            end
            
            [Yprime,Gprime] = DLSS(G, adj, adjrev, adj_urev, arg1, arg2);
            deltaCalc=Gr2;
            deltaCalcInd=Gr2ind;

        elseif sum(Gr3)>0

            for i=1:size(Gr3,2)
                xor=Gr3(1,i);
                sig=Gr3(2,i);
                xtar=Gr3(3,i); 

                adj{xor}(end+1)=xtar;
                adjrev{xtar}(end+1)=xor;
                if Sigma_u(sig)
                    adj_urev{xtar}(end+1)=xor;
                end
            end

            arg1 = ((Y+X0)>0);
            arg2 = ((Y+Xm)>0);
            [Yprime,Gprime] = DLSS(X, adj, adjrev, adj_urev, arg1, arg2);
            deltaCalc=Gr3;
            deltaCalcInd=Gr3ind;
            
        else
            skipped = true; Yprime = Y; Gprime=G;
        end
        
    end
end