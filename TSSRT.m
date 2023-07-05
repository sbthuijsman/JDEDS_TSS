function [Yprime, Gprime, skipped, deltaCalc] = TSSRT(X, ~, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, transminind, skipCalc, doGroup)
    deltaCalc=[];
    skipped=false;
    
    transdelta = transX(:,transminind);
    
    if skipCalc==false && doGroup==false
        %we are performing ITSS
        
        xor=transdelta(1);
        sigma=transdelta(2);
        xtar=transdelta(3);
        
        if Y(xor)==1 && Y(xtar)==1 && xor~=xtar

            adj{xor}(find(adj{xor}==xtar,1))=[]; %removes first element
            adjrev{xtar}(find(adjrev{xtar}==xor,1))=[];
            if Sigma_u(sigma)
                adj_urev{xtar}(find(adj_urev{xtar}==xor,1))=[];
            end

             [Yprime,Gprime] = DLSS(G, adj, adjrev, adj_urev, X0, Xm);
            
        elseif (G(xor)-Y(xor))==1 && G(xtar)==1 && xor~=xtar
            
            
            adj{xor}(find(adj{xor}==xtar,1))=[]; %removes first element
            adjrev{xtar}(find(adjrev{xtar}==xor,1))=[];
            if Sigma_u(sigma)
                adj_urev{xtar}(find(adj_urev{xtar}==xor,1))=[];
            end

             [Yprime,Gprime] = DLSS(G, adj, adjrev, adj_urev, (Y+X0)>0,(Y+Xm)>0);
            
        elseif (X(xor)-G(xor))==1 && (X(xtar)-G(xtar))==1 && Sigma_u(sigma)==1 && xor~=xtar
            
            
            adj{xor}(find(adj{xor}==xtar,1))=[]; %removes first element
            adjrev{xtar}(find(adjrev{xtar}==xor,1))=[];
            if Sigma_u(sigma)
                adj_urev{xtar}(find(adj_urev{xtar}==xor,1))=[];
            end

             [Yprime,Gprime] = DLSS(X, adj, adjrev, adj_urev,(Y+X0)>0,(G+Xm)>0);
            
        else
            Yprime = Y; Gprime=G;
        end

    elseif skipCalc==true && doGroup==false
        xor=transdelta(1);
        sigma=transdelta(2);
        xtar=transdelta(3);
        
        if Y(xor)==1 && Y(xtar)==1 && xor~=xtar
            skipped = true; Yprime = Y; Gprime=G;
             
        elseif (G(xor)-Y(xor))==1 && G(xtar)==1 && xor~=xtar
            skipped = true; Yprime = Y; Gprime=G;

        elseif (X(xor)-G(xor))==1 && (X(xtar)-G(xtar))==1 && Sigma_u(sigma)==1 && xor~=xtar
            skipped = true; Yprime = Y; Gprime=G;
            
        else
            Yprime = Y; Gprime=G;
        end
        
    elseif skipCalc==false && doGroup==true
        %perform grouped FRS, BRS, or SS call
        
        Gr1=[];
        Gr2=[];
        Gr3=[];
        
        for ti = transminind
            xor=transX(1,ti);
            sigma=transX(2,ti);
            xtar=transX(3,ti);
            
            if Y(xor)==1 && Y(xtar)==1 && xor~=xtar
                Gr1=[Gr1 ti];
            elseif (G(xor)-Y(xor))==1 && G(xtar)==1 && xor~=xtar
                Gr2=[Gr2 ti];
            elseif (X(xor)-G(xor))==1 && (X(xtar)-G(xtar))==1 && Sigma_u(sigma)==1 && xor~=xtar
                Gr3=[Gr3 ti];
            else
                0;
            end
        end
        
        if size(Gr1,2)>0

            for i=Gr1
                xor=transX(1,i);
                sig=transX(2,i);
                xtar=transX(3,i); 
                adj{xor}(find(adj{xor}==xtar,1))=[]; %removes first element
                adjrev{xtar}(find(adjrev{xtar}==xor,1))=[];
                if Sigma_u(sig)
                    adj_urev{xtar}(find(adj_urev{xtar}==xor,1))=[];
                end
            end

             [Yprime,Gprime] = DLSS(G, adj, adjrev, adj_urev, X0, Xm);
            deltaCalc=Gr1;
            
        elseif size(Gr2,2)>0

            for i=Gr2
                xor=transX(1,i);
                sig=transX(2,i);
                xtar=transX(3,i); 
                adj{xor}(find(adj{xor}==xtar,1))=[]; %removes first element
                adjrev{xtar}(find(adjrev{xtar}==xor,1))=[];
                if Sigma_u(sig)
                    adj_urev{xtar}(find(adj_urev{xtar}==xor,1))=[];
                end
            end 

             [Yprime,Gprime] = DLSS(G, adj, adjrev, adj_urev, (Y+X0)>0,(Y+Xm)>0);
            deltaCalc=Gr2;
            
        elseif size(Gr3,2)>0

            for i=Gr3
                xor=transX(1,i);
                sig=transX(2,i);
                xtar=transX(3,i); 
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
        
            [Yprime,Gprime] = DLSS(X, adj, adjrev, adj_urev,(Y+X0)>0,(G+Xm)>0);
            deltaCalc=Gr3;
            
        else
            skipped = true; Yprime = Y; Gprime=G;
        end
        
    end
end