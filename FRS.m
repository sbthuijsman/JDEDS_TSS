function disc = FRS(X,adj,Xinit)

    disc = Xinit.*X; %discovered (logical array)
    curlyr = find(Xinit.*X); %current BFS layer (indices)

    %perform RS using adjacency list (with linear complexity)
    while ~isempty(curlyr)
    nxtlyr = []; %next BFS layer
        for u = curlyr 
            for v = adj{u}
                if disc(v)==0 && X(v)==1
                    disc(v)=1;
                    nxtlyr = [nxtlyr v];
                end
            end
        end
        curlyr=nxtlyr;
    end

end
