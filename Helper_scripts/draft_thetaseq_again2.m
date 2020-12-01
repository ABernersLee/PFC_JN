function newmat = draft_thetaseq_again2(thisdir,ach,armposindex,curp,oc,pc)
EstBin = .02;
[MatIn,MatOut,Index,binspike] = make_PosteriorfromThetaCandEvents(thisdir,[ach(:,1)-EstBin/2 ach(:,1)+.125+EstBin/2],EstBin);

Matrix = MatIn+MatOut;
Mat = Matrix./(ones(size(Matrix,1),1)*sum(Matrix,1));

rat = NaN(size(ach,1),3);
wind = 10; %was 15
posexcl = true(size(armposindex,1),1);
for iarm = 1:3
    posexcl(find(armposindex(:,iarm)==1,1,'first'):find(armposindex(:,iarm)==1,1,'first')+wind) = false;
    posexcl(find(armposindex(:,iarm)==1,1,'last')-wind:find(armposindex(:,iarm)==1,1,'last')) = false;
end
postouse = find(posexcl);
newmat = NaN(wind*2+1,length(4:4*(unique(diff(Index)))),size(ach,1));
for c = 1:size(ach,1)
    Ind = 4*Index(c)+1:4*Index(c+1)-3;
    b = binspike(:,Ind);
    rat(c,2) = sum(sum(b,2)>0);
    if rat(c,2)>4 && ismember(curp(c),postouse)
        if ~isnan(oc(c,2))
            newmat(:,:,c) = Mat(curp(c)-wind:curp(c)+wind,Ind);
        elseif ~isnan(pc(c,2))
            newmat(:,:,c) = Mat(curp(c)+wind:-1:curp(c)-wind,Ind);
        end
        
    end
end
