function [CInd,CInd2local,CInd2global,maxt,M,Dist,MaxDist,epsilon,delta,InNeighborCount,GlobalDelta]=...
    ACX(TSs,embedDims,embedLags,epsilonPercents,eliminateSpan,option)

[TSlen,VarNum]=size(TSs);
lenEps=length(epsilonPercents);
maxt=TSlen;
for i=1:VarNum
    if maxt>TSlen-(embedDims(i)-1)*embedLags(i)
        maxt=TSlen-(embedDims(i)-1)*embedLags(i);
    end
end
M=zeros(maxt,max(embedDims),VarNum);
for i=1:VarNum
    for j=1:maxt
        for k=1:embedDims(i)
            M(j,k,i)=TSs(j+(k-1)*embedLags(i),i);
        end
    end
end
%%%%%%%%%%
if isfield(option,'permute') && option.permute
    disp('Permute manifolds...\n');
    tmpM=zeros(size(M));
    GrpNum=option.groupNum;
    for ii=1:VarNum
        Grps=zeros(GrpNum+5,2);
        wp=0;
        Gplen=round(maxt/GrpNum);
        i=1;
        while i<maxt
            wp=wp+1;
            Grps(wp,1)=i;
            Grps(wp,2)=min(maxt,i+Gplen-1);
            i=Grps(wp,2)+1;
        end
        Grps=Grps(randperm(wp),:);
        tmpwp=1;
        for i=1:wp
            ingrpldiff=Grps(i,2)-Grps(i,1);
            tmpM(tmpwp:tmpwp+ingrpldiff,:,ii)=M(Grps(i,1):Grps(i,2),:,ii);
            tmpwp=tmpwp+ingrpldiff+1;
        end
    end
    M=tmpM;
end
%%%%%%%%%%%
Dist=zeros(maxt,maxt,VarNum);
for i=1:VarNum
    for j=1:maxt
        for k=1:maxt
            Dist(j,k,i)=norm(M(j,:,i)-M(k,:,i),2);
        end
    end
end
MaxDist=zeros(VarNum,1);
for i=1:VarNum
    MaxDist(i)=max(max(Dist(:,:,i)));
end
epsilon=zeros(lenEps,VarNum);
for i=1:VarNum
    for j=1:lenEps
        epsilon(j,i)=epsilonPercents(j)*MaxDist(i);
    end
end
delta=zeros(lenEps,maxt,VarNum,VarNum);
InNeighborCount=zeros(lenEps,maxt,VarNum,VarNum);
AllInd=1:maxt;
for ic=1:VarNum
    for ir=1:VarNum
        if ic==ir
            continue;
        end
        for it=1:maxt-1
            for ie=lenEps:-1:1
                veps=epsilon(ie,ir);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                logicInd=(Dist(it+1,:,ir)<=veps);
                logicInd( max(1,it+1-eliminateSpan(ir)):min(maxt,it+1+eliminateSpan(ir)))=false;
                ind=AllInd(logicInd)-1;
                ind=ind((ind>0) & (ind <=maxt));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if isempty(ind)
                    delta(ie,it,ic,ir)=delta(ie+1,it,ic,ir);
                    InNeighborCount(ie,it,ic,ir)=0;
                else
                    delta(ie,it,ic,ir)=mean(Dist(it,ind,ic));
                    InNeighborCount(ie,it,ic,ir)=length(ind);
                end
            end
        end
    end
end
GlobalDelta=zeros(lenEps,VarNum,VarNum);
for ic=1:VarNum
    for ir=1:VarNum
        if ic==ir
            continue;
        end
        for ie=1:lenEps
            GlobalDelta(ie,ic,ir)=mean(delta(ie,1:maxt-1,ic,ir));
        end
    end
end
CInd=zeros(VarNum,VarNum);
for ic=1:VarNum
    for ir=1:VarNum
        if ic==ir
            continue;
        end
        wp=0;
        slopeveccont=zeros(lenEps,2);
        for ie0=1:lenEps-1
            ie1=ie0+1;
            wp=wp+1;
            slopeveccont(wp,1)=(GlobalDelta(ie1,ic,ir)-GlobalDelta(ie0,ic,ir))/ ...
                (log(epsilon(ie1,ir))-log(epsilon(ie0,ir)));
            slopeveccont(wp,2)=ie0;
        end
        [~,stind]=sort(slopeveccont(1:wp,1),'descend');
        slopeveccont=slopeveccont(stind,1:2);
        udsplen=round(wp/2);
        selind=union(slopeveccont(1:udsplen,2),1+slopeveccont(1:udsplen,2));
        sele=epsilon(selind,ir);
        seld=GlobalDelta(selind,ic,ir);
        rgres=regress(seld,[log(sele),ones(size(sele))]);
        CInd(ic,ir)=rgres(1);
    end
end
if isfield(option,'skipCInd2') && option.skipCInd2
    CInd2local=0;
    CInd2global=0;
else
    CInd2local=zeros(maxt,VarNum,VarNum);
    CInd2global=zeros(VarNum,VarNum);
    for ic=1:VarNum
        for ir=1:VarNum
            if ic==ir
                continue;
            end
            for it=1:maxt-1
                for ie0=1:lenEps
                    if InNeighborCount(ie0,it,ic,ir)~=0
                        break;
                    end
                end
                wp=0;
                slopeveccont=zeros(lenEps,2);
                for ie1=ie0:lenEps-1
                    wp=wp+1;
                    ie2=ie1+1;
                    slopeveccont(wp,1)=(GlobalDelta(ie2,ic,ir)-GlobalDelta(ie1,ic,ir))/ ...
                        (log(epsilon(ie2,ir))-log(epsilon(ie1,ir)));
                    slopeveccont(wp,2)=ie1;
                end
                if wp>=2
                    [~,stind]=sort(slopeveccont(1:wp,1),'descend');
                    slopeveccont=slopeveccont(stind,1:2);
                    udsplen=round(wp/2);
                    selind=union(slopeveccont(1:udsplen,2),1+slopeveccont(1:udsplen,2));
                    sele=epsilon(selind,ir);
                    seld=GlobalDelta(selind,ic,ir);
                    rgres=regress(seld,[log(sele),ones(size(sele))]);
                    CInd2local(it,ic,ir)=rgres(1);
                else
                    CInd2local(it,ic,ir)=0;
                end
            end
            CInd2global(ic,ir)=mean(eliminateExtremeValue(CInd2local(1:maxt-1,ic,ir),0.985));
        end
    end
end
end