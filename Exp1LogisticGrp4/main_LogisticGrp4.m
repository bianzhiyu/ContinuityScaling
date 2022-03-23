clear all;

MU21=[0,0.05,0.10,0.15];
Results=cell(4,1);
for ExpInd=1:4
    %%
    disclen=400;
    usedlen=5000;
    VarNum=2;
    rawData=zeros(VarNum,disclen+usedlen);
    rawData(:,1)=[cos(0.9);sin(1)];
    r=[3.8;3.7];
    mu=[0,0;  ...
        MU21(ExpInd),0];
    
    for i=1:(disclen+usedlen-1)
        if i<disclen/2
            rawData(:,i+1)=rawData(:,i).*(r-r.*rawData(:,i));
        else
            rawData(:,i+1)=rawData(:,i).*(r-r.*rawData(:,i)-mu*rawData(:,i));
        end
    end
    rawData=rawData(:,(disclen+1):(disclen+usedlen));
    %%
    TSs=rawData';
    %%
    embedDims=[3,3];
    embedLags=[1,1];
    epsilonPercents=exp(linspace(log(0.001),log(1),33))';
    eliminateSpan=[0;0];
    option=struct('skipCInd2',true);
    [CInd,CInd2local,CInd2global,maxt,M,Dist,MaxDist,epsilon,delta,InNeighborCount,GlobalDelta]=...
        ACX(TSs,embedDims,embedLags,epsilonPercents,eliminateSpan,option);
    
    [~,VarNum]=size(epsilon);
    TotalPermutTimes=25;
    CIndWhite=zeros(VarNum,VarNum,TotalPermutTimes);
    PTData=cell(TotalPermutTimes,1);
    option2=option;
    option2.permute=true;
    option2.groupNum=25;
    parfor PermutInd=1:TotalPermutTimes
        fprintf('PermutInd=%d\n',PermutInd);
        [CIndPT,~,~,PTData{PermutInd}.maxt,PTData{PermutInd}.M,DistPT,MaxDistPT,...
            PTData{PermutInd}.epsilon,PTData{PermutInd}.delta,...
            PTData{PermutInd}.InNeighborCount,PTData{PermutInd}.GlobalDelta]=...
            ACX(TSs,embedDims,embedLags,epsilonPercents,eliminateSpan,option2); 
%[CInd,CInd2local,CInd2global,maxt,M,Dist,MaxDist,epsilon,delta,InNeighborCount,GlobalDelta]=ACX(TSs,embedDims,embedLags,epsilonPercents,eliminateSpan,option)
        CIndWhite(:,:,PermutInd)=CIndPT;
    end
    CIndWhite(:,:,TotalPermutTimes+1)=CInd;
    TotalPermutTimes=TotalPermutTimes+1;
    
    Mn=zeros(VarNum,VarNum);
    StdDevi=zeros(VarNum,VarNum);
    Pval=zeros(VarNum,VarNum);
    for ic=1:VarNum
        for ir=1:VarNum
            if ic==ir
                continue;
            end
            tmppool=squeeze(CIndWhite(ic,ir,:));
            tmppool(length(tmppool)+1)=CInd(ic,ir);
            Mn(ic,ir)=mean(tmppool);
            StdDevi(ic,ir)=std(tmppool);
            Pval(ic,ir)=1-normcdf((CInd(ic,ir)-Mn(ic,ir))/StdDevi(ic,ir));
        end
    end
    Results{ExpInd}.mu=mu;
    Results{ExpInd}.TSs=TSs;
    Results{ExpInd}.maxt=maxt;
    Results{ExpInd}.M=M;
    Results{ExpInd}.epsilon=epsilon;
    Results{ExpInd}.delta=delta;
    Results{ExpInd}.GlobalDelta=GlobalDelta;
    Results{ExpInd}.InNeighborCount=InNeighborCount;
    Results{ExpInd}.CInd=CInd;
    Results{ExpInd}.Mn=Mn;
    Results{ExpInd}.StdDevi=StdDevi;
    Results{ExpInd}.Pval=Pval;
    Results{ExpInd}.CIndWhite=CIndWhite;
    Results{ExpInd}.PTData=PTData;
end
save('sav.mat','Results');
%%
clear all;
load('sav.mat');
fprintf('Output:\n')
for ExpInd=1:4
    fprintf('ExpInd=%d, 1c2=%.5f, 2c1=%.5f, mu21=%.2f, mu12=%.2f, p-val 1c2=%.4f, p-val 2c1=%.4f\n',...
        ExpInd,Results{ExpInd}.CInd(1,2),Results{ExpInd}.CInd(2,1),...
        Results{ExpInd}.mu(2,1),Results{ExpInd}.mu(1,2),...
        Results{ExpInd}.Pval(1,2),Results{ExpInd}.Pval(2,1));
end

plotoptions={'k-+','m--*','b:^','r-.d','g-s'};
lgd=cell(2,2,4);

FH=figure(1);
hold on;
box on;
set(gcf,'unit','centimeter','position',[2,2,18,15],'Color',[1,1,1]);
set(gca,'position',[0.2,0.2,0.7,0.7]);
set(gca,'FontSize',30,'FontName','Times New Roman');
for ExpInd=1:4
    plot(log(Results{ExpInd}.epsilon(:,2)),Results{ExpInd}.GlobalDelta(:,1,2),plotoptions{ExpInd}...
        ,'LineWidth',2,'MarkerSize',15);
    lgd{1,2,ExpInd}=sprintf('$\\mu_{21}=%.2f$',Results{ExpInd}.mu(2,1));
end
axis([-8,1,0,0.7]);
set(legend(lgd{1,2,1},lgd{1,2,2},lgd{1,2,3},lgd{1,2,4}),...
    'FontSize',20,...
    'Interpreter','latex',...
    'location','southeast');
set(xlabel('$\ln \varepsilon_{x_2}$'),...
    'Interpreter','latex');
set(ylabel('$\langle\delta^t_{x_1}\rangle_t$'),...
    'Interpreter','latex');
savefig(FH,'E1log_1c2.fig');
print(FH,'E1log_1c2.png','-dpng');
print(FH,'E1log_1c2.eps','-depsc');

FH=figure(2);
hold on;
box on;
set(gcf,'unit','centimeter','position',[2,2,18,15],'Color',[1,1,1]);
set(gca,'position',[0.2,0.2,0.7,0.7]);
set(gca,'FontSize',30,'FontName','Times New Roman');
for ExpInd=1:4
    plot(log(Results{ExpInd}.epsilon(:,1)),Results{ExpInd}.GlobalDelta(:,2,1),plotoptions{ExpInd}...
        ,'LineWidth',2,'MarkerSize',15);
    lgd{2,1,ExpInd}=sprintf('$\\mu_{12}=%.2f$',Results{ExpInd}.mu(1,2));
end
axis([-8,1,0,0.7]);
set(legend(lgd{2,1,1},lgd{2,1,2},lgd{2,1,3},lgd{2,1,4}),...
    'FontSize',20,...
    'Interpreter','latex',...
    'location','southeast');
set(xlabel('$\ln \varepsilon_{x_1}$'),...
    'Interpreter','latex');
set(ylabel('$\langle\delta^t_{x_2}\rangle_t$'),...
    'Interpreter','latex');
savefig(FH,'E1log_2c1.fig');
print(FH,'E1log_2c1.png','-dpng');
print(FH,'E1log_2c1.eps','-depsc');



