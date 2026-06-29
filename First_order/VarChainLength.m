CLList=[2 5 10 20 50 100];
for nrxni=1:length(CLList)
    nrxn=CLList(nrxni);
    clear SMat fVec fccMat dgMat kGood KGood giniVec dgMax
    RandSample;
    subplot(2,3,nrxni);
    %dscatter((dgMat(:,1)),fccMat(:,1));
    dscatter(dgMax(:),fccMat(:,1));
    xlabel('deltaG of R1');
    ylabel('FCC of R1');
    title(strcat('NRxn=',num2str(nrxn)));
end
colormap('redbluecmap');