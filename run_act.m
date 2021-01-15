function [PersonAttr] = run_act( taskActMatrix,connMatrix )
%run_act( taskActMatrix,connMatrix )
%Input variables:
%taskActMatrix - A regionXstateXsubject matrix with activation levels
%connMatrix - A regionXregionXstateXsubject matrix with connectivity values (e.g., DTI matrix)
%
%This function applies the "activity flow mapping" procedure on networks.
%Output variables:
%PersonAttr - The r-value comparing predicted to actual activation patterns
%for each network

region = size(taskActMatrix,1);
state = size(taskActMatrix,2);
subject =size(taskActMatrix,3);

PersonAttr=[];
%Get predictions
[r_overall, taskActualMatrix2, taskPredMatrix2, r_bytask, p_bytask, taskActualMatrix, taskPredMatrix, r_bysubj, r_avgfirst_bytask, r_avgfirst_mean] = actflowmapping_dti(taskActMatrix, connMatrix);
Actual=taskActualMatrix;
Pred=taskPredMatrix;
Regionp=[];
Regionr=[];
for i=1:region
    [r,p]=corrcoef(Pred(i,:),Actual(i,:))
    Regionp=[Regionp,p(1,2)];
    Regionr=[Regionr,r(1,2)];
end
 %%
        Att=[5 7 8 9 10 11 12 13 14 15 16 19 61 62 65 66 83 89];
        DMN=[3 4 6 23 24 25 26 27 28 31 32 35 36 67 68 85 86 90];
        Sen=[1 2 17 18 20 29 30 57 58 59 60 63 64 69 70 79 80 81 82 84];
        Sub=[21 22 33 34 37 38 39 40 41 42 71:78 87 88];
        Vis=[43:56];
        
%%  
Actual2=taskActualMatrix2;
Pred2=taskPredMatrix2;
%%
ActualAtt=zscore(Actual2(Att,:));
PredAtt=zscore(Pred2(Att,:));

for m=1:subject
    [r,p]=corrcoef(PredAtt(:,m),ActualAtt(:,m));
    PersonAttr=[PersonAttr,r(1,2)];
end

end

