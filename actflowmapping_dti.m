function [r_overall,  taskActualMatrix2, taskPredMatrix2, r_bytask, p_bytask, taskActualMatrix, taskPredMatrix, r_bysubj, r_avgfirst_bytask, r_avgfirst_mean] = actflowmapping_dti(taskActMatrix, connMatrix)

%actflowmapping_dti(taskActVector, connMatrix)
%
%Input variables:
%taskActMatrix - A regionXtaskXsubject matrix with activation levels
%connMatrix - A regionXregionXstateXsubject matrix with connectivity values (e.g., DTI matrix)
%
%This function takes a set of task-evoked activations (e.g., fMRI GLM beta activations) and a set of connections (e.g., DTI connectivity) and applies the "activity flow mapping" procedure. 
%See Cole et al. (2016; Nature Neuroscience) for more information. 
%
%Output variables:
%r_overall - The r-value comparing predicted to actual activation patterns averaged across all tasks
%p_overall - The p-value comparing predicted to actual activation patterns averaged across all tasks. This p-value is based on an across-subject paired t-test of the Fisher's z-transformed across-task average r-values.
%t_overall - The t-value associated with the p_overall value. Degrees of freedom = number of subjects - 1.
%r_bytask - The r-value comparing predicted to actual activation patterns for each task separately.
%p_bytask - The p-value associated with r_by task. This p-value is based on an across-subject paired t-test of the Fisher's z-transformed r-values for each task separately.
%taskActualMatrix - The actual task-evoked activation patterns. Z-normalized.Format: regionXtaskXsubject.
%taskPredMatrix - The predicted task-evoked activation patterns via the activity flow mapping approach.Z-normalized.Format: regionXtaskXsubject.
%taskActualMatrix2 - The actual task-evoked activation patterns. Format: regionXtaskXsubject.
%taskPredMatrix2 - The predicted task-evoked activation patterns via the activity flow mapping approach.Format: regionXtaskXsubject.
%r_bysubj - The r-value comparing predicted to actual activation patterns, separately for each task and each subject. Format: taskXsubject.
%r_avgfirst_bytask - The r-values (separately for each task) comparing predicted to actual activation patterns, computed after averaging the predicted activations across subjects and the actual activations across subjects. This tends to produce more accurate predictions, possibly due to higher signal-to-noise through averaging across subjects.
%r_avgfirst_mean - The across-task average of r_avgfirst_bytask.
%
%
%Adapted from Michael W. Cole (doi:10.1038/nn.4406)


numTasks=size(taskActMatrix,2);
numRegions=size(taskActMatrix,1);
numConnStates=size(connMatrix,3);
numSubjs=size(connMatrix,4);

%Setup for prediction
taskPredMatrix=zeros(numRegions,numTasks,numSubjs);
taskPredRs=zeros(numTasks,numSubjs);
taskActualMatrix=taskActMatrix;
regionNumList=1:numRegions;

for subjNum=1:numSubjs
    for taskNum=1:numTasks

        %Get this subject's activation pattern for this task
        taskActVect=taskActMatrix(:,taskNum,subjNum);

        for regionNum=1:numRegions

            %Hold out region whose activity is being predicted
            otherRegions=regionNumList;
            otherRegions(regionNum)=[];

            %Get this region's connectivity pattern
            if numConnStates > 1
                stateFCVect=connMatrix(:,regionNum,taskNum,subjNum);
            else
                %If using resting-state (or any single state) data
                stateFCVect=connMatrix(:,regionNum,1,subjNum);
            end

            %Calculate activity flow prediction
            taskPredMatrix(regionNum,taskNum,subjNum)=sum(taskActVect(otherRegions).*stateFCVect(otherRegions));

        end

        %Normalize values (z-score)
       
        %%
      taskActualMatrix2=taskActMatrix;
      taskPredMatrix2=taskPredMatrix;
        taskPredMatrix(:,taskNum,subjNum)=(taskPredMatrix(:,taskNum,subjNum)-mean(taskPredMatrix(:,taskNum,subjNum)))./std(taskPredMatrix(:,taskNum,subjNum));
        taskActualMatrix(:,taskNum,subjNum)=(taskActMatrix(:,taskNum,subjNum)-mean(taskActMatrix(:,taskNum,subjNum)))./std(taskActMatrix(:,taskNum,subjNum));

        %Calculate predicted to actual similarity for this task
        r=corrcoef(taskPredMatrix(:,taskNum,subjNum),taskActualMatrix(:,taskNum,subjNum));
        taskPredRs(taskNum,subjNum)=r(1,2);
    end
end

%Calculate average r, across-subject p-value
r_bytask=tanh(mean(atanh(taskPredRs),2));
p_bytask=ones(numTasks,1);
for taskNum=1:numTasks
    [~, p_bytask(taskNum)]=ttest(atanh(taskPredRs(taskNum,:)));
end
r_overall=tanh(mean(mean(atanh(taskPredRs),1),2));
[~, p_overall, ~, stats]=ttest(mean(atanh(taskPredRs(taskNum,:)),1));
%By subj
r_bysubj=taskPredRs;
t_overall=stats.tstat;

%Calculate average-then-compare results
r_avgfirst_bytask=zeros(numTasks,1);
for taskNum=1:numTasks
    r_avgfirst_bytask(taskNum)=corr(mean(taskPredMatrix(:,taskNum,:),3),mean(taskActualMatrix(:,taskNum,:),3));
end
r_avgfirst_mean=tanh(mean(atanh(r_avgfirst_bytask)));
