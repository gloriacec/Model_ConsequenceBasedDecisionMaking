function [stimuli,minFeedback,maxFeedback,lAlB,nT,nE]=f_loadStimuli(nH,alpha,beta)
% load stimuli
block=1; % block within horizon

stimuli=load(['hor',num2str(nH),'_block',num2str(block),'.txt']);
stimuli=stimuli(:,2:end);

stim=nan(size(stimuli));
temp1=zeros(size(stimuli)./[1,2]);
for i=1:size(temp1,2)
    ind=(1:2)+2*(i-1);
    temp1(:,i)=mean(stimuli(:,ind),2);
    stim(:,ind)=stimuli(:,ind)-temp1(:,i)+0.5;
end
lAlB=stim.*beta+alpha;

% stim difficulty
% stimDiff=abs(diff(stim(:,1:2),1,2));

nT=size(lAlB,1); % number of total trials
nE=nT/(nH+1); % number of episodes

% max reward per episode
allFeedback=zeros(nE,2^(nH+1));
for h=1:2^(nH+1)
    for t=1:nH+1
        l=ceil(h/2^(nH+1-t));
        allFeedback(:,h)=allFeedback(:,h)+stimuli(t:nH+1:end,l);
    end
end
minFeedback=min(allFeedback,[],2);
maxFeedback=max(allFeedback,[],2);
end