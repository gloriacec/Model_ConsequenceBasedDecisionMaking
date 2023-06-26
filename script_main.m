clearvars
%% choose horizon
disp('What horizon would you like to simulate?')
disp('Possible values 0,1, or 2 ')
prompt='nH=';
nH=input(prompt);
while ~ismember(nH,0:2)
    disp('You must insert either 0, 1, or 2')
    nH = input(prompt);
end

%% select alpha parameter
disp('Select the alpha parameter (within range [-0.03,0])')
disp('to tune the visual discrimination')
prompt='alpha=';
alpha=input(prompt);
while alpha<-0.03 || alpha>0
    disp('Alpha must be within range [-0.03,0]')
    alpha = input(prompt);
end
beta=0.055-2.5*alpha;

%% load stimuli
[stimuli,minFeedback,maxFeedback,lAlB,nT,nE]=f_loadStimuli(nH,alpha,beta);

%% time constant tau
disp('Time constant tau in Eq.1 and Eq.3')
disp('Range [10,100]')
prompt='tau=';
tau=input(prompt);
while tau<10 || tau>100
    disp('Tau must be within range [10,100]')
    tau = input(prompt);
end

%% decision threshold: difference between firing rates
disp('Decision threshold: difference between firing rates')
disp('Delta in Eq.1 and Eq.3. Range [0.01,0.035]')
prompt='Delta=';
th_decision=input(prompt);
while th_decision<0.01 || th_decision>0.035
    disp('Delta must be within range [0.01,0.035]')
    th_decision = input(prompt);
end

%% Noise for the firing rates
disp('Noise for the firing rates')
disp('sigma in Eq.1 and Eq.3. Range [0.001,0.01]')
prompt='sigma=';
sigma=input(prompt);
while sigma<0.001 || sigma>0.01
    disp('sigma must be within range [0.001,0.01]')
    sigma = input(prompt);
end

%% initial bias phi_0
disp('What is the initial bias? phi0 in Eq.4')
disp('Enter vector of dimension nH+1 with values in [0.3,0.7]')
disp('Values <0.5 indicate bias towards choosing small')
prompt='phi0=';
phi0=input(prompt);
while length(phi0)~=nH+1
    disp('Vector dimension must be nH+1')
    phi0=input(prompt);
end
while sum(phi0<0.3 | phi0>0.7)>0
    disp('each value of the initial bias must be within range [0.3,0.7]')
    phi0 = input(prompt);
end

%% Learning rate
disp('Learning rate k in Eq.4. Range [0,3]')
prompt='k=';
k=input(prompt);
while k<0 || k>3
    disp('sigma must be within range [0,3]')
    k = input(prompt);
end

%% Decisional uncertainty
disp('Decisional uncertainty sigma_psi in Eq.2. Range [0.2,1]')
prompt='sigma_psi=';
sigma_psi=input(prompt);
while sigma_psi<0.2 || sigma_psi>1
    disp('sigma must be within range [0.2,1]')
    sigma_psi = input(prompt);
end

%% Number of iterations
disp('The model is ready to run')
disp('How many iterations do you want to run?')
prompt='num_iter=';
num_iter=input(prompt);

%% run the model
script_runModel
disp('Done')

%% Plot some results
figure(1);clf;
dim_FigTot=[60 60 800 400];
set(gcf,'Position',dim_FigTot)
tl=tiledlayout(1,2,'Padding','none','TileSpacing','tight');

% Reaction time
nexttile(1)
ax=gca;
hold on
x=-250:2000;
cdf_x=x(2:end);
RT_iterFlat=reshape(RT_iter,[],1);
RT_cdf=histcounts(RT_iterFlat,x,'Normalization','cdf');
plot(cdf_x/1000,RT_cdf,'b-','LineWidth',2,'DisplayName','Model')

ylim([0,1])
ylabel('CDF')
xlabel('RT(s)')
box off; ax.FontSize=14;


% Performance
nexttile(2)
ax=gca;
hold on

perf=perf_iter;

sp=1:nE; % sample points
sw=0; % sample window
mm=nan(length(sp),1);
sd=nan(length(sp),1);
mm_p=nan(length(sp),1);
for i=1:length(sp)
    ind=max(sp(i)-sw,1):min(sp(i)+sw,size(perf,2));
    mm(i)=mean(perf(:,ind),'all');
    sd(i)=std(perf(:,ind),0,'all');
end
barra=0.5;
xshade=[sp sp(end:-1:1)];
yshade=[mm-barra*sd; mm(end:-1:1)+barra*sd(end:-1:1)];
shade=fill(xshade,yshade,'b');
shade.FaceAlpha=0.2;
shade.EdgeColor = 'none';     
plot(sp,mm,'b','LineWidth',2,'DisplayName','Model');
ylim([0,1])
ylabel('Performance')
xlabel('Episode')
box off; ax.FontSize=14;
%%