%% fixed parameters
strategy=1;
% 1: delta= (rew(T)-mean(status(T,:)))*(phi(E,TE)^2)*(phi(E,TE)-1)^2;
% 2: delta= F(E)-1
% 3: delta= (F(E)^4)*sign(rew(T)-mean(status(T,:)))*(phi(E,TE)^2)*(phi(E,TE)-1)^2;

dt=0.1; % integration step
time_trial=2000; % length of a trial in ms
t_think=500; % time before stimula, i.e., time to think which strategy to use
timeX=dt:dt:time_trial;

w_p=1.4; % strength of recurrent connection 
w_m=1.5; % strength of inhibitory connection between pools
F_max=40e-03;
kk=22e-03;
theta=15e-03;
f=@(x) F_max./(1+exp(-(x-theta)./kk));
tau_psi=10; % constant time for equation 2

% Decision time:
% 1. rA>rB+th
% 2. rA>rB+th for at least 100ms 
dec_time=1;

%% run
perf_iter=  nan(num_iter,nE);
phi_iter=   nan(num_iter,nE,nH+1);
psi_iter=   nan(num_iter,nT);
choice_iter=nan(num_iter,nT);
intended_iter=nan(num_iter,nT);
trErr_iter= nan(num_iter,nT);
RT_iter=    nan(num_iter,nT);

parfor iter=1:num_iter
    % initial conditions
    psi_T=nan(nT,1);
    time_decisions=nan(nT,1);
    choice=nan(nT,1);
    reward=nan(nT,1);
    rew=nan(nT,1);
    status=nan(nT,2);
    intended=nan(nT,1);
    trErr=zeros(nT,1);
    T=1; % trial
    l_A=lAlB(T,1);
    l_B=lAlB(T,2);
    status(T,:)=[stimuli(T,1),stimuli(T,2)];
    phi=nan(nE,nH+1);
    phi(1,:)=phi0;
    rewUpdate=nan(nE,nH+1);
    F=zeros(nE,1);

    % equations
    for E=1:nE
        for TE=1:nH+1
            decision=0;

            psi=nan(1,length(timeX));
            psi(1)=phi(E,TE);
            noise=nan(1,length(timeX));

            r_A=nan(length(timeX),1);
            r_B=nan(length(timeX),1);

            for t=1:length(timeX)-1
                % psi
                noise(t)=sigma_psi*randn*(1/(timeX(t)+1)^2);
                psi(t+1)=psi(t)+(dt/tau_psi)*(-4.*psi(t).*(psi(t)-1).*(psi(t)-0.5))+(1/tau_psi)*noise(t);

                % firing rates
                if ismember(t,0:t_think-1)
                    r_A(t)=0;
                    r_B(t)=0;
                end
                r_A(t+1)=r_A(t)+(dt/tau)*(-r_A(t)+f((1-psi(t))*l_B+psi(t)*l_A+w_p*r_A(t)-w_m*r_B(t)))+(1/tau)*sigma*randn;
                r_B(t+1)=r_B(t)+(dt/tau)*(-r_B(t)+f((1-psi(t))*l_A+psi(t)*l_B+w_p*r_B(t)-w_m*r_A(t)))+(1/tau)*sigma*randn;

                % is the decision taken?
                if decision==0 && abs(r_A(t)-r_B(t))>th_decision && t>(100/dt)
                    if dec_time==2 && all(abs(r_A((t-100/dt):t)-r_B((t-100/dt):t))>th_decision)
                        decision=1;
                        time_decisions(T)=t;
                        [~,choice(T)]=max([r_A(t),r_B(t)]); % choice= selected stimulus
                        reward(T)=status(T,choice(T));
                    elseif dec_time==1
                        decision=1;
                        time_decisions(T)=t;
                        [~,choice(T)]=max([r_A(t),r_B(t)]); % choice= selected stimulus
                        reward(T)=status(T,choice(T));
                    end
                end
            end

            psi_T(T)=round(psi(t));

            % if the trial ends and the decision is not taken 
            if decision==0
                time_decisions(T)=t;
                choice(T)=randi(2);
                reward(T)=0;
                trErr(T)=1;
            end

            % the selected stimulus was the one I was aiming for?
            [~,bigS]=max(status(T,:));
            if psi_T(T)==1
                intended(T)=(choice(T)==bigS);
            else
                intended(T)=(choice(T)~=bigS);
            end

            F(E)=F(E)+reward(T);

            % initialization for next trial
            % AND pre-frontal: elaborate consequence of chosen choise

            if TE==nH+1 % stimuli selection for next trial
                if E~=nE
                    l_A=lAlB(T+1,1);
                    l_B=lAlB(T+1,2);
                    status(T+1,:)=[stimuli(T+1,1),stimuli(T+1,2)];
                end
                rew(T)=reward(T);
            else
                % ######################
                % next=(2^TE)*(choice(T+1-TE)-1); % if TE==1
                % next=(2^TE)*(choice(T+1-TE)-1)+(2^(TE-1))*(choice(T+2-TE)-1); % if TE==2
                % next=(2^TE)*(choice(T+1-TE)-1)+(2^(TE-1))*(choice(T+2-TE)-1)+(2^(TE-2))*(choice(T+3-TE)-1); % if TE==3
                % ######################
                ind=(2^TE)*(choice(T+1-TE)-1);
                for tt=2:TE
                    ind=ind+(2^(TE+1-tt))*(choice(T+tt-TE)-1);
                end

                l_A=lAlB(T+1,ind+1);
                l_B=lAlB(T+1,ind+2);
                status(T+1,:)=[stimuli(T+1,ind+1),stimuli(T+1,ind+2)];
                rew(T)=mean(status(T+1,:));
            end

            % reward update
            rewUpdate(E,TE)=sign(rew(T)-mean(status(T,:)))*0.2;
            if strategy==1
                phi(E+1,TE)=phi(E,TE);
                phi(E+1,TE)=phi(E,TE)+k*(2*psi_T(T)-1)*rewUpdate(E,TE)*(phi(E,TE)^2)*(phi(E,TE)-1)^2;

            end

            T=T+1; % trial counter
        end

        if prod(reward(T-(nH+1):T-1))==0
            F(E)=0;
        else
            F(E)=(F(E)-minFeedback(E))/(maxFeedback(E)-minFeedback(E)); % normalize feedback by the max reward
        end
        if strategy==2
            temp=T-(nH+1):T-1;
            for tt=1:nH+1
                phi(E+1,tt)=phi(E,tt);
                phi(E+1,tt)=phi(E,tt)+k*(2*psi_T(temp(tt))-1)*(F(E)-1);
            end
        elseif strategy==3
            temp=T-(nH+1):T-1;
            for tt=1:nH+1
                phi(E+1,tt)=phi(E,tt);
                phi(E+1,tt)=phi(E,tt)+k*(F(E)^4)*(2*psi_T(temp(tt))-1)*sign(rew(temp(tt))-mean(status(temp(tt),:)))*(phi(E,tt)^2)*(phi(E,tt)-1)^2;
            end
        end

    end
    phi=phi(1:end-1,:);

%   ###############################
    perf_iter(iter,:)=F;
    phi_iter(iter,:,:)=phi;
    psi_iter(iter,:)=psi_T;
    choice_iter(iter,:)=choice;
    intended_iter(iter,:)=intended;
    trErr_iter(iter,:)=trErr;
    RT_iter(iter,:)=time_decisions.*dt;
end

%%