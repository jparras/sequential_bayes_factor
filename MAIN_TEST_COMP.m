%% Empirical test comparison
% Juan Parras, GAPS-UPM, September 2018
clear all; clc; close all;

%% Define parameters

% Common parameters
p0v=[0.1, 0.2, 0.3]; %Parameters to test
ptest=0:0.025:0.5;
ns=500; % Length of sequence (max test length also)
nmW=4; % MAx number of workers
nrep = 100; % Values for averaging

% UMP parameters
alpha_v=[0.01, 0.1];
nv=round(logspace(log10(5), log10(ns), 20)); % To obtain error

% BF parameters
weight=[1e3, 1e4];
dist=[0.01 0.05 0.1 0.2];
BF_limits=[1/3,3];

% For plotting purposes
leg_ump=cell(length(p0v), length(alpha_v));
leg_bf=cell(length(p0v),length(weight), length(dist));

%% Test UMP
ate_ump = zeros(length(p0v), length(alpha_v), length(nv));
error_plot_ump = zeros(length(p0v),length(alpha_v), length(nv), length(ptest));
for p0idx=1:length(p0v);
    p0=p0v(p0idx);
    for aidx = 1:length(alpha_v)
        alpha = alpha_v(aidx);
        leg_ump{p0idx,aidx}=['UMP, \alpha = ' num2str(alpha)];
        for nidx = 1:length(nv)
            n = nv(nidx);
            display(['Obtaining UMP: p0 = ' num2str(p0) ' alpha = ' num2str(alpha) ' and n ' num2str(n)])
            error=zeros(length(ptest),1);
            parfor (pidx=1:length(ptest), nmW)
            %for pidx=1:length(ptest)
                pt = ptest(pidx);
                for rep=1:nrep;  %For each repetition
                    x=binornd(1,pt,[1,n]); %Generate values
                    [dec,n_dec]=ump(x,p0,alpha,n);
                    if (dec==1 && pt<=p0) || (dec==0 && pt>p0)
                        error(pidx) = error(pidx) + 1/nrep;
                    end
                end
            end
            %Prepare output values
            error_plot_ump(p0idx, aidx, nidx,:) = error;
            ate_ump(p0idx, aidx, nidx) = mean(error);
        end
    end
end
%% Test BF

ate_bf = zeros(length(p0v),length(weight), length(dist));
arl_bf = zeros(length(p0v),length(weight), length(dist));
error_plot_bf = zeros(length(p0v), length(weight), length(dist), length(ptest));
arl_plot_bf = zeros(length(p0v), length(weight), length(dist), length(ptest));
for p0idx=1:length(p0v);
    p0=p0v(p0idx);
    for widx = 1:length(weight)
        for didx = 1:length(dist)
            w = weight(widx);
            d = dist(didx);
            display(['Obtaining BF: p0 = ' num2str(p0) ' w = ' num2str(w) ' and d ' num2str(d)])
            leg_bf{p0idx,widx, didx}=['BF, w = ' num2str(w) ',d = ' num2str(d)];
            error=zeros(length(ptest),1);
            arl=zeros(length(ptest),1);
            parfor (pidx=1:length(ptest), nmW)
                pt = ptest(pidx);
                for rep=1:nrep;  %For each repetition
                    x=binornd(1,pt,[1,ns]); %Generate values
                    f0=[w*p0 w*(1-p0)];
                    f1=[w*(p0+d) w*(1-(p0+d))];
                    [dec,n_dec]=bayes_factor(x,f0,f1, BF_limits);
                    arl(pidx) = arl(pidx) + n_dec/nrep;
                    if (dec==1 && pt<=p0) || (dec==0 && pt>p0)
                        error(pidx) = error(pidx) + 1/nrep;
                    end
                end
            end
            %Prepare output values
            error_plot_bf(p0idx, widx, didx,:) = error;
            arl_plot_bf(p0idx, widx, didx,:) = arl;
            ate_bf(p0idx, widx, didx) = mean(error);
            arl_bf(p0idx, widx, didx) = mean(arl);
        end
    end
end
%% Save data
save('bf_simulation')
%% Plot values
for p0vidx=1:length(p0v)
    figure();
    col = ['r','b','m'];
    for aidx = 1:length(alpha_v)
        loglog(nv,squeeze(ate_ump(p0vidx,aidx,:)), ['x-' col(aidx)], 'DisplayName', leg_ump{p0vidx,aidx}); grid on; hold all;
        legend('-DynamicLegend');
    end
    for widx = 1:length(weight)
        col = ['g','k','c'];
        for didx = 1:length(dist)
            shp=['o','*','s','^'];
            loglog(squeeze(arl_bf(p0vidx,widx,didx)),squeeze(ate_bf(p0vidx,widx,didx)), [col(widx) shp(didx)], 'DisplayName', leg_bf{p0vidx,widx,didx});
            hold all; grid on;
            legend('-DynamicLegend');
        end
    end
    axis([nv(1) nv(end) 0.01 1]);
    title(['p0 = ' num2str(p0v(p0vidx))])
    matlab2tikz(['bf_sim_ate_' num2str(p0vidx) '.tikz'],'showInfo',false,'height','\figureheight','width','\figurewidth')
    figure();
    for aidx = 1:length(alpha_v)
        for nidx = 1:length(nv)
            plot(ptest, squeeze(error_plot_ump(p0vidx,aidx, nidx,:)),'b');
            hold on; grid on;
        end
    end
    for widx = 1:length(weight)
        for didx = 1:length(dist)
            plot(ptest, squeeze(error_plot_bf(p0vidx,widx, didx,:)), 'r');
            hold on; grid on;
        end
    end
    title(['p0 = ' num2str(p0v(p0vidx))])
    matlab2tikz(['bf_sim_errp_' num2str(p0vidx) '.tikz'],'showInfo',false,'height','\figureheight','width','\figurewidth')
end