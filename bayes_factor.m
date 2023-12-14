function[dec,n]=bayes_factor(x,f0,f1, lim)
% f0 and f1 contain the beta parameters
lim=log(lim); %Convert to log units, to avoid overflow
n_hyp=size(f1,1); %Number of alt hypotheses

stop=0;
n=1;
S1=zeros(1,n_hyp+1); % S1(1) is the null hypothesis!
S2=zeros(1,n_hyp+1); % S2(1) is the null hypothesis!
S3=zeros(1,n_hyp+1); % S3(1) is the null hypothesis!
B=zeros(1,n_hyp);
output=-ones(1,n_hyp);

while stop==0
    sn=sum(x(1:n));
    % Update S3 term
    S3(1)=S3(1)+log(f0(1)+f0(2)+n-1);
    for i=1:n_hyp
        S3(i+1)=S3(i+1)+log(f1(i,1)+f1(i,2)+n-1);
    end
    % Update S1 and S2 terms
    if x(n)==1
        S1(1)=S1(1)+log(f0(1)+sn-1);
        for i=1:n_hyp
            S1(i+1)=S1(i+1)+log(f1(i,1)+sn-1);
        end
    else
        S2(1)=S2(1)+log(f0(2)+n-sn-1);
        for i=1:n_hyp
            S2(i+1)=S2(i+1)+log(f1(i,2)+n-sn-1);
        end
    end
    % Update Bayes factor
    B0=S1(1)+S2(1)-S3(1);
    for i=1:n_hyp
        B(i)=S1(i+1)+S2(i+1)-S3(i+1)-B0;
    end
    %Update output
    for i=1:n_hyp
        if B(i)<lim(1) && output(i)==-1
            output(i)=0; %Decide H0
        end
        if B(i)>lim(2) && output(i)==-1
            output(i)=1; %Decide H1
        end
    end
    % Decide to stop
    if isempty(find(output==1,1))==0 % An H1 has been decided
        stop=1;
        dec=1;
    elseif isequal(output,zeros(size(output))) % All decissions are H0
        stop=1;
        dec=0;
    else
        n=n+1;
    end
    if n==length(x) %End of sequence
        stop=1;
        % Decide using B
        if isempty(find(B>0, 1)) %% Decide H0: all tests have evidence to H0
            dec=0;
        else
            dec=1;
        end
    end
end
% ns=100;
% xp=linspace(0,1,ns);
% figure();
% plot(xp,betapdf(xp,f0(1),f0(2)),'b');
% hold on; grid on;
% for i=1:n_hyp
%     plot(xp,betapdf(xp,f1(i,1),f1(i,2)),'r');
% end
        