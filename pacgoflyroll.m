function out=pacgofly(in)
persistent w dt R iter wacc eacc interr f upalmacc derracc interacc usrcacc ...
    miubalast varba miuvarlast varvar prune grow miubamin stdbamin ...
    miuvarmin stdvarmin rules miuxelast miubacc miuvaracc rwin

X=[in(1);in(2);in(3)];  %input
n=length(X);            %n input
Xe=[X;1];               %extended input
tim=in(4);
yr=in(5);

%% initial condition
if tim==1
    %m=1;    %output
    R=1;    %initial number of rule
    %w=-0.01*ones(R,n+1);
    w=-0.01*ones(R,n+1);
    wacc=[];
    eacc=[];
    upalmacc=[];
    derracc=[];
    interacc=[];
    usrcacc=[];
    miubacc=[];
    miuvaracc=[];
    rules=[];
    rwin=[];
    dt=0.02;
    iter=1;
    interr=0;
    f=0.1*eye((n+1)*R,(n+1));
    varvar=0;
    varba=0;
    prune=0;
    grow=0;
    miubalast=0;
    miuvarlast=0;
    miuxelast=zeros(n+1,1);
end


%% ----distance calculation----%
distden=0;
for i=1:n+1
    distden=distden+sum(w(:,i).*w(:,i));
end
dist=abs((yr-w*Xe)./sqrt(1+distden));

%% ----degrees of membership calculation (fuzzyfication step)----%
for j=1:R
    miub(j)=exp(-0.5*(dist(j)/max(dist)));
end

%% ----Fuzzy Inference Rule---%
lambda=miub./sum(miub);
psij=lambda'*Xe';
% for j=1:R
%     psinor(j,:)=psij(j,:)./sum(psij(j,:));
% end

%% ----Upalm calculation (defuzzyfication step)----%
y=zeros(1,R);
for i=1:n+1
    for j=1:R
        y(j)=y(j)+(w(j,i)*psij(j,i)/sum(psij(j,:)));
    end
end

upalm=sum(y);

%% ----Usrc calculation (Sliding Mode Control)----%
% alpha1=0.01;
% alpha2=0.001;
% alpha3=0.00001;
% gamma1=alpha2/alpha1;             %sliding parameter 1
% gamma2=alpha3/alpha1;              %sliding parameter 2


err=in(1);
derr=in(2);
interr=interr+err*dt;
sl=5*(err+0.05*interr+0*derr);

% out=upalm+0.1*((2*sigmf(sl,[2 0]))-1);
usrc=100 *((2*sigmf(sl,[0.02 0]))-1);
out=usrc -1 * upalm;

%% Rule Significance Calculation
miuxe=miuxelast+(Xe-miuxelast)/iter;
RS=abs(w*miuxe);
[~,Rwin]=max(RS);

%%----weight update mechanism----%
% if grow==0 && prune==0
%     for j=1:R
%         for i=1:n+1
%             psit(((j-1)*(n+1)+i),1)=psij(j,i); % stacking psi
%         end
%     end
%
%     wdott=-0.1*f*psit*sl;
%     fdot=-1*f*(psit*psit')*f;
%     f=f+fdot*dt;
%
%     for j=1:R
%         for i=1:n+1
%             wdot(j,i)=wdott(((j-1)*(n+1)+i),1);   % unstacking weight
%         end
%     end
%
%     w=w+wdot.*dt;
% end

%% ----Winning Rule Update mechanism----%

psit=(psij(Rwin,:))'; % stacking psi

for i=1:n+1
    for j=1:n+1
        fwin(i,j)=f(((n+1)*(Rwin-1)+i),j);  %select f winner
    end
end

wdot=-0.1*fwin*psit*sl;
fdot=-1*fwin*(psit*psit')*fwin;
fwin=fwin+fdot*dt;

for i=1:n+1
    for j=1:n+1
        f(((n+1)*(Rwin-1)+i),j)=fwin(i,j); %save f update
    end
end

w(Rwin,:)=w(Rwin,:)+wdot'*dt; %update winning weight


%% Bias Calculation
Ey=sum(w*miuxe);
bias2=(yr-Ey)^2;
miuba=(miubalast+(bias2-miubalast))/iter;
varba=(varba+(bias2-miubalast)*(bias2-miuba))/iter;
stdba=sqrt(varba/iter);

if iter==1 || grow==0
    miubamin=miuba;
    stdbamin=stdba;
end

%% Rule Growing
Xi=1*(1.5*exp(-1*(bias2))+0.5);
RGleft=miuba+stdba;
RGright=miubamin+Xi*stdbamin;
if  RGleft >= RGright
    R=R+1;
    w(R,:)=(miuxe/max(abs(miuxe)))'.*mean(w(Rwin,:));
    %w(R,:)=w(Rwin,:);
    f=[f;fwin];
    %sumy=[sumy,0];
    grow=1;
    miubamin=miuba;
    stdbamin=stdba;
else
    grow=0;
end

miuxelast=miuxe;
miubalast=miuba;

%% Variance Calculation
Ey2=sum((w.^2)*(miuxe.^2));
variance=abs(Ey2-(Ey^2));
miuvar=(miuvarlast+(variance-miuvarlast))/iter;
varvar=(varvar+(variance-miuvarlast)*(variance-miuvar))/iter;
stdvar=sqrt(varvar/iter);

if iter>n+2 && R>=2
    if prune==0
        miuvarmin=miuvar;
        stdvarmin=stdvar;
    end
    
    %% Rule Pruning
    phii=0.4*(1.5*exp(-1*(variance))+0.5);
    left=miuvar+stdvar;
    right=miuvarmin+phii*stdvarmin;
    if left >= right
        miuvarmin=miuvar;
        stdvarmin=stdvar;
        [~,rulemin]=min(RS);
        index=1;
        for j=1:R
            if j~=rulemin
                wnew(index,:)=w(j,:);
                for a=1:(n+1)
                    for b=1:(n+1)
                        fnew((n+1)*(index-1)+a,b)=fwin(a,b);
                    end
                end
                index=index+1;
            end
        end
        w=wnew;
        R=R-1;
        f=fnew;
        %sumy=zeros(1,R);
        prune=1;
    else
        prune=0;
    end
    
end
miuvarlast=miuvar;

wacc(iter)=w(1,2);
eacc(iter)=err;
derracc(iter)=derr;
interacc(iter)=interr;
upalmacc(iter)=upalm;
usrcacc(iter)=usrc;
rmse=sqrt((eacc*eacc')/length(eacc));
rules(iter)=R;
miubacc(iter)=miuba;
miuvaracc(iter)=miuvar;
rwin(iter)=Rwin;
assignin('base','weight',w);
assignin('base','wacc',wacc);
assignin('base','eacc',eacc);
assignin('base','rmse',rmse);
assignin('base','upalmacc',upalmacc);
assignin('base','usrcacc',usrcacc);
assignin('base','derracc',derracc);
assignin('base','interacc',interacc);
% assignin('base','gamma1',gamma1);
% assignin('base','gamma2',gamma2);
assignin('base','facc',f);
assignin('base','rules',rules);
assignin('base','miubacc',miubacc);
assignin('base','miuvaracc',miuvaracc);
assignin('base','Rwin',rwin);
iter=iter+1;



