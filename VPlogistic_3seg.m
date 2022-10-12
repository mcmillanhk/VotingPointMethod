function [RCU, qs005, qs05, qs15, qs25, qs35, qs45, qs55, qs65, qs75, qs85, qs95, qs995]= VPlogistic_3seg(station)
% Voting point methodology of McMillan and Westerberg (2015) with logistic likelihood
% Created by David Ocio in May 2016
% Note that compared to the original version in the 2015 paper, the MCMC algorithm is different, the estimation 
% of the priors is simplified and all c parameters are set to 0 (eqn 5 in McMillan and Westerberg, 2015).
%
% ***** If using this code, please cite the following two papers *****
% McMillan, H. K., Westerberg, I. K. 2015. Rating curve estimation under epistemic uncertainty. 
% Hydrological Processes 29: 1873–1882. DOI: 10.1002/hyp.10419
% Ocio Moreno, D., Bulygina, N., Westerberg, I. K., Pappenberger, F., Buytaert, W.  2017. 
% The Role of Rating Curve Uncertainty in Real-Time Flood Forecasting. 
% Water Resources Research, 53, 4197-4213, doi:10.1002/2016WR020225.
%
%
% Input: String with the name of the station. A .mat file must exists with
% the name of that station containing the following variables
%   data: N1x2 matrix with the stage-discharge direct gaugings
%   ORC: N2x2 with the official rating curve
%   Qmean: average observed discharge
% Output:
%   RCU: matrix with the RCU. A 1001 CDF of the true discharge is obtained
%   at each stage from 0 to maximum of ORC with a precision of 1cm.
%   qsX: array with the X percentile of the true discharge
% To run the code with the example data, use command VPlogistic_3seg('VPMratingdata');
      
clc
load(station,'-mat')
% Prior distribution of parameters a1 b1 b2 b3 br1 br2 a2 a3 - specify prior bounds below
% The specification of the priors is important, particulary in sections
% where there is a large extrapolation of the rating curve. Priors can be
% set using standard values (McMillan and Westerberg, 2015) or based on
% hydraulic knowledge (Ocio et al., 2017).
Lbound=[3.65 1.57 1.0 1.05 1.17 3.24 50.4 42.4];
Ubound=[4.46 3.53 2.16 2.36 1.43 3.96 113.3 95.4];

%Calculation of Initial covariance
for i=1:10000
    theta(i,:)= Lbound(1:6)+ (Ubound(1:6)-Lbound(1:6)).*rand(1,6);
end
C=cov(theta);
N=1000000; %Number of samples. Decrease for quickness
theta=zeros(N,8);
p=zeros(N,1);
%Initial parameters taken as centre of priors
thetaOld(1)=(exp(Lbound(1))+exp(Ubound(1)))/2;
thetaOld(2:8)=(Lbound(2:8)+Ubound(2:8))./2;
naccept = 0; %Number of accepted samples
pOld = L(thetaOld,data,Qmean);
disp('Sampling with MCMC')
cont=1;
for t=1:N
    if t>N*cont/20
        fprintf('Progress: %i%%\n',cont*5)
        cont=cont+1;
    end
    %Recalculate covariance during burn-in period
    if t<N/3&fix(t/10000)==t/10000
        C=cov(theta(t-9999:t,1:6));
    end
    thetaOld(1)=log(thetaOld(1));
    out=1;
    while out>0
        %Sampling from jump distribution
        thetaNew(1:6)= mvnrnd(thetaOld(1:6),C,1);
        thetaNew(7)=(exp(thetaNew(1))*thetaNew(5)^thetaNew(2))/(thetaNew(5)^thetaNew(3)); %Estimation of a2 by continuity
        thetaNew(8)=(thetaNew(7)*thetaNew(6)^thetaNew(3))/(thetaNew(6)^thetaNew(4)); %Estimation of a3 by continuity
        out=0;
        %Check if new sample is within prior bounds
        for k=1:8
            if thetaNew(k)>Ubound(k)|thetaNew(k)<Lbound(k)
                out=1;
            end
        end
    end
    thetaOld(1)=exp(thetaOld(1));
    thetaNew(1)=exp(thetaNew(1));
    thetaNew(7)=(thetaNew(1)*thetaNew(5)^thetaNew(2))/(thetaNew(5)^thetaNew(3)); %Estimation of a2 by continuity
    thetaNew(8)=(thetaNew(7)*thetaNew(6)^thetaNew(3))/(thetaNew(6)^thetaNew(4)); %Estimation of a3 by continuity
    pNew=L(thetaNew,data,Qmean);
    %MCMC move
    r=min(1,pNew/pOld);
    if pOld==0&pNew==0
        r=0;
    end
    u = rand(1,1);
    if u < r
        theta(t,:) = thetaNew;
        p(t)=pNew;
        naccept = naccept + 1;
        pOld = pNew;
        thetaOld=thetaNew;
    else
        theta(t,:) = thetaOld;
        p(t)=pOld;
    end
end
disp('Progress: 100%')
disp(horzcat(['Ratio of aceptance= ',num2str(round(naccept/N,3))]))
%Plot the histogram of posterior distributions
figure(1);
subplot(3,2,1)
edges = [exp(Lbound(1)):0.05*(exp(Ubound(1))-exp(Lbound(1))):exp(Ubound(1))];
histogram(theta(:,1),edges,'Normalization','probability');
title('Parameter a1')
subplot(3,2,2)
edges=[Lbound(2):0.05*(Ubound(2)-Lbound(2)):Ubound(2)];
histogram(theta(:,2),edges,'Normalization','probability');
title('Parameter b1')
subplot(3,2,3)
edges=[Lbound(3):0.05*(Ubound(3)-Lbound(3)):Ubound(3)];
histogram(theta(:,3),edges,'Normalization','probability');
title('Parameter b2')
subplot(3,2,4)
edges=[Lbound(4):0.05*(Ubound(4)-Lbound(4)):Ubound(4)];
histogram(theta(:,4),edges,'Normalization','probability');
title('Parameter b3')
subplot(3,2,5)
edges=[Lbound(5):0.05*(Ubound(5)-Lbound(5)):Ubound(5)];
histogram(theta(:,5),edges,'Normalization','probability');
title('Parameter br1')
subplot(3,2,6)
edges=[Lbound(6):0.05*(Ubound(6)-Lbound(6)):Ubound(6)];
histogram(theta(:,6),edges,'Normalization','probability');
title('Parameter br2')

% Generate and plot credibility intervals for rating curve
disp('Obtaining and plotting RCU')
theta=theta(floor(N/3):end,:);
p=p(floor(N/3):end,:);
i=1:10:length(theta);
theta=theta(i',:);
p=p(i',:);
N=length(i);
MN=round(max(ORC(:,1)),2);
%Obtaining rating curves for each set of parameters
qs=zeros(N,MN*100+1);
for t=1:N
    l=0:0.01:MN;
    for k=1:(MN*100+1)
        if l(k)<theta(t,5)
            qs(t,k)=theta(t,1)*l(k)^theta(t,2);
        elseif l(k)<theta(t,6)
            qs(t,k)=theta(t,7)*l(k)^theta(t,3);
        else
            qs(t,k)=theta(t,8)*l(k)^theta(t,4);
        end
    end
end
%Obtaining percentiles
qs(:,1:end) = sort(qs(:,1:end),'ascend'); 
data = sortrows(data,1); 
qs005=qs(floor(N*0.005),:);
qs05=qs(floor(N*0.05),:);
qs15=qs(floor(N*0.15),:);
qs25=qs(floor(N*0.25),:);
qs35=qs(floor(N*0.35),:);
qs45=qs(floor(N*0.45),:);
qs55=qs(floor(N*0.55),:);
qs65=qs(floor(N*0.65),:);
qs75=qs(floor(N*0.75),:);
qs85=qs(floor(N*0.85),:);
qs95=qs(floor(N*0.95),:);
qs995=qs(floor(N*0.995),:);
%Obtaining RCU
for i=1:(MN*100+1)
    for j=1:1001
        RCU(j,i)=prctile(qs(:,i),(j-1)/10);
    end
end
%Obtaining aleatory uncertainty
qcut = min (data(:,2)/Qmean, 14.8);
qcut = max (qcut, 0.0782);
sig=(4.18*exp(-3.051*qcut)+3.531)/100;%logistic scale parameter
aunc(1,:)=data(:,2).*(1+sig*log(0.025/(1-0.025)));
aunc(2,:)=data(:,2).*(1-sig*log(0.025/(1-0.025)));
figure(2)
fillyy(l,qs005,qs995,[0.93 0.93 1]);
hold on
fillyy(l,qs05,qs95,[0.871 0.871 1]);
hold on
fillyy(l,qs15,qs85,[0.729 0.729 1]);
hold on      
fillyy(l,qs25,qs75,[0.6 0.6 1]) 
hold on
fillyy(l,qs35,qs65,[0.4 0.4 1])
hold on
fillyy(l,qs45,qs55,[0 0 1])
hold on
plot(data(:,1),data(:,2),'k.','Markersize',9)
hold on
plot(ORC(:,1),ORC(:,2),'r')
hold on
line([data(:,1)';data(:,1)'],aunc,'Color','k')
xlabel('Level (m)')
ylabel('Discharge (m3/s)')
xlim([0 MN])
ylim([0 1200])
legend('99%CI','90%CI','70%CI','50%CI','30%CI','10%CI','Gaugings','Official RC','Aleat. uncer.','location','northwest')
title(horzcat(['RCU at ',station]))
end

%%%%%%%%%%V
function pt = L(theta,data,Qmean)
% Find the likelihood of a rating curve with theta parameters using the
% logistic distribution with scale parameter as in Coxon et al (2015)
ns=size(data,1);
p=zeros(1,3);
fit=zeros(1,3);
num=zeros(1,3);
hfit=zeros(3,ns);
qfit=zeros(3,ns);
h=zeros(3,ns);
q=zeros(3,ns);
for i=1:ns
    % Calculate the logistic scale parameter
    qcut = min (data(i,2)/Qmean, 14.8);
    qcut = max (qcut, 0.0782);
    sig=(4.18*exp(-3.051*qcut)+3.531)/100;
    % Calculate the likelihood for each section of the curve
    if data(i,1)<theta(5)
        num(1)=num(1)+1;
        h(1,num(1))=data(i,1);
        q(1,num(1))=data(i,2);
        qsu(i,1)=theta(1)*(data(i,1)+0.005)^theta(2);
        qsm(i,1)=theta(1)*data(i,1)^theta(2);
        qsl(i,1)=theta(1)*(data(i,1)-0.005)^theta(2);
        Rgu=qsu(i,1)/data(i,2);
        Rgm=qsm(i,1)/data(i,2);
        Rgl=qsl(i,1)/data(i,2);
        if abs(Rgu-1)<abs(sig*log(0.025/(1-0.025)))|abs(Rgl-1)<abs(sig*log(0.025/(1-0.025)))
            fit(1)=fit(1)+1;
            hfit(1,fit(1))=data(i,1);
            qfit(1,fit(1))=data(i,2);
            p(1)=p(1)+4*exp(-(Rgm-1)/sig)/(1+exp(-(Rgm-1)/sig))^2;
        end
    elseif data(i,1)<theta(6)
        num(2)=num(2)+1;
        h(2,num(2))=data(i,1);
        q(2,num(2))=data(i,2);
        qsu(i,1)=theta(7)*(data(i,1)+0.005)^theta(3);
        qsm(i,1)=theta(7)*data(i,1)^theta(3);
        qsl(i,1)=theta(7)*(data(i,1)-0.005)^theta(3);
        Rgu=qsu(i,1)/data(i,2);
        Rgm=qsm(i,1)/data(i,2);
        Rgl=qsl(i,1)/data(i,2);
        if abs(Rgu-1)<abs(sig*log(0.025/(1-0.025)))|abs(Rgl-1)<abs(sig*log(0.025/(1-0.025)))
            fit(2)=fit(2)+1;
            hfit(2,fit(2))=data(i,1);
            qfit(2,fit(2))=data(i,2);
            p(2)=p(2)+4*exp(-(Rgm-1)/sig)/(1+exp(-(Rgm-1)/sig))^2;
        end
    else
        num(3)=num(3)+1;
        h(3,num(3))=data(i,1);
        q(3,num(3))=data(i,2);
        qsu(i,1)=theta(8)*(data(i,1)+0.005)^theta(4);
        qsm(i,1)=theta(8)*data(i,1)^theta(4);
        qsl(i,1)=theta(8)*(data(i,1)-0.005)^theta(4);
        Rgu=qsu(i,1)/data(i,2);
        Rgm=qsm(i,1)/data(i,2);
        Rgl=qsl(i,1)/data(i,2);
        if abs(Rgu-1)<abs(sig*log(0.025/(1-0.025)))|abs(Rgl-1)<abs(sig*log(0.025/(1-0.025)))
            fit(3)=fit(3)+1;
            hfit(3,fit(3))=data(i,1);
            qfit(3,fit(3))=data(i,2);
            p(3)=p(3)+4*exp(-(Rgm-1)/sig)/(1+exp(-(Rgm-1)/sig))^2;
        end
    end
end
% Apply the weights according to the H-Q portion of the curve actually fitted 
for i=1:3
    if fit(i)==0
        p(i)=0;
    else
        hfiti=hfit(i,1:fit(i));
        qfiti=qfit(i,1:fit(i));
        hi=h(i,1:num(i));
        qi=q(i,1:num(i));
        hmin=min(hi);
        hmax=max(hi);
        qmin=min(qi);
        qmax=max(qi);
        hminfit=min(hfiti);
        hmaxfit=max(hfiti);
        qminfit=min(qfiti);
        qmaxfit=max(qfiti);
        if i<3
            p(i)=p(i)*((log(hmaxfit)-log(hminfit))*(log(qmaxfit)-log(qminfit)))/((log(hmax)-log(hmin))*(log(qmax)-log(qmin)));
        else
            p(i)=p(i)*((hmaxfit-hminfit)*(qmaxfit-qminfit))/((hmax-hmin)*(qmax-qmin));
        end
        if hmin==hmax
            p(i)=0;
        end
    end
end
%Divide by number of points in each section 
p(1) = p(1)/sum(data(:,1)<theta(5)); 
p(2) = p(2)/sum(data(:,1)>=theta(5) & data(:,1)<theta(6)); 
p(3) = p(3)/sum(data(:,1)>=theta(6)); 
% Final likelihood product of the 3 sections
pt=p(1)*p(2)*p(3);
end

function out=fillyy(x,y1,y2,col)
% Function to help plotting shaded areas
x  = x(:)';
y1 = y1(:)';
y2 = y2(:)';
n   = length(x);
X = [ x(1),  x,  x(n),  fliplr(x)  ];
Y = [ y1(1), y2, y1(n), fliplr(y1) ];
h=fill(X,Y,col,'Linestyle','none');
if nargout>0
  out=h;
end
end
