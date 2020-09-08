%% Comparison of PRAME expression in MPNST v. other

clear all; close all; load 'Data.mat'

mpnst = find(strcmp(diagnosis,'MPNST'));
bpnst = intersect(find(~strcmp(diagnosis,'MPNST')),find(~strcmp(diagnosis,'X18')));

t = [sum(PRAME(mpnst)) length(mpnst)-sum(PRAME(mpnst));sum(PRAME(bpnst)) length(bpnst)-sum(PRAME(bpnst))]
[h,p,stats] = fishertest(t)
sum(PRAME(mpnst))/length(mpnst)
sum(PRAME(bpnst))/length(bpnst)

% 54/82 (65.85%) of MPNSTs are PRAME+
% 0/56 (0%) of others are PRAME+
% p=8.9056e-18

%% Compute sensitivity and specificity with confidence intervals
[phat,pci] = binofit(sum(PRAME(mpnst)), length(mpnst),0.05)
% phat = 0.6585, pci = [0.5455 0.7597]
[phat,pci] = binofit(length(find(PRAME(bpnst)==0)),length(bpnst),0.05)
%phat = 1.00, pci = [0.9362 1.0000]

%% Difference in NF1?
nfYes = intersect(find(NF1==1), mpnst);
nfNo = intersect(find(NF1==0), mpnst);
t = [sum(PRAME(nfYes)) length(nfYes)-sum(PRAME(nfYes)); sum(PRAME(nfNo)) length(nfNo)-sum(PRAME(nfNo))] 
[h,p,stats] = fishertest(t)
sum(PRAME(nfYes))/length(nfYes)
sum(PRAME(nfNo))/length(nfNo)

% 32/48 or 66.67% of MPNSTs arising in NF1 patients are positive for PRAME
% 22/34 or 64.71% of MPNSTs not arising in NF1 patients are positive for
% PRAME
% p = 1.0000

%% Kaplan-Meier PFS in PRAME+ v PRAME- MPNST

clear all; close all; load 'Data.mat'

mpnst = find(strcmp(diagnosis,'MPNST'));
y = PFS(mpnst);
for i = 1:length(mpnst)
    if recurrence(mpnst(i))==1
        cens(i) = 0;
    else
        cens(i) = 1;
    end
end

pos = find(PRAME(mpnst)==1);
neg = find(PRAME(mpnst)==0);
y1 = y(pos); y2 = y(neg);
c1 = cens(pos); c2 = cens(neg);

figure()
[f,x,flow,fup] = ecdf(y1,'censoring',c1,'function','survivor');
ax1 = stairs(x,f,'color','b');
hold on
stairs(x,flow,':b')
stairs(x,fup,':b')
[f,x,flow,fup] = ecdf(y2,'censoring',c2,'function','survivor');
ax2 = stairs(x,f,'color','r');
hold on
stairs(x,flow,':r')
stairs(x,fup,':r')
legend([ax1,ax2],{'PRAME Positive','PRAME Negative'})
xlabel('Progression-Free Survival (Months)')
ylabel('Probability')
legend('boxoff')
ax = gca; ax.Box = 'Off';

[b,logl,H,stats] = coxphfit(PRAME(mpnst),y,'Censoring',cens);
text(ax.XLim(2)*0.15,ax.YLim(2)*0.95,['p=' num2str(round(stats.p,4))])
%p=0.8638

%% Kaplan-Meier MFS in PRAME+ v PRAME- MPNST

clear all; close all; load 'Data.mat'

mpnst = find(strcmp(diagnosis,'MPNST'));
y = MFS(mpnst);
for i = 1:length(mpnst)
    if metastasis(mpnst(i))==1 
        cens(i) = 0;
    else
        cens(i) = 1;
    end
end

pos = find(PRAME(mpnst)==1);
neg = find(PRAME(mpnst)==0);
y1 = y(pos); y2 = y(neg);
c1 = cens(pos); c2 = cens(neg);

figure()
[f,x,flow,fup] = ecdf(y1,'censoring',c1,'function','survivor');
ax1 = stairs(x,f,'color','b');
hold on
stairs(x,flow,':b')
stairs(x,fup,':b')
[f,x,flow,fup] = ecdf(y2,'censoring',c2,'function','survivor');
ax2 = stairs(x,f,'color','r');
hold on
stairs(x,flow,':r')
stairs(x,fup,':r')
legend([ax1,ax2],{'PRAME Positive','PRAME Negative'})
xlabel('Metastasis-Free Survival (Months)')
ylabel('Probability')
legend('boxoff')
ax = gca; ax.Box = 'Off';

[b,logl,H,stats] = coxphfit(PRAME(mpnst),y,'Censoring',cens);
text(ax.XLim(2)*0.15,ax.YLim(2)*0.95,['p=' num2str(round(stats.p,4))])
%p=0.2489

%% Kaplan-Meier OS in PRAME+ v PRAME- MPNST

clear all; close all; load 'Data.mat'

mpnst = find(strcmp(diagnosis,'MPNST'));
y = OS(mpnst);
for i = 1:length(mpnst)
    if strcmp(outcome(mpnst(i)),'DOD') 
        cens(i) = 0;
    else
        cens(i) = 1;
    end
end

pos = find(PRAME(mpnst)==1);
neg = find(PRAME(mpnst)==0);
y1 = y(pos); y2 = y(neg);
c1 = cens(pos); c2 = cens(neg);

figure()
[f,x,flow,fup] = ecdf(y1,'censoring',c1,'function','survivor');
ax1 = stairs(x,f,'color','b');
hold on
stairs(x,flow,':b')
stairs(x,fup,':b')
[f,x,flow,fup] = ecdf(y2,'censoring',c2,'function','survivor');
ax2 = stairs(x,f,'color','r');
hold on
stairs(x,flow,':r')
stairs(x,fup,':r')
legend([ax1,ax2],{'PRAME Positive','PRAME Negative'})
xlabel('Overall Survival (Months)')
ylabel('Probability')
legend('boxoff')
ax = gca; ax.Box = 'Off';

[b,logl,H,stats] = coxphfit(PRAME(mpnst),y,'Censoring',cens);
text(ax.XLim(2)*0.15,ax.YLim(2)*0.95,['p=' num2str(round(stats.p,4))])
%p=0.1827

%% %% Kaplan-Meier PFS in PRAME-high v PRAME-low/negative MPNST

clear all; close all; load 'Data.mat'

mpnst = find(strcmp(diagnosis,'MPNST'));
y = PFS(mpnst);
for i = 1:length(mpnst)
    if recurrence(mpnst(i))==1
        cens(i) = 0;
    else
        cens(i) = 1;
    end
end

pos = find(PRAME_2(mpnst)==2);
neg = find(PRAME_2(mpnst)<=1);
y1 = y(pos); y2 = y(neg);
c1 = cens(pos); c2 = cens(neg);

figure()
[f,x,flow,fup] = ecdf(y1,'censoring',c1,'function','survivor');
ax1 = stairs(x,f,'color','b');
hold on
stairs(x,flow,':b')
stairs(x,fup,':b')
[f,x,flow,fup] = ecdf(y2,'censoring',c2,'function','survivor');
ax2 = stairs(x,f,'color','r');
hold on
stairs(x,flow,':r')
stairs(x,fup,':r')
legend([ax1,ax2],{'PRAME High','PRAME Low/Negative'})
xlabel('Progression-Free Survival (Months)')
ylabel('Probability')
legend('boxoff')
ax = gca; ax.Box = 'Off';

[b,logl,H,stats] = coxphfit(PRAME_2(mpnst),y,'Censoring',cens);
text(ax.XLim(2)*0.15,ax.YLim(2)*0.95,['p=' num2str(round(stats.p,4))])
%p=0.6188

%% Kaplan-Meier MFS in PRAME-high v PRAME-low/negative MPNST

clear all; close all; load 'Data.mat'

mpnst = find(strcmp(diagnosis,'MPNST'));
y = MFS(mpnst);
for i = 1:length(mpnst)
    if metastasis(mpnst(i))==1 
        cens(i) = 0;
    else
        cens(i) = 1;
    end
end

pos = find(PRAME_2(mpnst)==2);
neg = find(PRAME_2(mpnst)<=1);
y1 = y(pos); y2 = y(neg);
c1 = cens(pos); c2 = cens(neg);

figure()
[f,x,flow,fup] = ecdf(y1,'censoring',c1,'function','survivor');
ax1 = stairs(x,f,'color','b');
hold on
stairs(x,flow,':b')
stairs(x,fup,':b')
[f,x,flow,fup] = ecdf(y2,'censoring',c2,'function','survivor');
ax2 = stairs(x,f,'color','r');
hold on
stairs(x,flow,':r')
stairs(x,fup,':r')
legend([ax1,ax2],{'PRAME High','PRAME Low/Negative'})
xlabel('Metastasis-Free Survival (Months)')
ylabel('Probability')
legend('boxoff')
ax = gca; ax.Box = 'Off';

[b,logl,H,stats] = coxphfit(PRAME_2(mpnst),y,'Censoring',cens);
text(ax.XLim(2)*0.15,ax.YLim(2)*0.95,['p=' num2str(round(stats.p,4))])
%p=0.3626

%% Kaplan-Meier OS in PRAME-high v PRAME-low/negative MPNST

clear all; close all; load 'Data.mat'

mpnst = find(strcmp(diagnosis,'MPNST'));
y = OS(mpnst);
for i = 1:length(mpnst)
    if strcmp(outcome(mpnst(i)),'DOD') 
        cens(i) = 0;
    else
        cens(i) = 1;
    end
end

pos = find(PRAME_2(mpnst)==2);
neg = find(PRAME_2(mpnst)<=1);
y1 = y(pos); y2 = y(neg);
c1 = cens(pos); c2 = cens(neg);

figure()
[f,x,flow,fup] = ecdf(y1,'censoring',c1,'function','survivor');
ax1 = stairs(x,f,'color','b');
hold on
stairs(x,flow,':b')
stairs(x,fup,':b')
[f,x,flow,fup] = ecdf(y2,'censoring',c2,'function','survivor');
ax2 = stairs(x,f,'color','r');
hold on
stairs(x,flow,':r')
stairs(x,fup,':r')
legend([ax1,ax2],{'PRAME High','PRAME Low/Negative'})
xlabel('Overall Survival (Months)')
ylabel('Probability')
legend('boxoff')
ax = gca; ax.Box = 'Off';

[b,logl,H,stats] = coxphfit(PRAME_2(mpnst),y,'Censoring',cens);
text(ax.XLim(2)*0.15,ax.YLim(2)*0.95,['p=' num2str(round(stats.p,4))])
%p=0.3551