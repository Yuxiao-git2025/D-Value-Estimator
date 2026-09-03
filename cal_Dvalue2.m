% =========================================================================
% Calculate the correlation integral C(r) and estimate the fractal
% dimension using automatically selected log-log scaling range
% Copyright: XA.Yu (2026/08)
% Two interval selection methods are used:
%   Method-1: Sliding-Scoring-Determination method
%   Method-2: Slope-Drop-Determination method
% =========================================================================
% Inputs:
%   loc : [longitude, latitude]
%   r1  : Minimum scale, unit: meter
%   r2  : Maximum scale, unit: meter
%   Num : Number of logarithmically spaced scale points
%   plt : logical flag for plotting
%
% Outputs:
%   D       : estimated fractal dimension
%   Dstd    : standard error of D
%   rseq    : logarithmically spaced scale sequence
%   Cr      : correlation integral C(r)
%   Fit     : fitted log10(C) values over the selected fitting interval
%   Dist    : pairwise distances
%   Methods : selected method flag
%             1 = Method-1
%             2 = Method-2
%             0 = empirical 70% fitting interval
% function [D]=cal_Dvalue2(loc, r1, r2, Num, plt)
% =========================================================================
function [D,Dstd,Methods,finalR2,rseq,Cr,Fit,rij]=cal_Dvalue2(loc,r1,r2,Num,plt)


% ===============
diagnosis=false;  % Control whether print (**default**)
% diagnosis=true;  % Control whether print (**debug**)
% ===============

isEmpirical1 = false;   % Whether Method-1 uses empirical fallback
isEmpirical2 = false;   % Whether Method-2 uses empirical fallback
if nargin < 4 || isempty(Num); Num=30; end % 30~50 is Okay
if Num<30; warning('*****  Try more divisions  *****'); return; end
if nargin < 5 || isempty(plt); plt = false; end
if size(loc, 2) < 2;  error('# LOCATION IS [LON LAT (DEPTH)]'); end
if size(loc, 1) < 15; error('# NUMS ERROR'); end

%% Calculate pairwise distances
clear Neve Nij count rij rijm
Neve=size(loc,1);
Nij=Neve*(Neve-1)/2; % i.e., C(Nloc|2)
rij=zeros(Nij,1);
count=0;
% The total pairs of distance is calculated as: C(N|2), e.g, C(50|2)=1225
for i=1:Neve-1
    % Corresponding package is needed for function <DISTANCE>
    [rijm]=distance(loc(i+1:end, 2),loc(i+1:end, 1),loc(i, 2),loc(i, 1), ...
        wgs84Ellipsoid('meter')); % row1=longitude  row2=latitude
    rij(count+1:count+length(rijm))=rijm(:);
    count=count+length(rijm);
end
rij=rij(1:count);
% Zeros are not useful for logarithmic scale analysis
rij=rij(isfinite(rij) & rij>0);
if isempty(rij)
    fprintf('# Please Check Inputs \n');
    return;
end

%% This can be auto-determined itself but better give the specific value
if isempty(r2)
    r2=max(rij)*0.9;
else
    r2=min(r2, max(rij)*0.9);
end
if isempty(r1)
    r1=10; % 10 m is very small
else
    r1=max(r1, 10);
end
if r1 <= 1e-12; error('# NOTICE: R1>0'); end
if r2 <= r1; error('# NOTICE: R2>R1'); end

%% Calculate the correlation integral C(r)
rseq=logspace(log10(r1), log10(r2), Num); % Feature distance
Cr=zeros(size(rseq));
for k=1:length(rseq)
    rseqi=rseq(k);
    % Count unordered pairs with distance smaller than r 
    % (modeled as a Heaviside function H[r-rij])
    Nr=sum(rseqi>rij);
    % Normalize by the total number of unordered point pairs
    Cr(k)=2*Nr/(Neve*(Neve-1));
end


% C(r)=0 cannot be transformed using log10
% C(r)=1 usually represents the saturation region and is excluded
validC=isfinite(rseq) & isfinite(Cr) & Cr > 0 & Cr < 1; % logical index
Index=find(validC); % position index
if sum(validC)<10
    error('# Too few valid points are available for fitting');
end
% log10 transform now
x=log10(rseq(validC)); x = x(:);
y=log10(Cr(validC));   y = y(:);
Nva=length(x);



%% Calculate local slope and curvature
% Testing revealed that smoothing both the original and derivative data 
% yields more stable results
% Although this results in some loss of information, overall it is a good 
% approach to finding the optimal break-point
Nsm=3;          % Nums of smoothing-window 
Nsm=max(5,Nsm); % typically 3 or 5 for Num=(20~40);
ys=smoothdata(y,'movmean',Nsm);         % Smooth-(*1)
G1=gradient(ys, x);         % Scale-dependent local fractal dimension
G1s=smoothdata(G1, 'movmean', Nsm);     % Smooth-(*2)
G2=gradient(G1s,x);         % i.e., curvature of the log-log curve
G2s=smoothdata(G2, 'movmean', Nsm);   % Smooth-(*3)
% Minimum fitting window
fwin0=max(5,ceil(0.3*Nva)); % i.e., 30% percent

%% ========================================================================
%                             Methods-1: SSD
% ========================================================================
Ends=fwin0:Nva;
Ni=length(Ends);
R2_1        = NaN(Ni, 1);
Dstd_1      = NaN(Ni, 1);
RMSE_1      = NaN(Ni, 1);
SlopeVar_1  = NaN(Ni, 1);
SlopeDri_1  = NaN(Ni, 1);
CurvMean_1  = NaN(Ni, 1);
LenReward_1 = NaN(Ni, 1);

for ii=1:Ni
    Endii=Ends(ii);
    % Now fitting using the original data (i.e., logC(r) versus logr)
    x0=x(1:Endii);
    y0=y(1:Endii);
    [p0,s0]=polyfit(x0, y0, 1);
    y0p=polyval(p0, x0);
    res=y0-y0p;
    SSres=sum(res.^2);
    SStot=sum((y0-mean(y0)).^2);
    % Fitting quality and Std errors
    R2_1(ii)=1-SSres/SStot;
    RMSE_1(ii)=sqrt(mean(res.^2));
    Dstd_1(ii)=Fun_SlopeStd(s0);
    % Slope stability
    slope0=G1s(1:Endii);
    slopeRef=median(slope0,'omitnan');
    SlopeVar_1(ii)=std(slope0,'omitnan')/(abs(slopeRef));
    % Slope drift between the initial part and the tail part
    LL1=max(3,ceil(0.15*Endii));
    LL2=max(3,ceil(0.15*Endii));
    slope1=median(G1s(1:LL1),'omitnan');
    slope2=median(G1s(Endii-LL2+1:Endii),'omitnan');
    SlopeDri_1(ii)=abs(slope2-slope1)/(abs(slope1));
    % Mean absolute curvature
    CurvMean_1(ii)=mean(abs(G2s(1:Endii)),'omitnan')/(abs(slopeRef))*(x0(end)-x0(1));
    % Longer fitting interval is slightly rewarded
    LenReward_1(ii)=sqrt(Endii/Nva);
end
% Normalize penalty terms [0~1]
DstdN_1     = Fun_Norm01(Dstd_1);
RMSEN_1     = Fun_Norm01(RMSE_1);
SlopeVarN_1 = Fun_Norm01(SlopeVar_1);
SlopeDriN_1 = Fun_Norm01(SlopeDri_1);
CurvMeanN_1 = Fun_Norm01(CurvMean_1);

% Comprehensive score (weight can be adjusted)
% Larger R2 and longer interval, Smaller Dstd, RMSE, slope variation,
% slope drift, and curvature are preferred
Score=1.6 * R2_1 - 0.4 * RMSEN_1...
    - 0.8 * DstdN_1  ...
    - 0.3 * SlopeVarN_1 - 0.3 * SlopeDriN_1 - 0.2 * CurvMeanN_1 ...
    + 1.0 * LenReward_1;

% Avoid invalid scores
Score(~isfinite(Score))=-Inf;
[~,MaxScore]=max(Score);
End1=Ends(MaxScore);

% Empirical fallback
Empis=0.7; % 0.5-0.7 is Okay
EndEmpis=ceil(Empis*Nva);
if End1<=fwin0
    End1=EndEmpis;
    Method='Empirical-1';
    isEmpirical1=true;
end
%                       >> Fit using Method-1 <<
[D_1, Dstd_Final_1, Fit_1, R2_Final_1, RMSE_Final_1]=Fun_FitLogLog(x,y,(1:End1)');

%% ========================================================================
%                             Method-2: SDD
%                  (Usually consider it a complementary)
% ========================================================================
% What we do:
%   1. Find the first significant local maximum of D0s
%   2. If D0s drops more than a given ratio after this peak,
%      regard the peak as the turning point
%   3. Fit the prefix before this turning point
% What you need:
%   Adjust the parameter of dropRatio and dropCut
% Where:
%   peakRatio(<100):
%       Controls how significant a local peak should be compared with the
%       initial slope fluctuation
%   dropRatio(<100):
%       Controls how much the slope should drop after the peak before
%       declaring a turning behavior
%   dropCut(<100):
%       Used to determine the position of the true turning point
peakRatio=50; % A conservative estimation (*Not very useful*)
dropRatio=35;  % set 25%-50% is Okay (*need adjustment*)
dropCut=80;    % set 60%~80% is Okay (*need adjustment*)

% Search the first meaningful local peak
% 40% avoids selecting a later secondary peak
Endss=max(fwin0,ceil(0.4*Nva));
Endss=min(Endss,Nva-1);
MAXID=[];
if Endss>=3
    for ii=2:Endss
        if G1s(ii)>=G1s(ii-1) && G1s(ii)>G1s(ii+1)
            MAXID=[MAXID; ii]; %#ok<AGROW>
        end
    end
end

% Simple Significance Test for local peak
% This may not be entirely accurate, but generally the initial slope
% will gradually increase in the previous part
initWin=max(4, ceil(0.1 * Nva));
initWin=min(initWin, Nva);
initMean=mean(G1s(1:initWin));
initStd=std(G1s(1:initWin));
sigPeakIdx = [];
for ii=1:length(MAXID)
    idx=MAXID(ii);
    if G1s(idx)>=initMean+peakRatio*1e-2*initStd
        sigPeakIdx=idx;
        break;
    end
end
% If no significant local peak is found, use the maximum in the early half
if isempty(sigPeakIdx)
    earlyEnd=max(fwin0,ceil(0.5*Nva));
    [~,tmpPeakIdx]=max(G1s(1:earlyEnd));
    peakIndex=tmpPeakIdx;
else
    peakIndex=sigPeakIdx;
end
peakValue=G1s(peakIndex);
% Check whether D0s drops by more than dropRatio after the peak
% For positive slopes:
%                   CCrit=peakValue*(1-dropValue)
dropValue=peakValue*(1-dropRatio*1e-2);
dropIdx=find(G1s(peakIndex+1:end)<=dropValue,1,'first');
if ~isempty(dropIdx)
    % Select the turning point considering both peaking and flating
    End21=peakIndex;
    End22=peakIndex+dropIdx;
%     if peakValue>=1 && dropValue>=1
%         warning(['# You may need check the input paramters of Method-2 and' ...
%             'try some empirical trials']);
%     end
    End2=round(prctile(End21:1:End22,dropCut)); % percentile is used
else
    % Otherwise use empirical interval
    End2=EndEmpis;
    Method='Empirical-2';
    isEmpirical2=true;
end
% Empirical fallback
if End2<=fwin0
    End2=EndEmpis;
    Method='Empirical-2';
    isEmpirical2=true;
end
%                       >> Fit using Method-2 <<
[D_2, Dstd_Final_2, Fit_2, R2_Final_2, RMSE_Final_2]=Fun_FitLogLog(x,y,(1:End2)');

%% ========================================================================
% Select the final method according to the smaller sigma-value, i.e., Dstd
% ========================================================================
sig1=Dstd_Final_1;
sig2=Dstd_Final_2;
% sig1=R2_Final_1;
% sig2=R2_Final_2;
if ~isfinite(sig1);sig1 = Inf;end
if ~isfinite(sig2);sig2 = Inf;end
if sig1<=sig2
    Method='Method-1'; %#ok<NASGU>
    D=D_1;
    Dstd=Dstd_Final_1;
    Fit=Fit_1;
    idFitLocal=(1:End1)';
    finalR2=R2_Final_1;       
    finalRMSE=RMSE_Final_1;   %#ok<NASGU>
    finalEndIndex=End1;       %#ok<NASGU>
    if isEmpirical1
        Methods=0;
    else
        Methods=1;
    end
    if sig1==sig2
        Method='Both'; %#ok<NASGU>
    end
else
    Method='Method-2'; %#ok<NASGU>
    D=D_2;
    Dstd=Dstd_Final_2;
    Fit=Fit_2;
    idFitLocal=(1:End2)';
    finalR2=R2_Final_2;       
    finalRMSE=RMSE_Final_2;   %#ok<NASGU>
    finalEndIndex=End2;       %#ok<NASGU>
    if isEmpirical2
        Methods=0;
    else
        Methods=2;
    end
end
% Plot the results
NewIndex=Index(idFitLocal);
if plt
    Fit_Plot(rseq,Cr,D,Fit,NewIndex);
end
%% Optional diagnostic plots
if diagnosis
    fprintf('# %s is selected: End=%.0f \n', Method, finalEndIndex);
    % figure;
    nexttile; hold on;
    rseq=rseq(validC);
    plot(rseq(Ends), Score, 'Color',[0 0 0], 'LineWidth', 2.6);
    xline(rseq(End1), '--','Color',[0 0 0],'LineWidth',1.5);
    yline(Score((Ends==End1)), '--', 'Method-1','FontSize',16,'Color',[0 0 0],'LineWidth',1.5);
    plot(rseq(End1),Score((Ends==End1)),'Marker','pentagram','LineStyle','none', ...
        'MarkerSize',18,'MarkerFaceColor','none','MarkerEdgeColor',[0 0 0],'LineWidth',1.5);
    ylabel('Score');
    xlim('tight');
    yl=ylim;
    ylim([yl(1)-0.2 yl(2)+0.2]);
    yyaxis right;hold on;
    plot(rseq,G1s, 'k', 'LineWidth', 1.6,'Color',[0 .1 .9]);
    xline(rseq(End2), '-.','Color',[0 .1 .9],'LineWidth',1.5);
    yline(G1s(End2), '-.', 'Method-2','FontSize',16,'Color',[0 .1 .9], ...
        'LineWidth',1.5);
    plot(rseq(End2),G1s(End2),'Marker','pentagram','LineStyle','none', ...
        'MarkerSize',18,'MarkerFaceColor','none','MarkerEdgeColor',[0 .1 .9], ...
        'LineWidth',1.5);
    xlabel('Distance');
    ylabel('1st-Derivative');
    ax=gca;
    ax.YAxis(1).Color=[0 0 0];
    ax.YAxis(2).Color=[0 .1 .9];
    ax.XScale="log";
    if exist('Fun_defaultAxes','file')
        Fun_defaultAxes;
    end
    fprintf('# End1=%.0f \tEnd2=%.0f \n',End1,End2);

end
end





% =========================================================================
%                             Subfunctions
% =========================================================================
function Fit_Plot(rseq,C,D,Fit,Index)
loglog(rseq, C, 'LineWidth', 0.6, 'MarkerSize', 16, ...
    'MarkerFaceColor', [0.9412 0.9412 0.9412], ...
    'LineStyle', 'none', 'Marker', 'square', 'MarkerEdgeColor', 'k');
hold on;
loglog(rseq(Index), 10.^Fit, '-', 'LineWidth', 3.5,'Color',[0 .1 .9]);
% xline(rseq(id(end)),'LineWidth',2,'Color',[.4 .4 .4],'LineStyle',':');
plot(rseq(Index(end)),max(C(Index(end)+min(2,length(C)-length(Index))),10.^Fit(Index(end))),'Marker','v','MarkerSize',12,'MarkerFaceColor', ...
    [0 .2 .9],'LineStyle','none','MarkerEdgeColor','k','LineWidth',0.6);
if exist('Fun_defaultAxes','file')
    Fun_defaultAxes;
end
% xl=xlim; xlim([10, xl(2) * 1.0001]);
% yl=ylim; ylim([1e-2, yl(2) + 0.5]);
xlabel('Distance (m)');
ylabel('C(r)');
end


function [D, Dstd, Fit, R2, RMSE] = Fun_FitLogLog(x, y, idFit)
xFit = x(idFit);
yFit = y(idFit);
[p, S] = polyfit(xFit, yFit, 1);
D = p(1);
Fit = polyval(p, xFit);
residual = yFit - Fit;
SSres = sum(residual.^2);
SStot = sum((yFit - mean(yFit)).^2);

if SStot > 0
    R2 = 1 - SSres / SStot;
else
    R2 = NaN;
end
RMSE = sqrt(mean(residual.^2));
Dstd = Fun_SlopeStd(S);
end



function Dstd = Fun_SlopeStd(S)
if S.df > 0
    Rinv = S.R \ eye(size(S.R));
    Covp = (S.normr^2 / S.df) * (Rinv * Rinv.');
    se_p = sqrt(diag(Covp));
    Dstd = se_p(1);
else
    Dstd = NaN;
end
end



function z = Fun_Norm01(v)
z = NaN(size(v));
validv = isfinite(v);
if sum(validv) == 0
    z(:) = 0;
    return;
end
vValid = v(validv);
vMin = min(vValid);
vMax = max(vValid);
if abs(vMax - vMin) < eps
    z(validv) = 0;
else
    z(validv) = (vValid - vMin) / (vMax - vMin);
end
z(~validv) = 1;
end


