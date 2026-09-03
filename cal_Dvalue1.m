% =========================================================================
% Calculate the correlation integral C(r) and estimate the correlation
% fractal dimension according to Roy and Gupta (2015)
%
% Copyright: XA.Yu (2026/08)
%
% The fitting interval is selected according to the criteria described by
% Roy and Gupta (2015):
%   1. The fitting interval contains at least 51% continuous valid points
%   2. The preferred regression coefficient is R2>=0.99
%   3. R2>=0.98 is acceptable for active seismic regions
%   4. The first interval satisfying the threshold is selected
%   5. If neither threshold is satisfied, the maximum-R2 interval is used
% =========================================================================
% Inputs:
%   loc : [longitude,latitude] or [longitude,latitude,depth]
%   r1  : Minimum scale,unit:meter
%   r2  : Maximum scale,unit:meter
%   Num : Number of logarithmically spaced scale points
%   plt : Logical flag for plotting
%
% Outputs:
%   D       : Estimated correlation fractal dimension
%   Dstd    : Standard error of the fitted slope
%   Methods : Selected criterion flag
%             1=R2>=0.99
%             2=R2>=0.98
%             0=Maximum-R2 empirical fallback
%   rseq    : Logarithmically spaced scale sequence
%   Cr      : Correlation integral C(r)
%   Fit     : Fitted log10(C) values over the selected interval
%   rij     : Pairwise epicentral distances
% =========================================================================
function [D,Dstd,Methods,FinalR2,rseq,Cr,Fit,rij]=cal_Dvalue1(loc,r1,r2,Num,plt)

% ===============
diagnosis=false;  % Control whether print diagnostic results (**default**)
% diagnosis=true;  % Control whether print diagnostic results (**debug**)
% ===============

if nargin<4 || isempty(Num);Num=30;end
if Num<30;warning('***** Try some divisions *****');return;end
if nargin<5 || isempty(plt);plt=false;end
if size(loc,2)<2;error('# LOCATION IS [LON LAT (DEPTH)]');end
if size(loc,1)<15;error('# NUMS ERROR');end

%% Calculate pairwise epicentral distances
clear Neve Nij count rij rijm
Neve=size(loc,1);
Nij=Neve*(Neve-1)/2;
rij=zeros(Nij,1);
count=0;
% The total number of unique pairs is C(N,2)
% For example,C(50,2)=1225
for i=1:Neve-1
    rijm=distance(loc(i+1:end,2),loc(i+1:end,1),...
        loc(i,2),loc(i,1),wgs84Ellipsoid('meter'));
    rij(count+1:count+length(rijm))=rijm(:);
    count=count+length(rijm);
end
rij=rij(1:count);
% Zero and invalid distances cannot be used in logarithmic analysis
rij=rij(isfinite(rij) & rij>0);
if isempty(rij)
    error('# Please Check Inputs');
end
% Determine the analyzed scale range
if isempty(r2)
    r2=max(rij)*0.9;
else
    r2=min(r2,max(rij)*0.9);
end
if isempty(r1)
    r1=10;
else
    r1=max(r1,10);
end
if r1<=1e-12;error('# NOTICE: R1>0');end
if r2<=r1;error('# NOTICE: R2>R1');end

%% Calculate the correlation integral C(r)
rseq=logspace(log10(r1),log10(r2),Num);
Cr=zeros(size(rseq));
for k=1:length(rseq)
    % Heaviside counting function H(r-rij)
    Nr=sum(rij<rseq(k));
    % Normalize by the total number of unique event pairs
    Cr(k)=2*Nr/(Neve*(Neve-1));
end

%% Select valid correlation-integral points
% C(r)=0 is undefined after logarithmic transformation
% C(r)=1 represents the large-scale saturation region
validC=isfinite(rseq) & isfinite(Cr) & Cr>0 & Cr<1;
Index=find(validC);
if sum(validC)<10
    error('# Too few valid points are available for fitting');
end
x=log10(rseq(validC));
y=log10(Cr(validC));
x=x(:);
y=y(:);
Nva=length(x);

%% ========================================================================
%             Roy and Gupta (2015) Linear-Range Determination
% ========================================================================
% At least 51% continuous points are required
MinRatio=0.51;
fwin0=max(3,ceil(MinRatio*Nva));
% Preferred and acceptable R2 thresholds
R99=0.99;
R98=0.98;
% Prefix fitting intervals
Ends=fwin0:Nva;
Ni=length(Ends);
R2Seq=NaN(Ni,1);
DSeq=NaN(Ni,1);
DstdSeq=NaN(Ni,1);
RMSESeq=NaN(Ni,1);

for ii=1:Ni
    Endii=Ends(ii);
    idFit=(1:Endii)';
    [DSeq(ii),DstdSeq(ii),~,R2Seq(ii),RMSESeq(ii)]=...
        Roy_FitLogLog(x,y,idFit);
end

%% Select the first interval satisfying R2>=0.99
SelectID=find(isfinite(R2Seq) & R2Seq>=R99,1,'first');
if ~isempty(SelectID)
    Methods=1;
    Method='Roy-R2-0.99';
else
    % Select the first interval satisfying R2>=0.98
    SelectID=find(isfinite(R2Seq) & R2Seq>=R98,1,'first');
    if ~isempty(SelectID)
        Methods=2;
        Method='Roy-R2-0.98';
    else
        % Empirical fallback using the maximum R2
        ValidID=find(isfinite(R2Seq));
        if isempty(ValidID)
            error('# No valid fitting interval is available');
        end
        [~,LocalMaxID]=max(R2Seq(ValidID));
        SelectID=ValidID(LocalMaxID);
        Methods=0;
        Method='Roy-Maximum-R2';
    end
end
% Final fitting
EndFit=Ends(SelectID);
idFitLocal=(1:EndFit)';
[D,Dstd,Fit,FinalR2,FinalRMSE]=Roy_FitLogLog(x,y,idFitLocal);
NewIndex=Index(idFitLocal);
% Plots
if plt
    Roy_FitPlot(rseq,Cr,D,Fit,NewIndex);
end

%% Optional diagnostic results
if diagnosis
    fprintf('# %s is selected \n',Method);
    nexttile;
    hold on;
    rValid=rseq(validC);
    plot(rValid(Ends),R2Seq,'Color',[.2 .2 .2],'LineWidth',2.2);
%     yline(R99,'--','R^2=0.99',...
%         'Color',[0.8 0 0],...
%         'LineWidth',1.5,...
%         'FontSize',16);
%     yline(R98,'-.','R^2=0.98',...
%         'Color',[0 0.2 0.8],...
%         'LineWidth',1.5,...
%         'FontSize',16);
    xline(rValid(EndFit),'--',Method,'Color',[0 0 0],...
        'LineWidth',1.8,'FontSize',16);
    yline(R2Seq(SelectID),'--','Color',[0 0 0],'LineWidth',1.8,'FontSize',16);
    plot(rValid(EndFit),R2Seq(SelectID),...
        'Marker','pentagram',...
        'LineStyle','none',...
        'MarkerSize',18,...
        'MarkerFaceColor','none',...
        'MarkerEdgeColor','k',...
        'LineWidth',1.8);
    set(gca,'XScale','log');
    xlim([r1 +inf]);
    xlabel('Distance');
    ylabel('R2');
    if exist('Fun_defaultAxes','file')
        Fun_defaultAxes;
    end
end

end

% =========================================================================
%                             Subfunctions
% =========================================================================
function [D,Dstd,Fit,R2,RMSE]=Roy_FitLogLog(x,y,idFit)
xFit=x(idFit);
yFit=y(idFit);
[p,S]=polyfit(xFit,yFit,1);
D=p(1);
Fit=polyval(p,xFit);
residual=yFit-Fit;
SSres=sum(residual.^2);
SStot=sum((yFit-mean(yFit)).^2);
if SStot>0
    R2=1-SSres/SStot;
else
    R2=NaN;
end
RMSE=sqrt(mean(residual.^2));
Dstd=Roy_SlopeStd(S);
end



function Dstd=Roy_SlopeStd(S)
if S.df>0
    Rinv=S.R\eye(size(S.R));
    Covp=(S.normr^2/S.df)*(Rinv*Rinv.');
    se_p=sqrt(diag(Covp));
    Dstd=se_p(1);
else
    Dstd=NaN;
end
end



function Roy_FitPlot(rseq,C,D,Fit,Index)
loglog(rseq, C, 'LineWidth', 0.6, 'MarkerSize', 16, ...
    'MarkerFaceColor', [0.9412 0.9412 0.9412], ...
    'LineStyle', 'none', 'Marker', 'square', 'MarkerEdgeColor', 'k');
hold on;
loglog(rseq(Index), 10.^Fit, '-', 'LineWidth', 3.5,'Color',[0 .1 .9]);
% xline(rseq(id(end)),'LineWidth',2,'Color',[.4 .4 .4],'LineStyle',':');
plot(rseq(Index(end)),max(C(Index(end)+2),10.^Fit(Index(end))),'Marker','v','MarkerSize',12,'MarkerFaceColor', ...
    [0 .2 .9],'LineStyle','none','MarkerEdgeColor','k','LineWidth',0.6);
if exist('Fun_defaultAxes','file')
    Fun_defaultAxes;
end
xl=xlim; xlim([10, xl(2) * 1.0001]);
yl=ylim; ylim([1e-2, yl(2) + 0.5]);
xlabel('Distance (m)');
ylabel('C(r)');
end
