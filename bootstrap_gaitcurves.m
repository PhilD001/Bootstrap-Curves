function [c_b, cilo_b, cihi_b,c,cilo,cihi,mu_b,mu,se_b,se,meanmu_b,bias] = bootstrap_gaitcurves(data,nboots,alpha,display)

% bootstap_gaitcurves computes symmetric confidence intervals based on the bootstrap-t procedure
% 
% ARGUMENTS
%
% data          ...  n x m matrix data with n = num of trials/subjects and m = time points. 
%                    Usually 'm' is normalized to 100% for a complete gait cycle.
% nboots        ...  Number of bootstrap samplies. Default 1000.
% alpha         ...  Significance level. Default 0.05.
% display       ...  Choice to display results. Default 'no'. If using demo
%                    mode for user data, some editing might be needed.
%
% RETURNS
%
% mu            ...  Sample estimate of population mean mu.
% cilo          ...  Sample estimate of lower range of confidence interval.
% cihi          ...  Sample estimate of upper range of confidence interval.
% c             ...  Sample estimate of critical value of test statistic.
% se            ...  Sample estimate of standard error.
% mu_b          ...  Bootstrap estimates of muhat.
% cilo_b        ...  Bootstrap estimate of lower range of confidence interval
% cihi_b        ...  Bootstrap estimate of upper range of confidence interval
% c_b           ...  Bootstrap estimate of critical value of test statistic.
% meanmu_b      ...  Mean bootstrap estimate of muhat.
% bias          ...  Bias between sample and bootstrap estmates of population mean.
% se_b          ...  Bootstrap estimate of standard error (standard deviation of the 'b' 
%                    bootstrap population mean estimates).
%
% NOTES
% - a demo version can be tested by running the function without arguments.
%   This loads knee sagittal plane angle data for 27 subjects over the entire
%   gait cycle (0 - 100%), i.e., n = 27 and m = 101.
%
% - This function should produce Bootstrap intervals as described by 
%   Lenhoff et al. 1999. Gait and Posture, 9(1):10-17. and
%   Duhamel et al. 2004. Gait and Posture, 20(2): 204-212
%
% - Uses support from 'vline' and 'hline' for graphing by Brandon Kuczenski
%   available on the MatLab file exchange
%
%
%
% Created June 2012 by
% Philippe C. Dixon
% PhD Candidate
% Department of Engineering Science
% University of Oxford
% philippe.dixon@gmail.com



%--SET DEFAULTS--------------------------------------------------------------------------
%
if nargin==0
    nboots = 1000;
    alpha = 0.05;
    display = 'yes';
       
    file = which('bootstrap_gaitcurves.m');
    path = fileparts(file);
    data = load([path,slash,'boot_test_data.mat'],'-mat');
    data = data.mstk;
    ch = 'Knee Angle';
    
    disp('running demo mode using sample data')
end


if nargin==1
    nboots = 1000;
    alpha = 0.05;
    display = 'no';
end


if nargin==2
    alpha = 0.05;
    display = 'no';
end


if nargin==3
    display= 'no';
end




%--COMPUTE SAMPLE ESTIMATES---------------------------------------------------------------------
%
[n,cols] = size(data);

mu = mean(data);                                      % sample estimate of mean
shat = std(data);                                     % sample estimate of standard deviation
se = shat./sqrt(n);                                   % sample estimate of standard error

df = n-1;                                             % degrees of freedom of sample
c = tinv(1-alpha/2,df);                               % sample estimate of the critical t-value

cilo = mu-c*se;                                       % sample estimate of lower CI
cihi = mu+c*se;                                       % sample estimate of upper CI



%--COMPUTE BOOTSTRAP ESTIMATES------------------------------------------------------------------
%
mu_b = zeros(nboots,cols);
t = zeros(nboots,1);

for b = 1:nboots
    [~,indx] = datasample(data(:,1),n);               % indx of bth curve
    xb1 = data(indx',:);                              % bth bootstrap sample of data
    mu_b(b,:) = mean(xb1);                            % bth bootstrap estimate of the mean
end

se_b = std(mu_b);                                     % bootstrap standard error of the mean

for b = 1:nboots
    t(b) = max( abs ( mu_b(b,:) - mu )  ./ se_b );    % t-statistic
end



%--COMPUTE TEST STATISTIC AT ALPHA LEVEL--------------------------------------------------------
%
c_b = prctile(t,100*(1-alpha/2));



%--COMPUTE CONFIDENCE INTERVALS-----------------------------------------------------------------
%
cilo_b = mu - c_b*se_b;                               % lower confidence interval
cihi_b = mu + c_b*se_b;                               % upper confidence interval



%--COMPUTE OTHER QUANTITIES---------------------------------------------------------------------
%
meanmu_b = mean(mu_b);
bias = mu-meanmu_b;



%--DISPLAY RESULTS------------------------------------------------------------------------------
%
if isin(display,'yes')
    
    x = 1:1:cols;
    figure('name',ch)
    
    subplot(2,2,1)
    plot(x,mu,'r','linewidth',1.5);
    hold on
    plot(x,meanmu_b,'b','linewidth',1.5);
    
    plot(x,cilo_b,'b--','linewidth',1.5)
    plot(x,cihi_b,'b--','linewidth',1.5)
    plot(x,cilo,'r.-','linewidth',1.5)
    plot(x,cihi,'r.-','linewidth',1.5)
    text(x(5),max(mu),['tcrit sample = ' num2str(c)])
    text(x(5),max(mu)-0.2*max(mu),['tcrit bootstrap = ' num2str(c_b)])
    xlim([0,101])
    title('mean and CI estimates (sample and bootstrap)')
    legend('sample','boot','location','best')

    subplot(2,2,2)
    plot(x,bias,'k','linewidth',1.5)
    hline(0,'k')
    title('bias ( sample mu - bootstrap mu )')
    xlim([0,101])
    
    subplot(2,2,3)
    plot(x,se,'r','linewidth',1.5)
    hold on
    plot(x,se_b,'linewidth',1.5)
    title('Standard error estimates (sample and bootstrap)')
    xlim([0,101])
    legend('sample','boot','location','best')

    
    subplot(2,2,4)
    hist(mu_b(:,18),20)
    title('histogram for bootstrap estimates at t = 18 (peak knee flexion) ')
    hold on
    vline(mu(18),'r:')
    vline(meanmu_b(18),'b')
    vline(cilo_b(18),'b:')
    vline(cihi_b(18),'b:')
    vline(cilo(18),'r:')
    vline(cihi(18),'r:')
    
end


 
 
 
%==EMBEDDED FUNCTIONS=============================================================================

function r= isin(str,a)

% used to simplify code

r = ~isempty(strfind(str,a));


function [s,type] = slash

% ensures funcionality on windows and mac

type = computer;

switch type
    
        
    case {'PCWIN64','PCWIN'}
        s= '\';
        
    case {'MAC','MACI','MACI64'}
        s= '/';
        
    otherwise
        s = '\';
end



function hhh=hline(y,in1,in2)

% function h=hline(y, linetype, label)
% 
% Draws a horizontal line on the current axes at the location specified by 'y'.  Optional arguments are
% 'linetype' (default is 'r:') and 'label', which applies a text label to the graph near the line.  The
% label appears in the same color as the line.
%
% The line is held on the current axes, and after plotting the line, the function returns the axes to
% its prior hold state.
%
% The HandleVisibility property of the line object is set to "off", so not only does it not appear on
% legends, but it is not findable by using findobj.  Specifying an output argument causes the function to
% return a handle to the line, so it can be manipulated or deleted.  Also, the HandleVisibility can be 
% overridden by setting the root's ShowHiddenHandles property to on.
%
% h = hline(42,'g','The Answer')
%
% returns a handle to a green horizontal line on the current axes at y=42, and creates a text object on
% the current axes, close to the line, which reads "The Answer".
%
% hline also supports vector inputs to draw multiple lines at once.  For example,
%
% hline([4 8 12],{'g','r','b'},{'l1','lab2','LABELC'})
%
% draws three lines with the appropriate labels and colors.
% 
% By Brandon Kuczenski for Kensington Labs.
% brandon_kuczenski@kensingtonlabs.com
% 8 November 2001

if length(y)>1  % vector input
    for I=1:length(y)
        switch nargin
        case 1
            linetype='r:';
            label='';
        case 2
            if ~iscell(in1)
                in1={in1};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            label='';
        case 3
            if ~iscell(in1)
                in1={in1};
            end
            if ~iscell(in2)
                in2={in2};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            if I>length(in2)
                label=in2{end};
            else
                label=in2{I};
            end
        end
        h(I)=hline(y(I),linetype,label);
    end
else
    switch nargin
    case 1
        linetype='r:';
        label='';
    case 2
        linetype=in1;
        label='';
    case 3
        linetype=in1;
        label=in2;
    end

    
    
    
    g=ishold(gca);
    hold on

    x=get(gca,'xlim');
%     h=plot(x,[y y],linetype, 'LineStyle', ':');
    h=plot(x,[y y],linetype);
    if ~isempty(label)
        yy=get(gca,'ylim');
        yrange=yy(2)-yy(1);
        yunit=(y-yy(1))/yrange;
        if yunit<0.2
%             text(x(1)+0.02*(x(2)-x(1)),y+0.02*yrange,label,'color',get(h,'color'))
              text(x(1)+0.7*(x(2)-x(1)),y+0.00*yrange,label,'color',get(h,'color'), 'FontSize', 8, 'FontWeight', 'bold')
        else
%             text(x(1)+0.02*(x(2)-x(1)),y-0.02*yrange,label,'color',get(h,'color'))
              text(x(1)+0.7*(x(2)-x(1)),y-0.00*yrange,label,'color',get(h,'color'), 'FontSize', 8, 'FontWeight', 'bold')
        end
    end

    if g==0
    hold off
    end
    set(h,'tag','hline','handlevisibility','off') % this last part is so that it doesn't show up on legends
end % else

if nargout
    hhh=h;
end


function hhh=vline(x,in1,in2)
% function h=vline(x, linetype, label)
% 
% Draws a vertical line on the current axes at the location specified by 'x'.  Optional arguments are
% 'linetype' (default is 'r:') and 'label', which applies a text label to the graph near the line.  The
% label appears in the same color as the line.
%
% The line is held on the current axes, and after plotting the line, the function returns the axes to
% its prior hold state.
%
% The HandleVisibility property of the line object is set to "off", so not only does it not appear on
% legends, but it is not findable by using findobj.  Specifying an output argument causes the function to
% return a handle to the line, so it can be manipulated or deleted.  Also, the HandleVisibility can be 
% overridden by setting the root's ShowHiddenHandles property to on.
%
% h = vline(42,'g','The Answer')
%
% returns a handle to a green vertical line on the current axes at x=42, and creates a text object on
% the current axes, close to the line, which reads "The Answer".
%
% vline also supports vector inputs to draw multiple lines at once.  For example,
%
% vline([4 8 12],{'g','r','b'},{'l1','lab2','LABELC'})
%
% draws three lines with the appropriate labels and colors.
% 
% By Brandon Kuczenski for Kensington Labs.
% brandon_kuczenski@kensingtonlabs.com
% 8 November 2001

if length(x)>1  % vector input
    for I=1:length(x)
        switch nargin
        case 1
            linetype='r:';
            label='';
        case 2
            if ~iscell(in1)
                in1={in1};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            label='';
        case 3
            if ~iscell(in1)
                in1={in1};
            end
            if ~iscell(in2)
                in2={in2};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            if I>length(in2)
                label=in2{end};
            else
                label=in2{I};
            end
        end
        h(I)=vline(x(I),linetype,label);
    end
else
    switch nargin
    case 1
        linetype='r:';
        label='';
    case 2
        linetype=in1;
        label='';
    case 3
        linetype=in1;
        label=in2;
    end

    
    
    
    g=ishold(gca);
    hold on

    y=get(gca,'ylim');
    h=plot([x x],y,linetype);
    if length(label)
        xx=get(gca,'xlim');
        xrange=xx(2)-xx(1);
        xunit=(x-xx(1))/xrange;
        if xunit<0.8
            text(x+0.001*xrange,y(1)+0.01*(y(2)-y(1)),label,'color',get(h,'color')) %change vales*xrange to place text closer to line
        else
            text(x-.05*xrange,y(1)+0.1*(y(2)-y(1)),label,'color',get(h,'color'))
        end
    end     

    if g==0
    hold off
    end
    set(h,'tag','vline','handlevisibility','off')
end % else

if nargout
    hhh=h;
end

