function [N, Adj, kpercent, histArrayPercShort,histArrayShort,k, pval,Rij,s, possible] = links(calciumT, Threshold, Opts,fig)%ii, mm, phase, figs)
%% Network analysis
% Jennifer Briggs 12/2020
% This program was adapted from an original program by Vira Kravets (Apr
% 2020) who adapted the technique from Stozer 2013

%% THIS PROGRAM CALCULATES COVARIANCE BETWEEN ALL CELLS, ADJACENCY MATRIX, probability of links, and plots links (cov>threshold) 
%% INPUT PARAMETERS:
% calciumT = calcium time course;
% Threshold = Rth between 0 and 1 (usually somewhere from 0.8-0.99)
% Opts: structure with options
    % Reiqured opts: 
    % figs
    % If use print as subplot, also need:
    % Subplotnum
%fig sets the figure number

if ~exist('fig')
    fig = 1;
end
    
figs = Opts.figs;
if figs == 1
    try
        if Opts.printasSubplot == 1 
            mm = Opts.Subplotnum;
        end
    catch
        Opts.printasSubplot = 0;
    end
end


[data numcells]=size(calciumT);                                     


time = [1: size(calciumT,1)];

%% 2. Correlation matrix elements Rij (eq.1, stozer, 2013): correlation for each individual cell with the rest of the cells
% as well as P-values for each Rij, and the critical Rij values.
   
[Rij, pval]=corr(calciumT);
    
    
    if figs
    figure(fig)
    fig = fig+1;
    imagesc(Rij);
    title('Correlation Coefficient by Cell')
    colorbar
    figure
    imagesc(pval);
    title('Statistical Significance')
    colorbar
    end



%% 3. Making a link map
mean_coorRij = mean2(Rij);
disp(mean_coorRij)
%disp(mean2(Rij))
Adj = Rij;
Adj(Adj >= Threshold) = 1;
Adj(Adj < Threshold) = 0;

 Adj = Adj - diag(diag(Adj));    % replacing diagonal elemants with 0s to account for self-correlation
 mean_adj = mean2(Adj);
 X = mean2(Adj);
 Y = [X, 'avg coorelation coeff'];
disp(Y)

% Determine what proportion of connections meet threshold in proportion to
% all possible connections: 

all = ones(numcells);
all = all - diag(diag(all));
possible = mean2(all);

ratio_over_thresh = mean_adj/possible
Z = [ratio_over_thresh, 'ratio over threshold/possible'];
disp(Z)
disp(ratio_over_thresh)


if figs
    figure(fig);
    fig = fig+1;
    imagesc(Adj);
    title('Connection Map')
    colorbar
end


%% 4. Determine number of "links" based on cov threshold
for i=1:numcells
    N (i,1) = nnz(Adj(:,i));  % N is matrix containing # of links for each cell (nnz returns number of nonzero elements)
end

%% 5. Creating a "probability of a cell to have k links" matrix

histArray=zeros(size(N))'; % prealocate
histArray = [histArray 0];
% a forloop to count how many times you encounter in a particular value:
for n=1:length(N)
    histArray(1,N(n)+1)=histArray(1,N(n)+1)+1; % every time you meet the particular value, you add 1 into to corresponding bin
end

histArrayPerc=histArray.*100/sum(histArray); % converting into % via dividing # of cell with specific # of links by total # of links 

m=find(histArray);    % returns all indexes of non-zero elements of the array
maxNonZeroLinks=max(m);   % number of links to consider when normalizing probabilty
%maxNonZeroLinks = 20;
k=1:1:maxNonZeroLinks;            % index of # of links (starting with 0 links)
kpercent=k.*100/(maxNonZeroLinks);  % convert # of links into % of links
histArrayPercShort = histArrayPerc(1:maxNonZeroLinks);   % cropping the hisArray after the last non-0 value

histArrayShort = histArray(1:maxNonZeroLinks);
%histArrayShort = histArray(1:maxNonZeroLinks);
loghist = log(histArrayShort);

xlab = [1:length(histArrayPercShort)];
xlab(isinf(loghist))=[];
loghist(isinf(loghist))=[];

try
[s] = corrcoef(log(xlab),loghist);
s = s(2);
catch
    s = 0
end

if figs
    figure(fig)  % plot a bar graph of the probability of a cell to have k links
    fig = fig+1;
    bar(k,histArrayShort)
    title('Probability to have k links,(%)')
    xlabel('Number of links');
    ylabel('Number of cells')
end

if figs == 1
    fig2 = figure(fig),hold on;
    fig2.Position = [325 32 1132 892.5000];
    if Opts.printasSubplot == 1
        if Opts.multiseed == 1
            if mm == 1
                if Opts.multiseednum == 1
                    subplot(Opts.subplotdim+1), hold on
                end
                elseif mm == 2
                if Opts.multiseednum == 1
                    subplot(Opts.subplotdim+2), hold on
                end
            else
            subplot(Opts.subplotdim + mm)
            end
        else
        subplot(Opts.subplotdim + mm)
        end
    end

    % plot a bar graph of the probability of a cell to have k links
    bar(kpercent,histArrayPercShort)
    title(['Threshold = ' num2str(Threshold)],'FontSize', 15)
    %title(['Islet Size = ' num2str(numcells)])
    xlabel('Percent of links, (%)','FontSize', 15);
    ylabel('Pr(links)','FontSize', 15)
    % sgtitle(['Probability to have k links,(%)  ']) 

%% 6. Applying log function to probability and number of links
% plotting on logarithmic scale
    fig3 = figure(fig+1), hold on;
    fig3.Position = [325 32 1132 892.5000];

    if Opts.printasSubplot == 1
        if Opts.multiseed == 1
            if mm == 1
                if Opts.multiseednum == 1
                    subplot(Opts.subplotdim+1), hold on
                end
            elseif mm == 2
                if Opts.multiseednum == 1
                    subplot(Opts.subplotdim+2), hold on
                end
            else
            subplot(Opts.subplotdim + mm)
            end
        else
        subplot(Opts.subplotdim + mm)
        end
    end
plot(log(kpercent),log(histArrayPercShort),'o');
ll = log(histArrayPercShort);
xx = log(kpercent);
xx(isinf(ll)) = [];
ll(isinf(ll)) = [];
coef1 = polyfit(xx,ll , 1);
y1 = polyval(coef1, xx);
hold on, plot(xx,y1)
%title({'Probability on a log(10) scale' phase 'avg links = ' num2str(mean(k))});
title(['Threshold = ' num2str(Threshold)],'FontSize', 15);
%title(['Islet Size = ' num2str(numcells), 'Average Links = '])
xlabel('Log(number of links, (%))','FontSize', 15);
ylabel('Log(Pr(links))','FontSize', 15);
legend(['Power Law fit: R^2 = ' num2str(s)], 'FontSize', 15)
sgtitle(['Probability on a log(10) scale  '])
end

end
