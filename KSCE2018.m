%{
Developed by Byun, Ji-Eun
Hypothetical examples for KSCE 2018 paper

Output:
for each example (1-5),
1. PI_proposed: Probability of Indep. using reference prior & quantile-based discretization (proposed)
2. PI_uniform: Probability of Indep. using reference prior & uniform discretization
3. Pvalk: P-value of KCI test
=> All are 1000 x 6 ( 1000 exp's done for each data size  x 6 data size being examined)
%}

clear; close all;
%%
N_ = [25,50,100,200,400,1000];
nstart = 4;
import KCI.*

PI_proposed = zeros(1e3,length(N_));
PI_uniform = PI_proposed;
Pvalk = PI_proposed;
for kk = 1:5
    for jj = 1:length(N_)
        N = N_(jj);
        clearvars -except N_ N nstart jj kk PI_proposed PI_uniform Pvalk

        for ii = 1:length(PI_proposed)
            rng(ii)

            x = normrnd(0,1,N,1); y = normrnd(0,1,N,1);
            switch kk
                case 2
                    y = gamrnd( 2,2,[N,1] ); 
                case 3
                    y = sqrt(abs(x)*4) + 0.5*y;
                case 4
                    y = sin(pi*x) + sqrt(0.5)*y;
                case 5
                    y = 2*(1+2*pi*abs(x)) .* sin(2*pi*abs(x))+y;
                    x = 2*(1+2*pi*abs(x)) .* cos(2*pi*abs(x))+x;
            end

%             Bayesian Independence Test (proposed: reference prior + quantil-based discretization)
            PI_ = BUI( x,y );
            PI_proposed(ii,jj) = min(PI_); % if pval < 0.5: dependent; > independent
            
%             Bayesian Independence Test (proposed: reference prior + uniform discretization)
            PI_ = BUI_uni( x,y );
            PI_uniform(ii,jj) = min(PI_); % if pval < 0.5: dependent; > independent

%             KCI test (Zhang et el., 2012)
            [pval,stat] = UInd_KCItest( x,y );
            Pvalk(ii,jj) = pval; % if pval < threshold: dependent; > independent

            if mod(ii,100) == 1; sprintf( '%d th experiment',ii ); end

        end
    end
    save( sprintf('case%d',kk) )
end

%% Figure
fname = {'case1','case2','case3','case4','case5'};
sidx = [1,2,4,5,6];
titles = {'Ex. 1','Ex. 2','Ex. 3','Ex. 4','Ex. 5'};

figure;
for jjj = 1:length(fname)

    load(fname{jjj})
    errr = zeros(2,length(N_)); meanr = zeros(1,length(N_));
    errk = errr; meank = meanr;
    erru = errr; meanu = meanr;

    for kkk = 1:length(N_)

        meanr(kkk) = mean( PI_proposed(:,kkk) );
        errr(:,kkk) = abs( quantile( PI_proposed(:,kkk),[0.25;0.75] ) - meanr(kkk) );

    end
    subplot(2,3,jjj)
    errorbar( 1:length(N_),meanr,errr(1,:),errr(2,:) ,'s--','LineWidth',1.2 )
    axis([0,length(N_)+1,0,1])
    grid on
    set( gca,'XTick',1:length(N_) )
    set( gca,'XTickLabel',{'25','50','100','200','400','1,000'} )
    xlabel( 'No. of samples' );
    ylabel( 'p(H_0|D) (25%-q, mean, 75%-q)' );
    hold on
    plot( [0,length(N_)+1],0.5*ones(1,2),'k','LineWidth',2 )
    hold on
    title( titles{jjj} );
    
end

