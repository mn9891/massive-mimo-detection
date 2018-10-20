% MMSE - AMP 
% Amen Memmi

%% Code
tic
clc;clear; close all;
n=16;% # of transmitters
m=16;% # of receivers

SNRrange=1:11;
% snrLinear = 10^(0.1*SNRrange);
count=0;
for s=SNRrange
SNRdb=s;
    for monte=1:10000
    x=(2*randi([0,1],n,1)-ones(n,1))+sqrt(-1)*(2*randi([0,1],n,1)-ones(n,1));
    sigmas2=2;%signal variance in QPSK
    H=1/sqrt(2*m)*randn(m,n)+sqrt(-1)/sqrt(2*m)*randn(m,n);
    sigma2=2*n/m*10^(-SNRdb/10); %noise variance in control by SNR in DB
    w=sqrt(2*sigma2)*randn(m,1)+sqrt(-1)*sqrt(2*sigma2)*randn(m,1);
    y=H*x+w; %channel model
 
    %iterAMP is # of iterations in AMP
    iterAMP1=2;
    xhat1=AMP(y,H,sigma2,sigmas2,iterAMP1,m,n);
    iterAMP2=4;
    xhat2=AMP(y,H,sigma2,sigmas2,iterAMP2,m,n);
    iterAMP3=6;
    xhat3=AMP(y,H,sigma2,sigmas2,iterAMP3,m,n);
    
 
     x_mmse=(sigma2/sigmas2*eye(n)+H'*H)^(-1)*H'*y;
     x_mmse=sign(real(x_mmse))+sqrt(-1)*sign(imag(x_mmse));
    errorAMP1(monte)=sum(x~=xhat1);
    errorAMP2(monte)=sum(x~=xhat2);
    errorAMP3(monte)=sum(x~=xhat3);
    errorMMSE(monte)=sum(x~=x_mmse);
    
    end
    count=count+1;
serAMP1(count)=0.1*mean(errorAMP1);
serAMP2(count)=0.1*mean(errorAMP2);
serAMP3(count)=0.1*mean(errorAMP3);
serMMSE(count)=0.1*mean(errorMMSE);
end
figure(2)% plot the SER
semilogy(SNRrange,serAMP1,'-+', SNRrange,serAMP2,'-p',SNRrange,serAMP3,'-o',SNRrange, serMMSE,'-v');
grid on;
legend(['AMP iteration=' int2str(iterAMP1)], ['AMP iteration=' int2str(iterAMP2)], ['AMP iteration=' int2str(iterAMP3)], 'MMSE');
xlabel('SNR (dB)'); ylabel('SER');
title(['BER performance comparison in system m= ' int2str(m)  '  n=' int2str(n)]);
toc



