% AMP function
% Amen Memmi

function xhat=AMP(y,H,sigma2,sigmas2,iterAMP,m,n)
%   AMP detector in Massive MIMO

    r=zeros(m,1);
    xhat=zeros(n,1);
    alpha=sigmas2;%initial estimation variance
        for t=1:iterAMP
        r=y-H*xhat+(n/m)*sigmas2/(sigmas2+alpha)*r;
        alpha=sigma2+(n/m)*sigmas2*alpha/(sigmas2+alpha);
        xhat=(sigmas2/(sigmas2+alpha))*(H'*r+xhat);
        end
    xhat=sign(real(xhat))+sqrt(-1)*sign(imag(xhat));
end