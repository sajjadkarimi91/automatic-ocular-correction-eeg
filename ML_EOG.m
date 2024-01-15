function [Z, andis_EOG, andis_WEOG]=ML_EOG(EEG,fs,andissEOG)
%EEG is matrix of pure data that must be EOG corrected
%EEG is N*T (N is # of channels && T is # of time samples)
% fs is sampling frequency
%EEGpower is average power of EEG signal without EOG
N=size(EEG,1);
L1=round(fs/2);
Overlap=round(fs/20)*2;
WO=ones(N,1)*(linspace(0,1,Overlap));
WWO=ones(N,1)*(linspace(1,0,Overlap));
vv=[zeros(N,L1) double(EEG) zeros(N,L1)];
EEG=vv;
clear vv
T=size(EEG,2);
Z=zeros(N,T);
%landa must determined with fs
landa=0.997;
W=0.1*eye(N);
%Esitmation power of EEG by sliding window that has length 2L1+1
u=[];
andis_EOG=[];
andis_WEOG=[];
counter1=L1;
k_1=L1+1;
%Estimation of EEG average power
V=zeros(1,N);
for i=1:N
    [a, b]=hist(EEG(i,:),500);
    f = fit(b',a','gauss1');
    V(i)=f.c1^2;
end
EEGpower=(max(V)+mean(V(andissEOG)'));

clear a b f V
alfa=ones(N,1);
while counter1<T-L1
    counter1=counter1+1;
    %At first divide time to two part 1- time with EOG 2- time without EOG
    %Power estimation

    Mn=mean(EEG(andissEOG,counter1-L1:counter1+L1).');
    P=EEG(andissEOG,counter1-L1:counter1+L1)-Mn.'*ones(1,2*L1+1);
    sigma=max(diag((1/(2*L1+1))*P*(P.')));

    if sigma<=EEGpower
        Z(:,counter1)=EEG(:,counter1)-(mean(EEG(:,counter1-L1:counter1+L1).').');
    end

    k1=counter1;
    while((counter1<(T-L1))&&(sigma>EEGpower))
        counter1=counter1+1;
        Mn=mean(EEG(andissEOG,counter1-L1:counter1+L1).');
        P=EEG(andissEOG,counter1-L1:counter1+L1)-Mn.'*ones(1,2*L1+1);
        [sigma, maxvar]=max(diag((1/(2*L1+1))*P*(P.')));
    end
    k2=counter1;
    %after detection of segment with EOG and without EOG, correction method
    %begins 1- remove baseline 2- estimate u then estimate alfa and b
    if (k2-k1>fs/10)
        if(size(u,1)~=0)
            Z(:,k_1+1:Overlap+k_1)=(EEG(:,k_1+1:Overlap+k_1)-alfa*(u(end-Overlap+1:end).')).*(WWO)+Z(:,k_1+1:Overlap+k_1).*WO;
        end
        for i=k_1-L1+1:k1-L1-1
            W=UPW(W,Z(:,i),landa);
        end
        maxvar=andissEOG(maxvar);
        andis_EOG=[andis_EOG k1-L1:k2-L1];
        andis_WEOG=[andis_WEOG k_1-L1+1:k1-L1-1];

        gg=EEG(:,k1-Overlap:k2+Overlap);
        % Now some denoising used to better estimation of EOG component

        XX=gg;
        %u is vector k2-k1*1
        u=(XX(maxvar,:).')/norm(XX(maxvar,:));
        alfa_1=0;
        while norm(alfa-alfa_1)>0.0001
            %alfa is N*1
            alfa_1=alfa;
            u=((XX.')*W*alfa)/((alfa.')*W*alfa);
            alfa=(XX*u);
            alfa=normc(alfa);
        end
        Z(:,k1:k2)=(EEG(:,k1:k2)-alfa*(u(Overlap+1:end-Overlap).'));
        Z(:,k1-Overlap:k1-1)=(EEG(:,k1-Overlap:k1-1)-alfa*(u(1:Overlap).')).*(WO)+Z(:,k1-Overlap:k1-1).*WWO;
        k_1=k2;
        plot(alfa)
        hold on
    else

        for kkk=k1:k2
            Z(:,kkk)=EEG(:,kkk)-(mean(EEG(:,kkk-L1:kkk+L1).').');
        end
    end

end
aa=Z(:,L1+1:end-L1);
clear EEG Z
Z=aa;

function W=UPW(w,x,landa)
K=(w*x/(1+(x'*w*x)/landa))/landa;
W=(w-(K*(x')*w))/landa;


