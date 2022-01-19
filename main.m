clear all 
close all
clc

%% Reading data
T = readtable('covid_korea_n.csv');
dates_read = datetime(T{:,1},'InputFormat','ddMMMyy');
confirmed_read=T{:,2};
released_read=T{:,3};

%% Daily numbers
confirmed_read=diff(confirmed_read);
released_read=diff(released_read);

%% Proper data retrieval
C=confirmed_read(100:500);
R=released_read(100:500);

%% Visualized of useful data
figure
pause(10)
subplot(2,2,1)
plot(dates_read(100:500),C);
grid on
title('Daily confirmed')
ylabel('Number');
xlabel('Date');

subplot(2,2,2)
plot(dates_read(100:500),R);
grid on
title('Daily released')
ylabel('Number');
xlabel('Date');
pause(0.5)

%% Compartment Model identification
choice_L=50;
sp=subplot(2,2,3);
dx0 = 0;
dy0 = -0.05;
dwithx = 0.4;
dwithy = 0.1;
set(sp,'position',get(sp,'position')+[dx0,dy0,dwithx,dwithy])
for i=1:(500-100)/choice_L
    x1=C((i-1)*choice_L+1:choice_L*i);
    x2=R((i-1)*choice_L+1:choice_L*i);
    x1=x1';
    x2=x2';
    L=choice_L;             
    Fs = 1;            
    f = Fs*(0:(L/2))/L;
    w=2*pi*f;
    w=w';



    base1=abs(fft(ones(length(x1),1)));
    B1=base1(1:L/2+1);
    base2=abs(fft(x1));
    B2=base2(1:L/2+1);
    B2=B2';
    base3=abs(fft(x1.^2));
    B3=base3(1:L/2+1);
    B3=B3';


    case1=abs(fft(ones(length(x2),1)));
    C1=case1(1:L/2+1);
    case2=abs(fft(x2));
    C2=base2(1:L/2+1);
    C2=C2';

    case3=abs(fft(x2.^2));
    C3=case3(1:L/2+1);
    C3=C3';


    X1p=1j*(w).*B1-x1(1)*ones(length(B1),1);
    X2p=1j*(w).*C1-x2(1)*ones(length(C1),1);


    fun = @(p)sum(abs(X1p-p(1)+p(2)*B1-p(3)*C1+p(4)*B3-p(6)*C3)+abs(X2p-p(5)-p(2)*B1+p(3)*C1-p(4)*B3+p(6)*C3));
    p0 = [0,0,0,0,0,0];
    P = fminsearch(fun,p0);
    
    if abs(P(4)/P(2))<0.01

    s = [1 2 3 4];
    t = [2 1 2 1];
    weights = [abs(P(2)) abs(P(3)) abs(P(5)) abs(P(1))];
    G = digraph(s,t,weights);
    LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
    nLabels = {'C','R','',''};
    h=plot(G,'EdgeLabel',G.Edges.Weight,'LineWidth',LWidths,'NodeLabel',nLabels);
    labelText = {join([num2str(P(2)),' C']) join([num2str(P(3)),' R']) num2str(P(5)) num2str(P(1))};
    labeledge(h,s,t,labelText)
    title(join(['Compartment Model: ',datestr(dates_read((i-1)*choice_L+1)),' to ',datestr(dates_read(choice_L*i))]))
    pause(0.5)
    
    else
    s = [1 2 3 4];
    t = [2 1 2 1];
    weights = [abs(P(4)) abs(P(6)) abs(P(5)) abs(P(1))];
    G = digraph(s,t,weights);
    LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
    nLabels = {'C','R','',''};
    h=plot(G,'EdgeLabel',G.Edges.Weight,'LineWidth',LWidths,'NodeLabel',nLabels);
    labelText = {join([num2str(P(4)),' {C^2}']) join([num2str(P(6)),' {R^2}']) num2str(P(5)) num2str(P(1))};
    labeledge(h,s,t,labelText)
    title(join(['Compartment Model: ',datestr(dates_read((i-1)*choice_L+1)),' to ',datestr(dates_read(choice_L*i))]))
    pause(0.5)
    end
    
end


