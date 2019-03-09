clear all
close all
clc

%% Define MAS and tuning parameters
global alpha c L omega gamma lambda epsilon_eta N tau_miet ieta il itau delay

s = [1 1 2 2 3 3 4 5 5  7];
t = [2 8 3 7 4 6 5 6 8  8];
G = graph(s,t);

L = laplacian(G);
D = diag(degree(G));
N = length(L);

delta   = 0.05;
mu      = 0.05;
a       = 0.1*ones(N,1);
alpha   = 0.5*ones(N,1);
lambda  = 0.1*ones(N,1);
epsilon_eta = 0.05;

c       = (1-delta)*(1-a.*diag(D));
gamma   = sqrt(diag(D)./a+mu);

%% Compute the (tau_miet,tau_mad)-tradeoff curves
tau_end =-(alpha.*c)./(gamma.*sqrt(alpha.*c)).*atan((lambda.^2-1)./lambda.*(1./(sqrt(alpha.*c)+1./sqrt(alpha.*c))));

for i=1:2
    m=1;
    for lambda_loop = 0.1:0.1:0.4
        phi00   = 1/lambda_loop; 
        k = 1;
        for phi10   = 1.001:0.1:phi00
            phi0    = [phi00;phi10];
            options = odeset('RelTol',1e-4,'AbsTol',1e-6,...
           'Events',@(t,phi) event_phi(t,phi,gamma(i),lambda_loop,phi10));
            dphi   = @(t,phi)   [-gamma(i)*(1/(alpha(i)*c(i))*phi(1)^2+1);
                                 -gamma(i)/lambda_loop*(1/(alpha(i)*c(i))*phi(2)^2+1)];

            [t,phi,tau_mad(k),phie,ie]=   ode45(dphi,[0 tau_end(i)],phi0,options);
            
            dphi   = @(t,phi)   [-gamma(i)*(1/(alpha(i)*c(i))*phi(1)^2+1);
                                 0];

            tau_miet(k) =-(alpha(i).*c(i))./(gamma(i).*sqrt(alpha(i).*c(i))).*atan((lambda_loop*phi10-1./lambda_loop).*(1./(sqrt(alpha(i).*c(i))+phi10*1./sqrt(alpha(i).*c(i)))));
            k = k+1;
        end
        line_type = {'-','--',':','-.'};
        figure(2);subplot(1,2,i);hold on;plot(tau_miet,tau_mad,['k',line_type{m}]);drawnow
        xlabel('$\tau_{miet}^i$ [s]')
        ylabel('$\tau_{mad}^i$ [s]')
        
        grid on;
        clear tau_miet tau_mad tau_miet_2
        m = m+1;
    end
    legend('$\lambda = 0.1$','$\lambda = 0.2$','$\lambda = 0.3$','$\lambda = 0.4$')
end

%% Select a particular (tau_miet,tau_mad)-combination

tau_end =-(alpha.*c)./(gamma.*sqrt(alpha.*c)).*atan((lambda.^2-1)./lambda.*(1./(sqrt(alpha.*c)+1./sqrt(alpha.*c))));
lambda  = 0.2*ones(N,1);
for i=1:N
    phi00   = 1/lambda(i); 
    phi10   = phi00/2.5;
    phi0    = [phi00;phi10];
    options = odeset('RelTol',1e-4,'AbsTol',1e-6,'MaxStep',1e-4,...
   'Events',@(t,phi) event_phi(t,phi,gamma(i),lambda(i),phi10));
    dphi   = @(t,phi)   [-gamma(i)*(1/(alpha(i)*c(i))*phi(1)^2+1);
                         -gamma(i)/lambda(i)*(1/(alpha(i)*c(i))*phi(2)^2+1)];

    [t,phi,tau_mad(i,1),phie,ie]=   ode45(dphi,[0 tau_end(i)],phi0,options);

    dphi   = @(t,phi)   [-gamma(i)*(1/(alpha(i)*c(i))*phi(1)^2+1);
                         0];

    [t,phi,tau_miet(i,1),phie,ie]=   ode45(dphi,[tau_mad(i) tau_end(i)],phie,options);

    line_type = {'-','--',':','-.'};
    if i<3
    figure(2);subplot(1,2,i);hold on;plot(tau_miet(i),tau_mad(i),'o');drawnow
    xlabel('$\tau_{miet}^i$ [s]')
    ylabel('$\tau_{mad}^i$ [s]')

    grid on;

    end
end

% filename = ['Figures\MIET_MAD.tikz'];
% matlab2tikz(filename, 'height', '\figureheight', 'width', '\figurewidth' );

%% Define the flow and jump dynamics

omega    = @(tau) tau<=tau_miet;

Lbar    = [];
for i = 1:N
    Lbar    = blkdiag(Lbar,L(i,:));
end

% xi = (x,e,eta,tau,r,l)
ix          =    1:N;
ie          =    N+1:N^2+N;
iei         =    N+1:N+1:N^2+N;
ieta        =    N^2+N+1:N^2+2*N;
itau        =    N^2+2*N+1:N^2+3*N;
ir          =    N^2+3*N+1:N^2+4*N;
il          =    N^2+4*N+1:2*N^2+4*N;

F           = @(t,xi) [-L*xi(ix)-Lbar*xi(ie);
                        kron(ones(N,1),L*xi(ix)+Lbar*xi(ie));
                        (1-alpha).*c.*(L*xi(ix)+Lbar*xi(ie)).^2-(1-omega(xi(itau))).*gamma.^2.*(1+1./(alpha.*c).*phi10^2.*lambda.^2).*xi(iei).^2-epsilon_eta*xi(ieta);
                        ones(N,1);
                        zeros(N,1);
                        zeros(N^2,1)];
        
Gamma_i     = @(xi) diag(xi(ieta)<=0);
Gamma_bar_i = @(xi) diag(xi(ieta)>=0);
Gamma_ij    = @(xi,delay) diag(xi(il)>=1&kron(ones(N,1),xi(itau))>=delay);

G_0         = @(xi) [       xi(ix);
                            xi(ie);
                            Gamma_bar_i(xi)*xi(ieta)+Gamma_i(xi)*gamma.*phi10.*lambda.*xi(iei).^2;
                            Gamma_bar_i(xi)*xi(itau);
                            Gamma_i(xi)*xi(ix)+Gamma_bar_i(xi)*xi(ir);
                            xi(il) + kron(eye(N),Gamma_i(xi))*ones(N^2,1)];
                
G_1         = @(xi,delay) [ xi(ix);
                            Gamma_ij(xi,delay)*(kron(ones(N,1),xi(ir))-xi(ie)-kron(ones(N,1),xi(ix)))+xi(ie);
                            xi(ieta);
                            xi(itau);
                            xi(ir);
                            xi(il)-Gamma_ij(xi,delay)*ones(N^2,1)];
                        
%% Simulate the system

options = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1e-3,...
   'Events',@event_delays);

xi0_0   = [8;6;4;2;-2;-4;-6;-8;zeros(N^2,1);zeros(N,1);zeros(N,1);zeros(N,1);zeros(N^2,1)];
xi0     = xi0_0;
T_end   = 8;

te      = 0;
T       = [];
X       = [];

for i = 1:N
    TRIGGER{i} = [];
end

while te<T_end
    delay               =   kron(ones(N,1),tau_mad).*rand(N^2,1);
    [t,xi,te,xie,~]     =   ode45(F,[te T_end],xi0,options);

    if ~isempty(te)
        te = te(1);
        xie = xie(1,:);
        T               =   [T;t(t<te)];%
        X               =   [X;xi(t<te,:)];               
    else
        T               =   [T;t(1:end)];
        X               =   [X;xi(1:end,:)];
    end
    
    if ~isempty(te)
        for i = 1:N
            if xie(ieta(i))<0
                TRIGGER{i}  =   [TRIGGER{i} te];
            end
        end
        % Transmission jump
        xi0             =   G_0(xie');
        % Update jump
        xi0             =   G_1(xi0,delay);
    end
end

%% Plot the results

Color  = {'k','r','r','k','r','k','k','r'};
Line    = {'-','--','-.',':','-','--','-.',':'};
Marker = {'kx','ro','r+','k.','rx','ko','k+','r.'};

step = 250;

figure;
for i = 1:N
    figure(5);
    INTER{i} = TRIGGER{i}(2:end) - TRIGGER{i}(1:end-1);
    hold on;plot(TRIGGER{i}(2:end),INTER{i},Marker{i});
    
    
    figure(12);hold on;plot(T(1:step:end),X(1:step:end,ix(i)),[Color{i},Line{i}])

    figure(11);plot(T(1:step:end),X(1:step:end,ieta(i)),[Color{i},Line{i}]);hold on;
    xlabel('Time [s]')
    ylabel('$\eta_i$')
end

figure(5);hold on;plot([0 T_end],tau_miet(1)*ones(2,1),['-',Color{1}]);
xlabel('Time [s]')
ylabel('Inter-event times')
hold on;plot([0 T_end],tau_miet(2)*ones(2,1),['-',Color{2}]);
legend('$\mathcal{A}_1$','$\mathcal{A}_2$','$\mathcal{A}_3$','$\mathcal{A}_4$','$\mathcal{A}_5$','$\mathcal{A}_6$','$\mathcal{A}_7$','$\mathcal{A}_8$','$\tau_{miet}^1$','$\tau_{miet}^2$')
grid on

% filename = ['Figures\IET',num2str(i),'.tikz'];
% matlab2tikz(filename, 'height', '\figureheight', 'width', '\figurewidth' );

figure(12);
xlabel('Time [s]')
ylabel('$x_i$')
grid on
legend('$\mathcal{A}_1$','$\mathcal{A}_2$','$\mathcal{A}_3$','$\mathcal{A}_4$','$\mathcal{A}_5$','$\mathcal{A}_6$','$\mathcal{A}_7$','$\mathcal{A}_8$')

% filename = ['Figures\State_trajectories',num2str(i),'.tikz'];
% matlab2tikz(filename, 'height', '\figureheight', 'width', '\figurewidth' );
                     