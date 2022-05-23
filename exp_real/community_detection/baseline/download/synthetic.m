clear
clc


n = 1000;
% k = 10;

maxitr = 10;

%% vary rank

err1 = zeros(5,10);

for itr = 1:maxitr
for t = 1:10
    
    k = 5*t;

%     M = gamrnd(1/k,1,k,n);
    M = -log(rand(k,n)) .* (rand(k,n)<0.5);
    s = find(sum(M)==0);
    M(:,s) = 1;
    M = bsxfun(@times,M,1./sum(M));


%     B = eye(k) 
    A = M'*M;

    [M1,~] = GeoNMF(A,k);

    [M2,~] = SPOC(A,k);
    
    
    
%     M31 = anima_tensor( A, G, k, 1 );
%     M32 = anima_tensor( A, G, k, k );
    T4 = cpdgen({M(:,1:300)',M(:,301:600)',M(:,601:1000)'});
    Mc = cpd(T4,k);
    d1 = Mc{1}\ones(size(Mc{1},1),1);
    d2 = Mc{2}\ones(size(Mc{2},1),1);
    d3 = Mc{3}\ones(size(Mc{3},1),1);
    M4 = [ diag(d1)*Mc{1}', diag(d2)*Mc{2}', diag(d3)*Mc{3}' ];
    M4 = bsxfun(@times,M4,sign(sum(M4)));
    M4 = max(M4,eps);
    M4 = bsxfun(@times,M4,1./sum(M4));
    
    
    
    G0 = (1/k)^3*ones(k,k,k)/6;
    for l=1:k
        G0(l,:,:) = squeeze(G0(l,:,:)) + (1/k)^2*eye(k)/6;
        G0(:,l,:) = squeeze(G0(:,l,:)) + (1/k)^2*eye(k)/6;
        G0(:,:,l) = squeeze(G0(:,:,l)) + (1/k)^2*eye(k)/6;
        G0(l,l,l) = G0(l,l,l) + 2*(1/k)^2/6 + 2/k/6;
        G0(l,:,:) = squeeze(G0(l,:,:)) - eye(k)/(k+2)/(k+1)/k;
        G0(:,l,:) = squeeze(G0(:,l,:)) - eye(k)/(k+2)/(k+1)/k;
        G0(:,:,l) = squeeze(G0(:,:,l)) - eye(k)/(k+2)/(k+1)/k;
    end
    G0 = G0 - 1/(k+2)/(k+1)/k;
    T5 = tmprod(G0,{M(:,1:300)',M(:,301:600)',M(:,601:1000)'},1:3);
    Mc = cpd(T5,k);
    d1 = Mc{1}\ones(size(Mc{1},1),1);
    d2 = Mc{2}\ones(size(Mc{2},1),1);
    d3 = Mc{3}\ones(size(Mc{3},1),1);
    M5 = [ diag(d1)*Mc{1}', diag(d2)*Mc{2}', diag(d3)*Mc{3}' ];
    M5 = bsxfun(@times,M5,sign(sum(M5)));
    M5 = max(M5,eps);
    M5 = bsxfun(@times,M5,1./sum(M5));
    
    
    Y = randn(k-1,n)*A;
    [M3,~] = CDMVS(Y);
    
    
    M0 = M;
    M1 = M1';
%     M0 = M(:,G);
%     M1 = M1(G,:)';
%     M2 = M2(:,G);
    
    P = Hungarian(-M1*M0'); M1 = P'*M1;
    err1(1,t) = err1(1,t) + norm(M1(:)-M0(:),1)/n;
    P = Hungarian(-M2*M0'); M2 = P'*M2;
    err1(2,t) = err1(2,t) + norm(M2(:)-M0(:),1)/n;
    
%     P = Hungarian(-M31*M0'); M31 = P'*M31;
%     err(3,t) = norm(M31(:)-M0(:),1)/n;
%     P = Hungarian(-M32*M0'); M32 = P'*M32;
%     err(4,t) = norm(M32(:)-M0(:),1)/n;
    
    P = Hungarian(-M3*M0'); M3 = P'*M3;
    err1(3,t) = err1(3,t) + norm(M3(:)-M0(:),1)/n;
    
    P = Hungarian(-M4*M0'); M4 = P'*M4;
    err1(4,t) = err1(4,t) + norm(M4(:)-M0(:),1)/n;
    P = Hungarian(-M5*M0'); M5 = P'*M5;
    err1(5,t) = err1(5,t) + norm(M5(:)-M0(:),1)/n;

end
end

err1 = err1/maxitr;

%% change B

k = 10;

err2 = zeros(5,10);
for itr = 1:maxitr
for t = 1:10

%     M = gamrnd(1/k,1,k,n);
    M = -log(rand(k,n)) .* (rand(k,n)<0.5);
    s = find(sum(M)==0);
    M(:,s) = 1;
    M = bsxfun(@times,M,1./sum(M));

    
    B = tril(rand(k))*t/20;
    B = eye(k) + B + B';
    A = M'*B*M;

    [M1,~] = GeoNMF(A+n*eye(n),k);

    [M2,~] = SPOC(A,k);
    
    
    G = 1:100;
    
%     M31 = anima_tensor( A, G, k, 1 );
%     M32 = anima_tensor( A, G, k, k );
    T4 = cpdgen({M(:,1:300)'*B,M(:,301:600)'*B,M(:,601:1000)'*B});
    Mc = cpd(T4,k);
    Mc{1} = Mc{1}\A(1:300,1:300);
    Mc{2} = Mc{2}\A(301:600,301:600);
    Mc{3} = Mc{3}\A(601:1000,601:1000);
    d1 = ones(1,size(Mc{1},2))/Mc{1};
    d2 = ones(1,size(Mc{2},2))/Mc{2};
    d3 = ones(1,size(Mc{3},2))/Mc{3};
    M4 = [ diag(d1)*Mc{1}, diag(d2)*Mc{2}, diag(d3)*Mc{3} ];
    M4 = bsxfun(@times,M4,sign(sum(M4)));
    M4 = max(M4,eps);
    M4 = bsxfun(@times,M4,1./sum(M4));
    
    
    
    G0 = (1/k)^3*ones(k,k,k)/6;
    for l=1:k
        G0(l,:,:) = squeeze(G0(l,:,:)) + (1/k)^2*eye(k)/6;
        G0(:,l,:) = squeeze(G0(:,l,:)) + (1/k)^2*eye(k)/6;
        G0(:,:,l) = squeeze(G0(:,:,l)) + (1/k)^2*eye(k)/6;
        G0(l,l,l) = G0(l,l,l) + 2*(1/k)^2/6 + 2/k/6;
        G0(l,:,:) = squeeze(G0(l,:,:)) - eye(k)/(k+2)/(k+1)/k;
        G0(:,l,:) = squeeze(G0(:,l,:)) - eye(k)/(k+2)/(k+1)/k;
        G0(:,:,l) = squeeze(G0(:,:,l)) - eye(k)/(k+2)/(k+1)/k;
    end
    G0 = G0 - 1/(k+2)/(k+1)/k;
    T5 = tmprod(G0,{M(:,1:300)'*B,M(:,301:600)'*B,M(:,601:1000)'*B},1:3);
    Mc = cpd(T5,k);
    Mc{1} = Mc{1}\A(1:300,1:300);
    Mc{2} = Mc{2}\A(301:600,301:600);
    Mc{3} = Mc{3}\A(601:1000,601:1000);
    d1 = ones(1,size(Mc{1},2))/Mc{1};
    d2 = ones(1,size(Mc{2},2))/Mc{2};
    d3 = ones(1,size(Mc{3},2))/Mc{3};
    M5 = [ diag(d1)*Mc{1}, diag(d2)*Mc{2}, diag(d3)*Mc{3} ];
    M5 = bsxfun(@times,M5,sign(sum(M5)));
    M5 = max(M5,eps);
    M5 = bsxfun(@times,M5,1./sum(M5));
    
    
    Y = rand(k-1,n)*A;
    [M3,~] = CDMVS(Y);
    
    
    M0 = M;
    M1 = M1';
%     M0 = M(:,G);
%     M1 = M1(G,:)';
%     M2 = M2(:,G);
    
    P = Hungarian(-M1*M0'); M1 = P'*M1;
    err2(1,t) = err2(1,t) + norm(M1(:)-M0(:),1)/n;
    P = Hungarian(-M2*M0'); M2 = P'*M2;
    err2(2,t) = err2(2,t) + norm(M2(:)-M0(:),1)/n;
    
%     P = Hungarian(-M31*M0'); M31 = P'*M31;
%     err(3,t) = norm(M31(:)-M0(:),1);
%     P = Hungarian(-M32*M0'); M32 = P'*M32;
%     err(4,t) = norm(M32(:)-M0(:),1);
    
    P = Hungarian(-M3*M0'); M3 = P'*M3;
    err2(3,t) = err2(3,t) + norm(M3(:)-M0(:),1)/n;
    
    P = Hungarian(-M4*M0'); M4 = P'*M4;
    err2(4,t) = err2(4,t) + norm(M4(:)-M0(:),1)/n;
    P = Hungarian(-M5*M0'); M5 = P'*M5;
    err2(5,t) = err2(5,t) + norm(M5(:)-M0(:),1)/n;

end
end

err2 = err2/maxitr;


% save('results_Jan22.mat','err1','err2')





subplot(1,2,1)
semilogy(5:5:50,err1(3,:),'-o','LineWidth',2)
hold on
semilogy(5:5:50,err1(1,:),'-x','LineWidth',2)
semilogy(5:5:50,err1(2,:),'-*','LineWidth',2)
semilogy(5:5:50,err1(4,:),'-s','LineWidth',2)
semilogy(5:5:50,err1(5,:),'-d','LineWidth',2)
hold off
xlabel('k: number of community')
ylabel('$\frac{\|M-M^\natural\|_1}{n}$','Interpreter','latex')
axis([5,50,1e-16,1e5])

subplot(1,2,2)
semilogy(0.05:0.05:0.5,err2(3,:),'-o','LineWidth',2)
hold on
semilogy(0.05:0.05:0.5,err2(1,:),'-x','LineWidth',2)
semilogy(0.05:0.05:0.5,err2(2,:),'-*','LineWidth',2)
semilogy(0.05:0.05:0.5,err2(4,:),'-s','LineWidth',2)
semilogy(0.05:0.05:0.5,err2(5,:),'-d','LineWidth',2)
hold off
xlabel('$\max(B_{ij}),i\neq j$','Interpreter','latex')
axis([0.05,0.5,1e-16,1e5])

legend('CD-MVSI','GeoNMF','SPOC','CPD (correct \alpha)','CPD (incorrect \alpha)')