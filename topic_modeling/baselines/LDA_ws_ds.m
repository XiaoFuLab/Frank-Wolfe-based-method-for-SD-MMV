function [WS DS] = LDA_ws_ds(D)

% input 
% D document word-by-doc

[M,N]=size(D);

T = N*M; % total number of samples

mm = 1;
WS=[];
DS=[];
% for m=1:M
%    for n=1:N
%        disp(['running at ',num2str(m*n),' of total ',num2str(T)])
%        zz = D(m,n);
%        WSinc  = m*ones(zz,1);
%        WS = [WS;WSinc];
%        DSinc = n*ones(zz,1);
%        DS = [DS;DSinc];
%    end
% end


   for n=1:N
%        disp(['running at ',num2str(n),' of total ',num2str(N)])
       
       nnz_id = find(D(:,n)>0);
       
       for ii=1:length(nnz_id)
           
           zz = D(nnz_id(ii),n);
           WSinc  = nnz_id(ii)*ones(zz,1);
           WS = [WS;WSinc];
           DSinc = n*ones(zz,1);
           DS = [DS;DSinc];
       end
   end
