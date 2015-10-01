%% plot_spe_data.m

clear all; close all; clc;

% loading data;
disp('Carregando dados...')
matrixdata = load('../dat/reshaped_spe_phi_k_new.dat','-ascii');
ids = [ matrixdata(:,1) matrixdata(:,2) matrixdata(:,3) ]; 
phi = matrixdata(:,4);
kx  = matrixdata(:,5);
ky  = matrixdata(:,6);
kz  = matrixdata(:,7);

splot = 2; 
 
ia = input('Escolha Imin do Volume de Interesse:\n');
ib = input('Escolha Imax do Volume de Interesse:\n');
ja = input('Escolha Jmin do Volume de Interesse:\n');
jb = input('Escolha Jmax do Volume de Interesse:\n');
ka = input('Escolha Zmin do Volume de Interesse:\n');
kb = input('Escolha Zmax do Volume de Interesse:\n');

tic 
vals =  [ia ib ja jb ka kb];
if ~isempty( find(vals<=0, 1) )
    error('Limites devem ser valores nao-negativos!'); 
elseif ( ib > 60 || jb > 220 || kb > 85 )
    error('Um dos limites foi excedido: Imax, ou Jmax, ou Kmax.'); 
end

I = ia:ib;
J = ja:jb;
K = ka:kb;
[KK,JJ,II] = meshgrid(K,J,I);

%{ 
   Alocacao: 
   indices J,K sao trocados pois meshgrid expande K em colunas
   e J em linhas. Entao, para alocar corretamente, 
   devemos chamar 'zeros(J,K,.) para J linhas e K colunas.
   
%}
%PHI = zeros(


lines = [];
for k=ka:kb
    for j=ja:jb 
        for i=ia:ib            
            line = find( ( ids(:,1) == k ) & ( ids(:,2) == j ) & ( ids(:,3) == i )  );
            lines = [lines; line];
            
        end   
    end
end

count = 1; 
for j = 1:length(J)
    for k = 1:length(K);
        for i = 1:length(I);
            PHI(j,k,i) = phi( lines(count) ); 
            KX(j,k,i)  = kx ( lines(count) );
            KY(j,k,i)  = ky ( lines(count) );
            KZ(j,k,i)  = kz ( lines(count) );
            count = count + 1;
        end
    end
end

% todas as fatias do VOI
islice = I(1):I(end);
jslice = J(1):J(end);
kslice = K(1):K(end);

switch splot
    case 1
        subplot(2,2,1)
        title('Porosity')
        slice(KK,JJ,II,PHI,kslice,jslice,islice)
        colorbar
    case 2
        title('Porosity')
        slice(KK,JJ,II,PHI,jslice,kslice,islice)
        colorbar
end
toc