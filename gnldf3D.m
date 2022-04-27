function Imf = gnldf3D(In,N,dt,diffunc)
% Arguments:
%       I           - Input 3D image
%       N           - Number of iterations.
%       dt          - Constant integration, 0<dt<1/4 for a 2D problem with 
%                     4 point neighbors. (Max. value of 0.25 for numerical
%                     stability)
%       diffunc     - 'wregion' = Perona Malik diffusion equation 1
%                     'hgrad' = Perona Malik diffusion equation 2
%                     'sigmoid' = Sigmoid function
%
%
% Returns:
% Imf      - Filtered Image.
%
%
% Usage Example:
% Imf = gnldf(Imf,5,0.25,'wregion');
% 
% Any comments and suggestions to m.arafathussain@gmail.com
% Acknowledgement: This code is patly adopted from a code of Eric Michel Glez 
%************************************************************************
%************************************************************************


if (nargin<1) 
   error('not enough arguments (At least 2 should be given)');
   fprintf(1,'Not enough arguments');
   return;
end

if ~exist('N')
   N=1;
end

if ~exist('dt')
   dt = 0.25;
end

if ~exist('diffunc')
   diffunc = 'wregion';
end

% if ndims(In)==3
%   error('Sorry, DIFFILTERING only operates on 2D grey-scale images');
% end;


%Setting initial conditions

I = In;   
[height,width,elevation] = size(I);


max_i = max(max(I(:)));
I(~isfinite(I)) = max_i; % Eliminate Infinite values

Imf = zeros(height,width,elevation+2);
Imf(:,:,1) = I(:,:,1);
Imf(:,:,elevation+2) = I(:,:,elevation);
Imf(:,:,2:elevation+1) = I;

%------------------------------------------
% Automatic Computation of parameter "d"
%------------------------------------------
    ImAux = padarray(Imf,[1 1],'symmetric','both');    
    ImAux2 = ImAux;%medfilt3(ImAux);

    dxd = (ImAux2(2:height+1,3:width+2,2:elevation+1)-ImAux2(2:height+1,1:width,2:elevation+1))/2;
    dyd = (ImAux2(3:height+2,2:width+1,2:elevation+1)-ImAux2(1:height,2:width+1,2:elevation+1))/2;
    dzd = (ImAux2(2:height+1,2:width+1,3:elevation+2)-ImAux2(2:height+1,2:width+1,1:elevation))/2;

    GI_mag = single((dxd.*dxd+dyd.*dyd+dzd.*dzd).^(1/2));


    % Using Median Absolute Deviation (MAD)
    d = 1.4826*median(abs(GI_mag(:)-median(GI_mag(:))))+eps;

    clear GI_mag;
    clear dxd; clear dyd, clear dzd;
%------------------------------------------



for ni = 1:N
    
    if ni > N
        break;
    end;
    
    % Construct ImAux which is the same as Im but has an extra padding of zeros around it.
    ImAux2 = padarray(Imf,[1 1],'symmetric','both');
    

    Ax = (ImAux2(2:height+1,3:width+2,2:elevation+1) + ImAux2(2:height+1,1:width,2:elevation+1))/2;
    Ay = (ImAux2(3:height+2,2:width+1,2:elevation+1) + ImAux2(1:height,2:width+1,2:elevation+1))/2;
    Az = (ImAux2(2:height+1,2:width+1,3:elevation+2) + ImAux2(2:height+1,2:width+1,1:elevation))/2;

    dE = ImAux2(2:height+1,3:width+2,2:elevation+1)   - Imf(:,:,2:elevation+1);       % difference East
    dW = ImAux2(2:height+1,1:width,2:elevation+1)     - Imf(:,:,2:elevation+1);       % difference West 
    
    dN = ImAux2(1:height,2:width+1,2:elevation+1)     - Imf(:,:,2:elevation+1);       % difference North
    dS = ImAux2(3:height+2,2:width+1,2:elevation+1)   - Imf(:,:,2:elevation+1);       % difference South
    
    dU = ImAux2(2:height+1,2:width+1,1:elevation)     - Imf(:,:,2:elevation+1);       % difference Up
    dD = ImAux2(2:height+1,2:width+1,3:elevation+2)   - Imf(:,:,2:elevation+1);       % difference Down

    
    clear ImAux2;

    Dx = abs(dE-dW) - d;       %In order to smooth small edges, more SNR                        
    Dx(Dx<=0)=0;

    Dy = abs(dN-dS) - d;       %In order to smooth small edges, more SNR                        
    Dy(Dy<=0)=0;
    
    Dz = abs(dU-dD) - d;       %In order to smooth small edges, more SNR                        
    Dz(Dz<=0)=0;

    signo_x = (dE+dW)./(abs(dE+dW)+eps);
    signo_y = (dN+dS)./(abs(dN+dS)+eps);
    signo_z = (dU+dD)./(abs(dU+dD)+eps);

    Gx = Imf(:,:,2:elevation+1) + signo_x.*Dx/2;
    Gy = Imf(:,:,2:elevation+1) + signo_y.*Dy/2;
    Gz = Imf(:,:,2:elevation+1) + signo_z.*Dz/2;

    clear signo_x; clear signo_y; clear signo_z;

    Px  = abs(Gx - Ax);
    Py  = abs(Gy - Ay);
    Pz  = abs(Gz - Az);

    clear Ax; clear Gx;
    clear Ay; clear Gy;
    clear Az; clear Gz;

    %Diffusivity function        
    Cx = single(feval(diffunc,Dx,Px));
    Cy = single(feval(diffunc,Dy,Py));
    Cz = single(feval(diffunc,Dz,Pz));

    clear Dx; clear Py; 
    clear Dy; clear Py; 
    clear Dz; clear Pz;
    
    %Solved by Explicit scheme
    Imf(:,:,2:elevation+1) = Imf(:,:,2:elevation+1) + dt*((dE+dW).*Cx + (dN+dS).*Cy + (dU+dD).*Cz);
       
    
    clear Cx; clear Cy; clear Cz;
    clear dE; clear dW; clear dN; clear dS; clear dU; clear dD;


end



return;




%--------------------------------------------
% SOME DIFFUSIVITY FUNCTIONS                |
%--------------------------------------------

%-----------------------------------------
% Favours wide regions over smaller ones |
%-----------------------------------------
%Perona-Malik 1
function y = wregion(D,P)
    y = 1./(1+ (D./(abs(P) + eps)).^2);
return;

%-----------------------------------------------------
% Favours high contrast edges over low contrast ones |
%-----------------------------------------------------
%Perona-Malik 2
function y = hgrad(D,P)
    y = exp(-(D./(abs(P) + eps)).^4);
return;


%---------------------------------------
% Sigmoid diffusivity function     |
%---------------------------------------
%General Function
function y = sigmoid(D,P)
    y = 1./(1+exp(-(abs(P)-abs(D))*50));
return;


%-----------------------------------------------------
% Favours high contrast edges over low contrast ones |
%-----------------------------------------------------
%Black and Saphiro, no sirve
function y = tukey(D,P)
    y = zeros(size(D));
    if(P > D)
        y = 1/2*(1-(P/D).^2).^2;
    end
        
return;



