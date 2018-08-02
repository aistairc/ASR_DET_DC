% Copyright 2018 National Institute of Advanced Industrial Science and Technology (AIST)
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%    http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

function main()

M = 8*8;

A = double(imread('yaleB01_P00A-005E-10.pgm'));
[H,W]=size(A);

AN = A + 0.5*(128.0*rand(H,W));

Render(1,A,'(a)')
Render(2,AN,'(b)')

Y = Encode(AN,M,H,W);

param.maxiter = 10;
param.alpha = 1e+40;
param.Data = Y;
param.mu = 0;
param.tol = 1e-16;
param.pp = 0;

[X] = AnalysisDenoising_det(Y,param);
B = Decode(X,M,H,W);
Render(3,B,'(c)')

Y2 = Encode(A,M,H,W);
PSNRIn0 = PSNRIn(Y,Y2);
display(PSNRIn0);
PSNRIn1 = PSNRIn(X,Y2);
display(PSNRIn1);

end


function Y = Encode(A,M,H,W)

Ms = sqrt(M);
Hn = H/Ms;
Wn = W/Ms;

Y = zeros(M,H*W/M);
for h = 1:Hn
    for w = 1:Wn
        Y(:,Wn*(h-1)+w) = reshape(A(Ms*(h-1)+1:Ms*(h),Ms*(w-1)+1:Ms*(w)),M,1);
    end    
end

end


function B = Decode(X,M,H,W)

Ms = sqrt(M);
Hn = H/Ms;
Wn = W/Ms;

B = zeros(H,W);
for h = 1:Hn
    for w = 1:Wn
        B(Ms*(h-1)+1:Ms*(h),Ms*(w-1)+1:Ms*(w)) = reshape(X(:,Wn*(h-1)+w),Ms,Ms);
    end    
end

end


function PSNRIn = PSNRIn(X,Y)

PSNRIn = 20*log10(255/sqrt(mean((X(:)-Y(:)).^2)));

end


function Render(i,A,lab)

subplot(1,3,i)
colormap(gray(256))
image(A);
set(gca,'XTick',[],'YTick',[])
title(lab)

end
