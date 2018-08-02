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

function [X] = AnalysisDenoising_det(Y,param)


TrueOmega=Y(1:64,1:64);
TrueDictionary=eye(64,64);
 U0=rand(size(TrueOmega));
 U0=normcols(U0);
 V0=rand(64,504);
 V0=normrows(V0);
Y0=Y;
X=Y;
W0=TrueDictionary;
Omega=TrueOmega+param.pp*U0;
H0=V0;



A_temp=U0*V0;
B_temp=param.Data'*A_temp;
C_temp=trace(B_temp);
D_temp=(norm(A_temp,'fro'))^2;
E_temp=C_temp/D_temp;



Y=Y0;


N=(size(W0,2));
sparsityIter=[];
costtimeIter = [];
iterNums = [];
SNR_1=[];
erro=[];
tStart=tic;

% calcuate the initial H_k
 if det(H0*H0')==0;
      H_k=zeros(size(H0*H0',1),size(H0,2));
      H_k(:)=0;    
    elseif det(H0*H0')>0;
        temp=H0*H0';
        temp_det=det(temp);
        temp_inv=inv(temp);
        temp_inv_H=temp_inv*H0;
        H_k=2*param.alpha* temp_det*temp_inv_H;
 end
 %%
 
 
 H=H0;
 W=W0;
 

%%-----------main loop---------------------------------

for iter = 1:param.maxiter
  
%---------- update W ---------------------------------


 V_1=ADCA_DIC2(Y,H,W,param.mu);

W=W0; 
W=normrows(W);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------- improve H -------using w_(k) to updata H-------

%%=========algorithm by using HALS-DCA=======%%

M=Y';V=W';U=H';
R_1=V*V';
MV=M*V';
VV=R_1;


for jj=1:N
    Eab_1=R_1;
    Eab_1(jj,jj) = 0;
    U(:,jj) = max(  (   (   MV(:,jj)-U*Eab_1(:,jj) + 0.5* (H_k(jj,:))'  )/VV(jj,jj)  ), param.tol  );

end

H=U';
H=normrows(H);


%=========Calcuating H_k=========

if det(H*H')==0;
   H_k=zeros(size(H*H',1),size(H,2));
   H_k(:)=0;    
elseif det(H*H')>0;
   temp=H*H';
   temp_det=det(temp);
   temp_inv=inv(temp);
   temp_inv_H=temp_inv*H;
   H_k=2*param.alpha* temp_det*temp_inv_H;
end

  
%% Evaluate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 Erro = norm(Y-W*H,'fro')^2/(norm(Y,'fro')^2); 
    SNR=10*log10(1/Erro);
    SNR_1=[SNR_1 SNR];
    erro=[erro Erro];
     sparsity = Sparsity_Hoyer(H);


    sparsityIter=[sparsityIter sparsity];

    costtime=toc(tStart);
    costtimeIter = [costtimeIter costtime];
    
    Omega=H*pinv(Y0);
    Omega = ProjUNColBall(Omega')';    
    %Omega=TrueOmega;
    
    %Omega=pinv(W);
    [ratio(iter),ErrorBetweenDictionaries] = I_findDistanseBetweenDictionaries(normcols(TrueOmega'),normcols(Omega'));
%     if ratio==100
%         break;
%     end
%     if rem(iter,100)==0
%         disp(['The FASTNMF-L1 algorithm retrived ',num2str(ratio),' atoms from the original dictionary']);
%     end
% ratioIter = [ratioIter ratio];
   
 X=pinv(TrueDictionary+0.5*(Omega')*Omega)*(Y+0.5*Omega'*H);
 
end

parm.erro=erro;
parm.SNR=SNR_1;
parm.sparsityIter=sparsityIter;
parm.costtimeIter = costtimeIter;
parm.ratioIter = ratio;
parm.iterNums = iterNums;

%% BP method%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A=param.TrueDictionary;
% H_l1=[];
% parm.sparsity_l1=[];
% tic
% 
% for t=1:1000
% H_l1(:,t)= SolveBP(A, param.Data(:,t), 50,100);
% sparsity_l1 = Sparsity_Hoyer(H_l1);
% parm.sparsity_l1=[parm.sparsity_l1 sparsity_l1];
% end
% 
% plot(parm.sparsity_l1)
% toc
end






function V_1=HALS_DIC(Y,H,W,mu) % algorithm from myself by using quadritic function

A=Y';U=H';V_0=W;
 
MM=U'*U; NN=A'*U;

[n,r]=size(V_0);

NablaF=V_0*MM-NN;

V_1=[];

for i=1:n
    for j=1:r
        V_1(:,j)=max(V_0(:,j)-NablaF(:,j)/(max(norm(U(:,j),'fro'),mu)),1e-8);
    end
end

end



function V_1=ADCA_DIC(Y,H,W,mu)
A=Y';
U=H';
V_0=W;
 
MM=U'*U;
NN=A'*U;
[n,r]=size(V_0);

NablaF=V_0*MM-NN;
 
DI=diag(MM);
Delta_V=[];
k=0;
flag=0;
for i=1:n
    for j=1:r
        
       V_1(i,j)=min(V_0(i,j),NablaF(i,j)/max(DI(j,1),mu)); %The updated value.

       Delta_V(i,j)=V_0(i,j)-V_1(i,j); %displacement.
       
       delta_F(i,j)= NablaF(i,j)*Delta_V(i,j)-0.5*MM(j,j)*(Delta_V(i,j))^2;
    end
end
       %F=norm(A-U*V_0','fro')-param.alpha*det(U'*U); %the value of object function. 

for i=1:n
      
      for j=1:r 
       deltaF_max=max(delta_F(i,:));
      
       k=1;
       
       V_1max=max(V_1(i,:));
       
       V_temp=0.5*(V_1max+V_1(i,:));% the updated t(v) and
       
       deltaJ =(norm(A(:,j)-U(:,j)*V_temp(1,j),2))^2-(norm(A(:,j)-U(:,j)*V_0(i,j)',2))^2;
      
       
       if deltaJ>deltaF_max
           
       flag=1;
       
       end
       if flag==1, k=1;       
       
       D(i,j)=delta_F(i,j);
        
        else 
     
       D(i,j)=0;
        end
      end
end

D_1=max(D*MM,param.mu);
    
for i=1:n
    for j=1:r
        if D(i,j)==0
           V_l(i,j)=V_0(i,j);
        end
        if D(i,j)>0
           V_1(i,j)=V_0(i,j)-D(i,j).*NablaF(i,j)/D_1(i,j);
        end
       
    end
end
end



function V_1=SHALS_Dic(Y,H,W_0,tol)% method form Tang
    E = H*H';
    YH = Y*H';
    HH = E;
    for i=1:size(W_0,2)
        Eab = E;
        Eab(i,i) = 0;
%         W(:,i) = max( (VH(:,i)-W*Eab(:,i))/HH(i,i), tol);
        nw = sqrt(sum(((YH(:,i)-W_0*Eab(:,i))/HH(i,i)).^2));
        V_1(:,i) = max( ((YH(:,i)-W_0*Eab(:,i))/HH(i,i))/nw, tol);
    end
end

function V_1=ADCA_DIC2(Y,H,W,mu) % method form DC extension
A=Y';
U=H';
V_0=W;
 
MM=U'*U;
NN=A'*U;
[n,r]=size(V_0);

NablaF=V_0*MM-NN;
 
D=V_0;
D_1=max(D*MM,mu);
    
for i=1:n
    for j=1:r
        if D(i,j)==0
           V_l(i,j)=V_0(i,j);
        end
        if D(i,j)>0
           V_1(i,j)=V_0(i,j)-D(i,j).*NablaF(i,j)/D_1(i,j);
        end
       
    end
end

end
