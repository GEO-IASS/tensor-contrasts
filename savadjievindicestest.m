% testing savadjiev's dispersion and curving indices
%% STEP 0, generate test data


%size of volume
k = 20;
l = 20;
m = 20;

%size of tensor
a = 3;
b = 3;

% test tensor eigenvalues
l1 = 2;
l2 = 1;
l3 = 1;

b1 = [1 0 0];
b2 = [0 1 0];
b3 = [0 0 1];

b1 = b1./norm(b1,2);
b2 = b2./norm(b2,2);
b3 = b3./norm(b3,2);
E = [b1' b2' b3'];

D0 = E*diag([l1 l2 l3])*E';

% rotations
Mx = @(t) [1 0 0; 0 cos(t) -sin(t); 0 sin(t) cos(t)];
My = @(t) [cos(t) 0 -sin(t); 0 1 0; sin(t) 0 cos(t)];
Mz = @(t) [cos(t) sin(t) 0; -sin(t) cos(t) 0; 0 0 1];
rot = @(d,tx,ty,tz) Mx(tx)*My(ty)*Mz(tz)*D0*(Mx(tx)*My(ty)*Mz(tz))';

% generate test function
%test = (1:k)/k*pi; % uniform inc x

%test = (1+cos((1:k)/k*pi))*pi/2; %ridge
%test = repmat(test', [1, l, m]);

%test = atan2(repmat((1:k)'-k/2, 1, l),repmat((1:l)-l/2, k, 1)); %curving
test = pi/2+atan2(repmat((1:k)'-k/2, 1, l),repmat((1:l)-l/2, k, 1)); %dispersion
%test = repmat(test, [1 1 m]);

test = shiftdim(test',-1); %shift XY-test to YZ-test
test = repmat(test, [m 1 1]);

figure(1); clf;
imagesc(test(:,:,10))
colorbar
colormap('gray')

%tensor field;
F = zeros(k,l,m,a,b);

%randomly generate
for x = 1:k
    for y = 1:l
        for z = 1:m
            %F(x,y,z,:,:) = rot(D0, 0, 0, test(x,y,z)); % xy plane rotation
            F(x,y,z,:,:) = rot(D0, test(x,y,z), pi/2, 0); % yz plane rotation
            %F(x,y,z,:,:) = squeeze(F(x,y,z,:,:)) + rand(3,3)/10;
        end
    end
end


%% Optional, read input from xml file


[outVol param]=ReadXml('DispersionXZ.xml');

k = size(outVol,1);
l = size(outVol,2);
m = size(outVol,3);
F = zeros(k,l,m,3,3);

F(:,:,:,1,1) = outVol(:,:,:,1);
F(:,:,:,2,1) = outVol(:,:,:,2);
F(:,:,:,1,2) = outVol(:,:,:,2);
F(:,:,:,3,1) = outVol(:,:,:,3);
F(:,:,:,1,3) = outVol(:,:,:,3);
F(:,:,:,2,2) = outVol(:,:,:,4);
F(:,:,:,2,3) = outVol(:,:,:,5);
F(:,:,:,3,2) = outVol(:,:,:,5);
F(:,:,:,3,3) = outVol(:,:,:,6);


%% STEP 1, find tensor gradient field
% find tensor gradient field

DF = zeros(k,l,m,a,b,3);
lambda = 0; %was 2
for x1 = 1:a
    for x2 = 1:b %can be reduced
        bs = bsarray(F(:,:,:,x1,x2), 'degree', 3, 'lambda', lambda);
        dxbs = partial(bs,1);
        dybs = partial(bs,2);
        dzbs = partial(bs,3);
        dx = indirectFilter(dxbs);
        dy = indirectFilter(dybs);
        dz = indirectFilter(dzbs);
        DF(:,:,:,x1,x2,1) = dx;
        DF(:,:,:,x1,x2,2) = dy;
        DF(:,:,:,x1,x2,3) = dz;
        %DF(1:end-1,:,:,x1,x2,1) = F(2:end,:,:,x1,x2) - F(1:end-1,:,:,x1,x2);
        %DF(:,1:end-1,:,x1,x2,2) = F(:,2:end,:,x1,x2) - F(:,1:end-1,:,x1,x2);
        %DF(:,:,1:end-1,x1,x2,3) = F(:,:,2:end,x1,x2) - F(:,:,1:end-1,x1,x2);
    end
end

%symmetry
%DF(:,:,:,1,2) = DF(:,:,:,2,1);
%DF(:,:,:,1,3) = DF(:,:,:,3,1);
%DF(:,:,:,2,3) = DF(:,:,:,3,2);


%% test contrasts
figure(1)
clf;
hold all;
tr = zeros(k,l,m);
th = zeros(k,l,m);
ph = zeros(k,l,m);

for x = 1:k
    for y = 1:l
        for z = 1:m
            tr(x,y,z) = trace(squeeze(F(x,y,z,:,:)));
            
            [E,L] = eig(squeeze(F(x,y,z,:,:)));
            L = diag(L);
            [~,I] = sort(L,1,'descend');
            e1 = E(:,I(1)); %pick max eig vec
            
            [TH,PHI,~] = cart2sph(e1(1),e1(2),e1(3));
            th(x,y,z) = TH;
            ph(x,y,z) = PHI;
            
            e1 = e1./norm(e1,2);
            plot3([x-e1(1)/2, x+e1(1)/2],[y-e1(2)/2, y+e1(2)/2],[z-e1(3)/2, z+e1(3)/2], 'k');
            
        end
        
        
    end
end

xlabel('x'); ylabel('y'); zlabel('z')
%view([1 0 0]) %view yz plane
view([0 1 0]) %view xz plane
%view([0 0 1]) %view zy plane

axis tight

%%

%showDiffusionEllipsoidxyz(0,0,0,squeeze(F(10,10,20,:,:)))

%clf;
%showDiffusionEllipsoidxyz(0,0,0,rot(D0, 0, 0, test(20,1,1)))


%% STEP 2, calculate tensor parameters

phi1 = zeros(k,l,m,a,b);
phi2 = zeros(k,l,m,a,b);
phi3 = zeros(k,l,m,a,b);

dphi1 = zeros(k,l,m,3);
dphi2 = zeros(k,l,m,3);
dphi3 = zeros(k,l,m,3);

for x = 1:k
    for y = 1:l
        for z = 1:m
            phi1(x,y,z,:,:) = phinorm(1, squeeze(F(x,y,z,:,:)));
            phi2(x,y,z,:,:) = phinorm(2, squeeze(F(x,y,z,:,:)));
            phi3(x,y,z,:,:) = phinorm(3, squeeze(F(x,y,z,:,:)));
            
            dphi1(x,y,z,:) = contr(squeeze(DF(x,y,z,:,:,:)), squeeze(phi1(x,y,z,:,:)));
            dphi2(x,y,z,:) = contr(squeeze(DF(x,y,z,:,:,:)), squeeze(phi2(x,y,z,:,:)));
            dphi3(x,y,z,:) = contr(squeeze(DF(x,y,z,:,:,:)), squeeze(phi3(x,y,z,:,:)));
            
            
            %if x == 10 
            %    return;
            %end
        end
    end
end




%% STEP 3, calculate contrasts

C = zeros(k,l,m);
D = zeros(k,l,m);

for x = 1:k
    for y = 1:l
        for z = 1:m
            C(x,y,z) = curving( squeeze(dphi2(x,y,z,:,:)), squeeze(dphi3(x,y,z,:,:)), squeeze(F(x,y,z,:,:)) );
            D(x,y,z) = dispersion( squeeze(dphi2(x,y,z,:,:)), squeeze(dphi3(x,y,z,:,:)), squeeze(F(x,y,z,:,:)) );
        end
    end
end

%% visualization

figure(2)
clf;
imagesc(C(:,:,5)')
colormap('gray')
%axis ij
%caxis([min(min(C(:)),min(D(:))) max(max(C(:)),max(D(:)))])
colorbar

figure(3)
clf;
imagesc(D(:,:,5)')
colormap('gray')
%axis ij
%caxis([min(min(C(:)),min(D(:))) max(max(C(:)),max(D(:)))])
colorbar

%% print

figure(1)
%print -dpng pevsRidge2


figure(2)
%print -dpng curvingRidge2


figure(3)
%print -dpng dispersionRidge2

%%

figure(2)
%print -dpng curvingCurving3FD

figure(3)
%print -dpng dispersionCurving3FD


