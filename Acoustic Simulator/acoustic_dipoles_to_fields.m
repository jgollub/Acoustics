% Copyright 2013 Evolv Technologies, Inc.  Development of this software was
% supported in part by the U.S. Government under contract number HSHQDC-12-C-00049
function [scene] = acoustic_dipoles_to_fields(panel, locs, speed, verbose)
if nargin < 4
    verbose = 1; % turn on/off print statements
    if nargin < 3
        speed = 'slow'; % use the less memory intensive solver
    end
end
% if a cell array of panels is handed in, loop through each
if iscell(panel)
    scene = cell(size(panel));
    for in=1:length(panel)
        scene{in} = dipoles_to_fields(panel{in}, locs, speed, verbose);
    end
    return;
end

tic;

tol = 1/1000; % tolerance in neglecting weak dipoles
c0 = 343.2; % speed of light
% Z0 = 376.7303; % impedance of free space
% u0 = 4*pi*1e-7; % permeability of free-space

if ~isfield(panel, 'feedLocs')
    error('Panel has not been fed');
end

% extract locations of voxels in the scene
x=locs(:,1);
y=locs(:,2);
z=locs(:,3);
clear locs;

% get important panel params and remove the rest from memory
TotalNumMeasurements = length(panel.f)*panel.numfeeds;
scene.f              = repmat(panel.f, [1, panel.numfeeds]);
scene.numfeeds       = panel.numfeeds;
scene.origin         = locate(panel);
w   = 2*pi*scene.f;
k   = w/c0;
m_x = panel.dipoles.x(:,:);
m_y = panel.dipoles.y(:,:);
m_z = panel.dipoles.z(:,:);
px  = panel.x;
py  = panel.y;
pz  = panel.z;
[Ni, Nj] = size(panel.y);
clear panel;

% fields initialization
P = zeros(TotalNumMeasurements, length(x));
% E_x = zeros(size(E_y));
% E_z = zeros(size(E_y));

% compute the fields at the scene
if strcmp(speed, 'slow')
    % loop through frequencies/switches
    parfor fi = 1:TotalNumMeasurements
        % select only those dipoles with significant energy compared to the
        % strongest dipole at this frequency
        max_m   = max(abs(m_x(fi,:)).^2 + abs(m_y(fi,:)).^2 + abs(m_z(fi,:)).^2, [], 2);
        indices =  (abs(m_x(fi,:)).^2 + abs(m_y(fi,:)).^2 + abs(m_z(fi,:)).^2) > tol*max_m;
        % loop through dipoles in the panel
        for i1 = find(indices)
            % calculate the distance vector between the dipole `i1' and every
            % voxel in the scene
            cx   = x - px(i1);
            cy   = y - py(i1);
            cz   = z - pz(i1);
            r    = sqrt(cx.^2+cy.^2+cz.^2);
            % calculate part of the magnetic dipole Green's function
%             G_md = (-1j*w(fi)*u0)./(4*pi).*(-1j*k(fi)./r-1./(r.^2)).*exp(-1j*k(fi).*r);
              G_md = 1/(4*pi*r) .*exp(-1j*k(fi).*r)
            % add the contribution from dipole `i1' to every voxel in the
            % scene: `m x r' written out explicitly
%             E_x(fi,:) = E_x(fi,:) + (G_md.*(m_y(fi,i1).*cz - m_z(fi,i1).*cy)./r).';
%             E_y(fi,:) = E_y(fi,:) + (G_md.*(m_z(fi,i1).*cx - m_x(fi,i1).*cz)./r).';
%             E_z(fi,:) = E_z(fi,:) + (G_md.*(m_x(fi,i1).*cy - m_y(fi,i1).*cx)./r).';
           
            P(fi) = P(fi) + G_md;
        end
    end
% elseif strcmp(speed, 'fast') % uses too much memory at the moment, do not use
%     cx = (big_sum(x, -px(:)'));
%     cy = (big_sum(y, -py(:)'));
%     cz = (big_sum(z, -pz(:)'));
%     r=sqrt(cx.^2+cy.^2+cz.^2);
%     parfor fi = 1:TotalNumMeasurements
%         % select only those dipoles with significant energy
%         max_m   = max(abs(m_x(fi,:)).^2 + abs(m_y(fi,:)).^2 + abs(m_z(fi,:)).^2, [], 2);
%         indices =  (abs(m_x(fi,:)).^2 + abs(m_y(fi,:)).^2 + abs(m_z(fi,:)).^2) > tol*max_m;
%         for i1 = find(indices)
%             G_md = (-1j*w(fi)*u0)./(4*pi).*(-1j*k(fi)./r(:,i1)-1./r(:,i1).^2).*exp(-1j*k(fi).*r(:,i1));
%             
%             E_x(fi,:) = E_x(fi,:) + (G_md.*(m_y(fi,i1).*cz(:,i1) - m_z(fi,i1).*cy(:,i1))./r(:,i1)).';
%             E_y(fi,:) = E_y(fi,:) + (G_md.*(m_z(fi,i1).*cx(:,i1) - m_x(fi,i1).*cz(:,i1))./r(:,i1)).';
%             E_z(fi,:) = E_z(fi,:) + (G_md.*(m_x(fi,i1).*cy(:,i1) - m_y(fi,i1).*cx(:,i1))./r(:,i1)).';
%         end
%     end
% end
tstop = toc;
if verbose
    fprintf('scene complete in %g minutes.\n', tstop/60);
end
clear m_x m_y m_z;

scene.build_time = tstop;

% store the fields and voxel locations
% scene.E_x=E_x;
% scene.E_y=E_y;
% scene.E_z=E_z;

scene.P=P;

scene.Xmat=x;
scene.Ymat=y;
scene.Zmat=z;