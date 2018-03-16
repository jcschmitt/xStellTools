function [ output_args ] = plot_fieldlines(data,varargin)
%PLOT_FIELDLINES(data,[plottype])  Plots the data from read_mgrid
%   The PLOT_FIELDLINES routine plots data read by READ_FIELDLINES.  There 
%   are various plotting options.
%   Options:
%       'basic':        Poincare plot on the first cutplane.
%       'color':        Specify color ('color','r')
%       'cutplane':     Poincare plot on specific cutplane ('cutplane',5)    
%       'camera':       Make a camera image by binning poincare points.
%       'camera_AEV30': W7-X AEV30 view (from AEQ21)
%       'strike_2D':    Strike pattern
%       'wall_strike':  Strucutre strike heat map
%
%   Usage:
%       line_data=read_fieldlines('fieldlines_test.h5');
%       plot_fieldlines(line_data);
%
%   See also read_fieldlines.
%
%   Created by: S. Lazerson (lazerson@pppl.gov)
%   Version:    1.22
%   Date:       11/19/14


% Initialize some variables
plottype=0;
nphi=1;
npoinc = data.npoinc;
nsteps = data.nsteps;
nlines = data.nlines;
line_color='k';
% Handle varargin
if nargin > 1
    for i=1:nargin-1
        switch varargin{i}
            case 'basic'
                plottype=0;
            case 'camera'
                plottype=1;
            case 'camera_frame'
                plottype=2;
            case 'camera_AEV30_old'
                plottype=3;
            case 'camera_AEV30'
                plottype=4;
            case 'camera_AEV30_3D'
                plottype=5;
            case 'cutplane'
                i=i+1;
                nphi=varargin{i};
            case 'strike_2D'
                plottype=6;
            case 'wall_strike'
                plottype=7;
            case 'color'
                i=i+1;
                line_color=varargin{i};
        end
    end
end

switch plottype
    case{0}
        line_dex = nphi:npoinc:nsteps;
        for i=1:10:nlines
            hold on;
            plot(data.R_lines(i,line_dex),data.Z_lines(i,line_dex),'.','Color',line_color,'MarkerSize',0.1);
            hold off;
        end
        if isfield(data,'Rhc_lines')
            line_dex = nphi:npoinc:size(data.Rhc_lines,2);
            for i=1:size(data.Rhc_lines,1)
                hold on
                plot(data.Rhc_lines(i,line_dex),data.Zhc_lines(i,line_dex),'.r');
                hold off
            end
            hold on
            plot(data.Rhc_lines(1,nphi),data.Zhc_lines(1,nphi),'+r');
            hold off
        end
        axis equal
    case{1}
        camera=[1024 1024];
        line_dex = nphi:npoinc:nsteps;
        r=data.R_lines(:,line_dex);
        z=data.Z_lines(:,line_dex);
        r=reshape(r,[numel(r) 1]);
        z=reshape(z,[numel(z) 1]);
        % Use vessel as constraint
        r=[4.298; r; 6.341];
        z=[-0.76; z; 0.76];
        syn=hist3([r z],camera);
        xb=linspace(min(r),max(r),size(syn,1));
        yb=linspace(min(z),max(z),size(syn,1));
        pixplot(xb,yb,syn)
        caxis([0 max(mean(syn))]);
        colormap bone;
        axis equal;
        set(gca,'Color','black');
        %axis tight;
    case{2}
        camera=[1300 1030];
        cx=[0.53 0.91];
        cy=[-0.19 0.19];
        line_dex = nphi:npoinc:nsteps;
        r=data.R_lines(:,line_dex);
        z=data.Z_lines(:,line_dex);
        r=reshape(r,[numel(r) 1]);
        z=reshape(z,[numel(z) 1]);
        z = z(r>0.541);
        r = r(r>0.541);
        syn=hist3([r z],camera);
        xb=linspace(min(r),max(r),size(syn,1));
        yb=linspace(min(z),max(z),size(syn,1));
        pixplot(xb,yb,syn)
        caxis([0 max(mean(syn))]);
        colormap bone;
        axis equal;
        xlim(cx);
        ylim(cy);
    case{3} % AEV30 W7-X
        camera=[1392 1024];
        n_cam = [-0.91206066 -0.37672560 -0.16193573];
        x_cam = [1.69477886 6.12453262 0.64880745];
        a_cam = 32.21906432;
        u_cam = [-0.24104390  0.17308427  0.95495532]; %camroll(14.59570350);
        dex = and(data.R_lines>0,data.PHI_lines <= 2*pi);
        dex = and(dex,data.R_lines>data.raxis(1));
        dex = and(dex,data.R_lines<data.raxis(end));
        dex = and(dex,data.Z_lines>data.zaxis(1));
        dex = and(dex,data.Z_lines<data.raxis(end));
        X   = data.X_lines(dex);
        Y   = data.Y_lines(dex);
        Z   = data.Z_lines(dex);
        % A failed attempt to do a weighting
        %phi = data.PHI_lines;
        %n   = data.nsteps-2;
        %for i=1:data.nlines
        %    phi(i,2:data.nsteps)=1:-1/double(n):0;
        %end
        %phi = phi(dex);
        %phi = phi(2:end);
        % New Stuff
        [x_im,  y_im] = points_to_camera(X(2:end),Y(2:end),Z(2:end),...
            'camera',camera,...
            'fov',a_cam,'camera_pos',x_cam,'camera_normal',n_cam,...
            'camera_up',u_cam);
        x_max = max(x_im);
        y_max = max(y_im);
        x_min = min(x_im);
        y_min = min(y_im);
        syn=hist3([x_im y_im],'nbins',[round(x_max-x_min) round(y_max-y_min)]);
        xb=linspace(x_min,x_max,size(syn,1));
        yb=linspace(y_min,y_max,size(syn,2));
        pixplot(xb,yb,syn)
        caxis([0 max(mean(syn))]);
        set(gcf,'Units','pixels','Position',[1 1 camera]);
        set(gca,'Units','pixels','Color','black','Position',[1 1 camera]);
        xlim([1 camera(1)]);
        ylim([1 camera(2)]);
        colormap hot;
        hold on; plot([333 1336],[376 622],'r','LineWidth',4.0); % plot z=0
    case{4} % AEV30 W7-X
        camera=[1392 1024];
        n_cam = [-0.91206066 -0.37672560 -0.16193573];
        x_cam = [1.69477886 6.12453262 0.64880745];
        a_cam = 32.21906432;
        u_cam = [-0.24104390  0.17308427  0.95495532]; %camroll(14.59570350);
        syn = zeros(camera);
        for i = 2:data.nsteps
            dex = and(data.R_lines(:,i)>0,data.PHI_lines(:,i) <= 2*pi);
            dex = and(dex,data.R_lines(:,i)>data.raxis(1));
            dex = and(dex,data.R_lines(:,i)<data.raxis(end));
            dex = and(dex,data.Z_lines(:,i)>data.zaxis(1));
            dex = and(dex,data.Z_lines(:,i)<data.raxis(end));
            X   = data.X_lines(dex,i);
            Y   = data.Y_lines(dex,i);
            Z   = data.Z_lines(dex,i);
            if isempty(X), continue; end;
            [x_im,  y_im] = points_to_camera(X(2:end),Y(2:end),Z(2:end),...
                'camera',camera,...
                'fov',a_cam,'camera_pos',x_cam,'camera_normal',n_cam,...
                'camera_up',u_cam);
            dex   = x_im <= camera(1);
            x_im  = x_im(dex);
            y_im  = y_im(dex);
            dex   = y_im <= camera(2);
            x_im  = x_im(dex);
            y_im  = y_im(dex);
            dex   = x_im >= 1;
            x_im  = x_im(dex);
            y_im  = y_im(dex);
            dex   = y_im >= 1;
            x_im  = x_im(dex);
            y_im  = y_im(dex);
            if isempty(x_im), continue; end;
            x_max = max(x_im);
            y_max = max(y_im);
            x_min = min(x_im);
            y_min = min(y_im);
            syn_temp=hist3([x_im y_im],'nbins',[round(x_max-x_min) round(y_max-y_min)]);
            xb=linspace(x_min,x_max,size(syn_temp,1));
            yb=linspace(y_min,y_max,size(syn_temp,2));
            syn(round(xb),round(yb))=syn_temp./double(i)+syn(round(xb),round(yb));
        end
        % New Stuff
        pixplot(syn)
        caxis([0 max(mean(syn))]);
        set(gcf,'Units','pixels','Position',[1 1 camera]);
        set(gca,'Units','pixels','Color','black','Position',[1 1 camera]);
        xlim([1 camera(1)]);
        ylim([1 camera(2)]);
        colormap hot;
        hold on; plot([333 1336],[376 622],'r','LineWidth',4.0); % plot z=0
    case{5} % AEV30 W7-X (3D plot)
        % OLD STUFF
        n_cam = [-0.91206066 -0.37672560 -0.16193573];
        x_cam = [1.69477886 6.12453262 0.64880745];
        a_cam = 32.21906432;
        u_cam = [-0.24104390  0.17308427  0.95495532]; %camroll(14.59570350);
        dex = and(data.R_lines>0,data.PHI_lines <= 2*pi);
        X   = data.X_lines(dex);
        Y   = data.Y_lines(dex);
        Z   = data.Z_lines(dex);
        campos(x_cam);
        camtarget(x_cam+n_cam);
        camva(a_cam);
        camup(u_cam);
        %axis equal;
        axis off;
        camproj('perspective');
        hold on;
        plot3(X,Y,Z,'.');
    case{6} % Strike_2D
        r=data.R_lines(:,2);
        z=data.Z_lines(:,2);
        r0 = 9.8;
        rr=r-r0;
        u = atan2(z,rr);
        u = u(isfinite(u));
        v = mod(data.PHI_lines(:,2),data.phiaxis(end));
        v = v(isfinite(v));
        u(u<0) = u(u<0) + 2*pi;
        v(v<0) = v(v<0) + data.phiaxis(end);
        %plot(v,u,'.')
        syn=hist3([v u],[100 100]);
        xb=linspace(min(v),max(v),size(syn,1));
        yb=linspace(min(u),max(u),size(syn,1));
        pixplot(xb,yb,syn)
        colormap hot;
        axis square;
        xlim([min(v) max(v)]);
        ylim([min(u) max(u)]);
    case{7} % Wall strike heat map
        dex1 = data.wall_faces(:,1);
        dex2 = data.wall_faces(:,2);
        dex3 = data.wall_faces(:,3);
        V0   = data.wall_vertex(dex3,:)-data.wall_vertex(dex1,:);
        V1   = data.wall_vertex(dex2,:)-data.wall_vertex(dex1,:);
        FNx  = V1(:,2).*V0(:,3)-V1(:,3).*V0(:,2);
        FNy  = V1(:,3).*V0(:,1)-V1(:,1).*V0(:,3);
        FNz  = V1(:,1).*V0(:,2)-V1(:,2).*V0(:,1);
        heat = 2*double(data.wall_strikes)./sqrt(FNx.*FNx+FNy.*FNy+FNz.*FNz);
        output_args{1}=patch('Vertices',data.wall_vertex,'Faces',data.wall_faces,'FaceVertexCData',heat,'LineStyle','none','CDataMapping','scaled','FaceColor','flat');
        
end
end
