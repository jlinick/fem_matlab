clear

% params
lx = 50e3 ; % length x
ly = 50e3 ; % length y
Emod = 9e9 ; % Young's modulus (Pa)
ht = 200 ; % height (m)
pois = 0.31 ; % Poisson's ratio
%h0 = depth; %5 ; % height of lake (m)
rho_load = 1023.6 ;% density of lake (kg/m3)
rho_top = 919 ; % density of ice (kg/m3)
rho_bot = 1027 ; % density of ocean (kg/m3)
gr = 9.8 ; % gravity b/c I used g and it screwed me up

start_time = 0.;
end_time = 32.;
delta_t = 1.;

t1 = 0.:0.5:30;
t2 = 30.:0.02:31;
times = [t1 t2];%cat(2, 0.:0.25:30, 30.:0.05:31);%[];
depths = [];
times_log = [];

%iterator = 0;
%for time = start_time:delta_t:end_time
for iterator = 1:length(times)
    time = times(iterator);
    %times(end+1) = time;
    h0 = get_depth(time); % height of lake (m)
    times_log(end+1) = time;
    depths(end+1) = h0;
    
    % params (numerical)
    nip = 9 ;   % # of integration points
    ndof = 3 ;  % # of degrees of freedom per node
    nod = 4 ;   % # of nodes per element
    nelx = 80 ; % # of nodes in x-direction
    nely = 80 ; % # of nodes in y-direction

    % setting zero to horizontal forces
    ks = (rho_bot-rho_top)*gr ; % drho g (density diff)
    Px = 0;% in-plane force in x direction (Pa m)
    Pxy = 0 ; % in-plane shear force (Pa m)
    Py = 0 ; % in-plane force in y direction (Pa m)

    length_scale = ht ; % length
    stress_scale = Emod ; % stress
    lx = lx/length_scale ;
    ly = ly/length_scale ;
    Emod = Emod/stress_scale ;
    ht = ht/length_scale ;

    ks = ks*length_scale/stress_scale ;
    Px = Px/(stress_scale*length_scale) ;
    Py = Py/(stress_scale*length_scale) ;
    Pxy = Pxy/(stress_scale*length_scale) ;

    Pm = [ Px Pxy ; Pxy Py ] ;
    x0 = lx/2;%0.1*lx ;
    y0 = ly/2 ;
    sigmax = 0.05*lx ;
    sigmay = 0.075*ly ;

    % calc params
    nx = nelx+1 ; % nodes in x direction
    ny = nely+1 ; % nodes in y direction
    nn = nx*ny ; % total number of nodes
    nels = nelx*nely; % total number of elements
    ntot = ndof*nod ; % element degrees of freedom
    dx = lx/nelx ; % element size in x direction
    dy = ly/nely ; % element size in y direction

    % mesh
    g_coord = zeros(2,nn) ;
    n = 1 ;
    for i = 1:nx
        for j=1:ny
            g_coord(1,n) = (i-1)*dx ;
            g_coord(2,n) = (j-1)*dy ;
            n = n + 1;
        end
    end

    % nodal connections
    gnumbers = reshape(1:nn,[nx ny]) ;
    iel = 1 ;

    for j=1:nelx
        for i=1:nely
            g_num(1,iel) = gnumbers(j,i) ; % node 1
            g_num(2,iel) = gnumbers(j,i+1) ; % node 2
            g_num(3,iel) = gnumbers(j+1,i+1); % node 3
            g_num(4,iel) = gnumbers(j+1,i) ; % node 4
            iel = iel + 1 ;
        end
    end

    % create connection arrays
    sdof = 0 ;
    nf = zeros(ndof,nn) ;
    for n = 1:nn
        for i=1:ndof
            sdof = sdof + 1 ;
            nf(i,n) = sdof ;
        end
    end

    % equation number for elements
    g = zeros(ntot,1) ;
    g_g = zeros(ntot,nels);
    for iel=1:nels ;
        num = g_num(:,iel) ;
        inc=0 ;
        for i=1:nod; for k=1:ndof ; inc=inc+1 ; g(inc)=nf(k,num(i)) ; end ;end
        g_g(:,iel) = g ;
    end

    % B.C.s
    bx0 = find(g_coord(1,:)==0);
    bxn = find(g_coord(1,:)==lx);
    by0 = find(g_coord(2,:)==0);
    byn = find(g_coord(2,:)==ly);
    bcdof = [ nf(1,bxn) nf(3,bxn) nf(2,byn) nf(2,byn) ] ;
    bcval = zeros(1,length(bcdof)) ;

    % integrate

    % local coordinates of integration points
    points(1:3:7,1) = -sqrt(0.6);
    points(2:3:8,1) = 0;
    points(3:3:9,1) = sqrt(0.6);
    points(1:3,2) = sqrt(0.6);
    points(4:6,2) = 0 ;
    points(7:9,2) = -sqrt(0.6);
    w = [ 5./9. 8./9. 5./9.] ;
    v = [ 5./9.*w ; 8./9.*w ; 5./9.*w ] ;
    wts = v(:) ;

    % initialize matrices
    bv = zeros(sdof,1); % rhs load vector
    displ = zeros(sdof,1); % displacement solution
    lhs = sparse(sdof,sdof) ; % stiffness matrix
    xg = reshape(g_coord(1,:),ny,nx); % x-coorindates of nodes
    yg = reshape(g_coord(2,:),ny,nx); % y-coorindates of nodes
    % surface load function defined on nodes
    q0 = h0*rho_load*gr/stress_scale ;
    qw = q0*exp(-((xg-x0).^2/(2*sigmax^2)+(yg-y0).^2/(2*sigmay^2))) ;
    % compute gradients of load
    [FX,FY] = gradient(qw,dx,dy) ;
    % storage of the three load components in vector on nodes
    Qg(nf(1,:)) = qw(:) ;
    Qg(nf(2,:)) = FY(:) ;
    Qg(nf(3,:)) = -FX(:) ;

    % integrate by element
    for iel=1:nels % sum over elements
        num = g_num(:,iel) ;
        coord = g_coord(:,num)' ;
        g = g_g(:,iel) ;
        KM = zeros(ntot,ntot) ;
        KP = zeros(ntot,ntot) ;
        MM = zeros(ntot,ntot) ;
        f = zeros(ntot,1) ;
        % form C matrix for the current element
        is = 1 ; ie = 3 ; % start and end storage indices
        for i=1:nod % loop over nodes
            x = coord(i,1) ; % x coordinates of node
            y = coord(i,2) ; % y coordinates of node
            C(is:ie,:) = [ 1 x y x^2 x*y y^2 x^3 x^2*y ...
            x*y^2 y^3 x^3*y x*y^3 ; ... % w
            0 0 1 0 x 2*y 0 x^2 ...
            2*x*y 3*y^2 x^3 3*x*y^2 ; ...% dwdy
            0 -1 0 -2*x -y 0 -3*x^2 -2*x*y ...
            -y^2 0 -3*x^2*y -y^3 ]; % -dwdx
            is = ie+1 ; % start storage index
            ie = is+2 ; % end storage index
        end
        iC = inv(C); % invert C matrix
        for k=1:nip % integration loop
            xi = points(k,1); % local x-coordinate of integration point
            b = max(coord(:,1)); 
            a = min(coord(:,1)); 
            m = (b-a)/2 ; 
            n = (b+a)/2 ; 
            x = xi*m+n ; % integration point in real x coordinates
            eta = points(k,2); 
            d = max(coord(:,2)); % max y
            c = min(coord(:,2)); % min y
            o = (d-c)/2 ; % half element width in y
            p = (d+c)/2 ; % element center y
            y = eta*o+p ; % integration point
            h = ht + x/lx*ht ; %
            D = Emod*h^3/12/(1-pois^2) ; % bending stiffness
            Dm = D*[1 pois 0 ; pois 1 0 ; 0 0 (1-pois)/2] ; % bending matrix
            Q = [0 0 0 2 0 0 6*x 2*y 0 0 6*x*y 0 ; 0 0 0 0 0 2 0 0 2*x 6*y 0 6*x*y ; 0 0 0 0 2 0 0 4*x 4*y 0 6*x^2 6*y^2 ];
            H = [1 x y x^2 x*y y^2 x^3 x^2*y x*y^2 y^3 x^3*y x*y^3]; % wtf
            R = [0 1 0 2*x y 0 3*x^2 2*x*y y^2 0 3*x^2*y y^3 ; 0 0 1 0 x 2*y 0 x^2 2*x*y 3*y^2 x^3 3*x*y^2] ;
            fun = H * iC ; % shape functions
            B = Q * iC ; % kinematic matrix
            A = R * iC ; % A matrix
            q = fun*Qg(g)' ; % interpolate loads to integration point
            KM = KM + B' * Dm * B * wts(k) ; % elastic stiffness matrix
            KP = KP + A' * Pm * A * wts(k) ; % inplane force matrix
            MM = MM + fun' * fun * ks * wts(k) ; % mass matrix (restoring force)
            f = f + fun'*q*wts(k) ; % load vector
        end
        % element storage of matrix coefficients
        lhs(g,g) = lhs(g,g) + KM + KP + MM ;
        bv(g) = bv(g) - f ;
    end


    % boundary conditions
    lhs(bcdof,:) = 0 ;
    tmp = spdiags(lhs,0) ; tmp(bcdof)=1 ; lhs=spdiags(tmp,0,lhs);
    bv(bcdof) = bcval ;
    % inverse solution
    displ = lhs \ bv ;

    % plots
    f = figure('visible','off');
    subplot(1,2,1);
    %hold on;
    sg = reshape(displ(nf(1,:)),ny,nx);
    colormap parula;
    shading interp;
    surf(xg*length_scale - x0*length_scale,yg*length_scale - y0*length_scale,sg*length_scale)
    zlabel('deflection (m)')
    xlabel('y (m)')
    ylabel('x (m)')
    view(0,0);
    title('2D Finite Element Model of Ice Shelf Loading')
    %zlim([-60.,10.])
    %xlim([-60 60])
    %ylim([-60,60])
    hold on;
    
    subplot(1,2,2);
    %f2 = figure(2);
    scatter(times_log,depths, 'filled');
    daystr = num2str(time,'%2.2f');
    xline(time, '-', {'day', daystr}, 'LineWidth', 1.) 
    title('Lake Depth (m)')
    xlim([0.,31.])
    ylim([0., 5.]);
    xlabel('time (days)')
    ylabel('depth (m)')
    %set(f2, 'visible', 'on');
    %hold off;
     set(gcf, 'PaperUnits', 'inches');
    x_width=10 ;y_width=5;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]);
    snum = num2str(iterator,'%04.f');
    filename = strcat('results/ice-shelf.', snum, '.png');
    %filename = strcat('animation/test.', num2str(time,'%2.2f'), '.png');
    saveas(f, filename, 'png');
    close all;
    
    %return;
    
end

function [depth] = get_depth(time)
    if time < 30
        depth = 5. * sin(time / (6. * pi));
    elseif (30 <= time) && (time < 31)
        depth = 5 - 5 * (time - 30.);
    else
        depth = 0.;
    end
end