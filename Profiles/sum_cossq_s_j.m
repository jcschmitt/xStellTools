function pcurr = sum_cossq_sqrts_j(x, ac_array)

% from VMEC2000/Sources/Initialization_Cleanup/profile_functions.f
% !      CASE ('sum_cossq_sqrts')
% ! 20180219, Joachim Geiger
% ! This is the implementation for using sqrt(s) as coordinate.
% ! The same idea as for sum_cos_sq.
% ! The analytical integration give different integrals.
% !  Sum of ncssq cos**2 terms put at different radial locations
% !  and with a windowing around the maxima as current density:
% !  ac(0) holds ncssq
% !     ac(1)*(H(x)-H(x-dx))*(cos(pi*(x-xi(1))/(2*dx))**2 
% !    +ac(2)*(H(x-(xi(2)-dx))-H(x-(xi(2)+dx)))*cos(pi*(x-xi(2))/(2*dx))**2
% !    + ... 
% !    +ac(ncssq)*(H(x-(xi(ncssq)-dx))-H(x-(xi(ncssq))))*cos(pi*(x-xi(ncssq))/(2*dx))**2
% !  windowing is done with the Heavyside-functions producing the
% !  Boxcar Function: H(x-a)-H(x-b) is, for a<b, 1 for a<x<b and zero otherwise.
% !  use of variables:
% !     half of interval width : dx=1/(ncssq-1.0)   .
% !     location of ith: xi=(i-1)*dx = (i-1)/(ncssq-1)
% !     example: ncssq=5, dx=1/4.0=0.25, locations:0.0, 0.25, 0.5, ...
% !       Note the locations are with respect to sqrt(s), the s-value is given in brackets
% !     interval for first  wave:  0.00 to 0.25, maximum at 0.00  (0.0)
% !     interval for second wave:  0.00 to 0.50, maximum at 0.25  (0.0625)
% !     interval for third  wave:  0.25 to 0.75, maximum at 0.50  (0.25)
% !     interval for fourth wave:  0.50 to 1.00, maximum at 0.75  (0.5625)
% !     interval for fifth  wave:  0.75 to 1.00, maximum at 1.00  (1.0)
% !  Analytic integration now wrt sqrt(s)=x results in a more complicated sum of
% !  sine-functions and cosine-functions with and a linear increase:
% !  int(x(i)-dy to x) cos(pi*(t-xi(i))/(2*dx))**2 x dx =
% !       =  0.25*(x-(x(i))+dx**2/(2*pi**2)+
% !        + (dx/(2*pi))*x*sin(pi(x-x(i))/dx)) +
% !        + (dx**2/(2*pi**2))*x*cos(pi(x-x(i))/dx))
% !  For the very first part integration start at mid-interval,
% !  i.e. x=0, x(i)=0:
% !  int(0 to x) cos(pi*(t)/(2*dx))**2 x dx =
% !       =0.25*x+(dx/(2*pi)) * x  * sin(pi(x)/dx)
% !              +(dx**2/(2*pi**2))*(cos(pi(x)/dx)-1)
% !  When adding up fully integrated intervals the first needs
% !  also to be accounted for differently.
% !


%      ncssq=int(ac(0))

% matlab arrays start at one.  the one in vmec comes from a namelist and
% the fortran array has its index base set to '0'.
ncssq = (ac_array(1));

%       if (ncssq .lt. 1 .or. ncssq .gt. 20) then
%         write(6,*) "Number of coeff.s for sum_cos_sq :",  
%      &             "1<= sum_cos_sq <= 21!"
%         stop 'Check input!'
%       endif
if ( (ncssq < 1) || (ncssq > 20) ) 
    disp('<----Number of coeffs for sum_cos_sq :1<= sum_cos_sq <= 21!');
    disp('<----Check input!')
    pcurr = NaN;
    return
end

%       sqx=sqrt(x)
%       delx=1.0_dp/(ncssq-1.0_dp)
%       delxsq=delx*delx
%       pisq=pi*pi
%       bsta = 2.0_dp
%       bend = 2.0_dp

% Instead of 'sqx', I will use 'this_x' and 'this_sqx'.
%sqx = sqrt(x);
delx = 1.0 / (ncssq-1.0);
delxsq = delx * delx;
pisq = pi * pi;
%bsta = 2.0;
%bend = 2.0;
      
%       do i=1,ncssq
%         xi(i)=(i-1)*delx               ! location of maximum curr. dens.
%         bsta(i)=max(0.0_dp,xi(i)-delx) ! interval start of curr.dens.part
%         bend(i)=min(xi(i)+delx,1.0_dp) ! interval end of curr.dens.part
%       enddo
xi = zeros(1,ncssq);
bsta = 2.0 * ones(1,ncssq);
bend = 2.0 * ones(1,ncssq);

for ii = 1:ncssq
    xi(ii) = (ii-1) * delx;               %! location of maximum curr. dens.
    bsta(ii) = max([0.0, xi(ii) - delx]); %! interval start of curr.dens.part
    bend(ii) = min([xi(ii)+delx, 1.0]); %! interval end of curr.dens.part
end

      
%       ni=1           ! first interval initially.
% Moved to inside the loop
%ni = 1;      %     ! first interval initially.

%       do i=1,ncssq   ! check the interval of x in ( bsta(i) : bend(i) ]
%         if( sqx .le. bend(i) .and. sqx .gt. bsta(i))then
%           ni=i
%         endif
%       enddo

% loop over each x and evaluate
pcurr = zeros(size(x));

for jj = 1:length(x)
    this_x = x(jj);
    %this_sqx = sqrt(this_x);
    ni = 1;      %     ! first interval initially.
    
    for ii =1:ncssq  % ! check the interval of x in ( bsta(i) : bend(i) ]
        if ( (this_x <= bend(ii)) && (this_x > bsta(ii)) )
            ni = ii;
        end
    end
    
    
    %       pcurr=0.0
    pcurr(jj) = 0.0;
    
    
    %       do i=1,ni
    %         if( sqx .gt. bsta(i) .and. sqx .le. bend(i) ) then
    %           if(i .eq. 1) then   ! the first has no linear term.
    %             pcurr=pcurr + ac(i)*(0.25_dp*x+                                    &
    %      &          (delx/(2*pi)*sqx*sin(pi*sqx/delx))+                            &
    %      &          delxsq/(2*pisq)*(cos(pi*sqx/delx)-1.0_dp))
    %           else
    %             pcurr=pcurr + ac(i)*(0.25_dp*x+                                    &
    %      &          delx/(2*pi)*sqx*sin(pi*(sqx-xi(i))/delx)+                      &
    %      &          delxsq/(2*pisq)*cos(pi*(sqx-xi(i))/delx)-                      &
    %      &          0.25_dp*bsta(i)**2 + delxsq/(2*pisq))
    %           endif
    %         else
    %           if(i .eq. 1) then   !the first is only a half interval
    % !           if(x .gt. 0.0_dp )
    % !    &          pcurr=pcurr+ac(i)*(0.25_dp-1/pisq)*delxsq
    %             pcurr=pcurr+ac(i)*(0.25_dp-1/pisq)*delxsq
    %           else
    %             pcurr=pcurr+ac(i)*delx*xi(i)
    %           endif
    %         endif
    %       enddo


    % In region ni == 1 (axis), only the first region applies
    % In region ni = ncssq (LCFS), only the last region applies
    % In regions 1 < ni < ncssq, regions #(ni-1) and #(ni) both apply
    % Note that the coefficient are located in ac_array starting at
    % index=2, and the last index is ncssq+1
    if (ni == 1) %   !the first is only a half interval
        pcurr(jj) = 2*ac_array(ni+1);
    elseif (ni == 2) %   !the first is only a half interval
        pcurr(jj) = 2*ac_array(ni) * ...
                       (cos(pi*(this_x - xi(ni-1)) / (2*delx)))^2 + ...
                    ac_array(ni+1) * ...
                        (cos(pi*(this_x - xi(ni)) / (2*delx)))^2;
%     elseif (ni == (ncssq-1))
%         pcurr(jj) = ac_array(ni) * ...
%                        (cos(pi*(this_sqx - xi(ni-1)) / (2*delx)))^2 + ...
%                     2 * ac_array(ni+1) * ...
%                         (cos(pi*(this_sqx - xi(ni)) / (2*delx)))^2;
    elseif (ni == ncssq)
        pcurr(jj) = ac_array(ni) * ...
                       (cos(pi*(this_x - xi(ni-1)) / (2*delx)))^2 + ...
                    2 * ac_array(ni+1) * ...
                        (cos(pi*(this_x - xi(ni)) / (2*delx)))^2;
%        pcurr(jj) = 2*ac_array(ni+1);        
    else
        pcurr(jj) = ac_array(ni) * ...
                       (cos(pi*(this_x - xi(ni-1)) / (2*delx)))^2 + ...
                    ac_array(ni+1) * ...
                        (cos(pi*(this_x - xi(ni)) / (2*delx)))^2;
    end
    
    
end




