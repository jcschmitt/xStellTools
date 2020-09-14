function pcurr = sum_cossq_s_free_j(x, ac_array)

% from VMEC2000/Sources/Initialization_Cleanup/profile_functions.f
% !      CASE ('sum_cossq_s_free_j')
% ! 20180219, Joachim Geiger

ncssq = length(ac_array) / 3;

%ncssq = 4;

if length(ac_array) < 21
    ac_array(21) = 0;
end

%bsta = 2.0;
%bend = 2.0;
      
xi = zeros(1,ncssq);
wd = ones(1,ncssq);
bsta = 2.0 * ones(1,ncssq);
bend = 2.0 * ones(1,ncssq);

for ii = 1:ncssq
    xi(ii) = ac_array(2+3*(ii-1));               %! location of maximum curr. dens.
    wd(ii) = ac_array(3+3*(ii-1));
    bsta(ii) = max([0.0, xi(ii) - wd(ii)/2]); %! interval start of curr.dens.part
    bend(ii) = min([xi(ii)+wd(ii)/2, 1.0]); %! interval end of curr.dens.part
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

    for ii =1:ncssq  % ! check the interval of x in ( bsta(i) : bend(i) ]
        this_coeff = ac_array(1 + (3 * (ii - 1) )) ;
        if ~(this_coeff == 0)
            if ( (this_x >= bsta(ii)) && (this_x <= bend(ii)) )
                amplitude = ac_array(1+3*(ii-1));
                x_offset =  ac_array(2+3*(ii-1));
                width =     ac_array(3+3*(ii-1));
                pcurr(jj) = pcurr(jj) +  ...
                    amplitude * ...
                    (cos( pi * (this_x - x_offset) / width ) ) ^2;
            end
        end
                
    end
    
    
end




