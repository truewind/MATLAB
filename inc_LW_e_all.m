% Calculates all sky atmospheric emissivity based on a specified method for
% clear sky and a specified method for cloudy skies
%
% RELEASE NOTES
%   Written by Mark Raleigh (raleigh@ucar.edu), Oct 2013)
%   Version 2.0 Overhauled by Mark Raleigh (Feb 2015) to have structure inputs
%
% SYNTAX
%   [LWdwn,e_all] = inc_LW_e_all(M_INPUTS, M_PARAMS, M_OPTIONS)
%
% INPUTS
%   M_INPUTS = structure with the following variables (must have exact name)
%       TIME = time matrix (per time_builder.m format), same number of rows as other inputs
%       Ta = air temperature, K or C (will assume C if mean is less than 40 C)
%       eo = vapor pressure (kPa) - optional if RH is provided and option Method_vpr>0
%       RH = fractional relative humidity
%       c_or_s = time series of cloud fraction (0-1 range) for Method_all from 1-8, or solar index (0-1 range) for Method_all from 9-10
%
%   M_PARAMS = structure with the coefficients used for the method.  If you
%   want to use the defaults, set M_PARAMS.P1 = []; or simply exclude the
%   variable in the structure.  Variables include:
%       STA_Elev = elevation (m)
%       P1_clr = 1xn array of parameter values for PARAMETER 1 (method-specific) for clear sky method
%       P2_clr = 1xn array of parameter values for PARAMETER 2 (method-specific) for clear sky method
%       ...
%       Pn_clr = 1xn array of parameter values for PARAMETER n (method-specific) for clear sky method
%
%       P1_all = 1xn array of parameter values for PARAMETER 1 (method-specific) for all sky method (cloud correction method)
%       P2_all = 1xn array of parameter values for PARAMETER 2 (method-specific) for all sky method (cloud correction method)
%       ...
%       Pn_all = 1xn array of parameter values for PARAMETER n (method-specific) for all sky method (cloud correction method)
%
%   M_OPTIONS = structure with options for the
%       Method_vpr = if vapor pressure (eo) is not supplied in M_INPUTS,
%       need to specify how to calculate it based on Ta and RH. Methods:
%           0 = do not calculate (eo is supplied)
%           1 = calculate with Dozier and Marks approach
%           2 = calculate with Clausius-Clapeyron e_sat in mb (hPa) from Murray 1967
%       Method_clr = enter code for clear-sky longwave method, where:
%           1 = Angstrom (1918)
%           2 = Brunt (1932)
%           3 = Brutsaert (1975)
%           4 = Garratt (1992)
%           5 = Idso and Jackson (1969) (Idso-1)
%           6 = Idso (1981) (Idso-2)
%           7 = Iziomon et al. (2003)
%           8 = Keding (1989)
%           9 = Niemela et al. (2001)
%           10 = Prata (1996)
%           11 = Satterlund (1979)
%           12 = Swinbank (1963)
%           13 = Dilley and O'Brien (1998)
%       Method_all = enter code for all-sky longwave method, where:
%           %%% Cloud cover based methods
%           1 = Brutsaert (1982)
%           2 = Iziomon et al. (2003)
%           3 = Jacobs (1978)
%           4 = Keding (1989)
%           5 = Maykut and Church (1973)
%           6 = Sugita and Brutsaert (1993)
%           7 = Unsworth and Monteith (1975)
%           8 = Kimball et al. (1982)
%           %%% Clearness/solar index based methods
%           9 = Crawford and Duchon (1999)
%           10 = Lhomme et al. (2007)
%
% OUTPUTS
%   LWdwn = incoming longwave (W m^-2)
%   e_all = all-sky emissivity of the atmosphere
%
% NOTES
%   The list of LW methods is populated based on Flerchinger et al. (2009).
%   Additional empirical LW models may be added to this list.

function [LWdwn,e_all] = inc_LW_e_all(M_INPUTS, M_PARAMS, M_OPTIONS)
% function [LWdwn,e_all] = inc_LW_e_all(Ta, RH, c_or_s, STA_Elev, Method_clr, Method_all)

%% Checks

if M_OPTIONS.Method_clr-floor(M_OPTIONS.Method_clr)~=0 || M_OPTIONS.Method_clr < 1 || M_OPTIONS.Method_clr>13
    error('Invalid M_OPTIONS.Method_clr')
end

if M_OPTIONS.Method_all-floor(M_OPTIONS.Method_all)~=0 || M_OPTIONS.Method_all < 1 || M_OPTIONS.Method_all>10
    error('Invalid M_OPTIONS.Method_all')
end

% Now check to see if Ta is in Celsius, and if so, convert to K
if nanmean(nanmean(M_INPUTS.Ta)) < 40
    % then assume the input was in Celsius
    M_INPUTS.Ta = M_INPUTS.Ta + 273.15;
end

% Now check to make sure RH is in fractional
if isfield(M_INPUTS, 'RH')==1
    if nanmean(nanmean(M_INPUTS.RH)) > 1
        disp('Assuming RH was input in percent.  Converting to fractional')
        M_INPUTS.RH = M_INPUTS.RH/100;
    end
end

%% Common variables

stefan = 5.67 * (10^-8);        % Stefan-Boltzmann constant (J/s/m^2/K^4)
Ta = M_INPUTS.Ta;

if M_OPTIONS.Method_vpr==0
    if isfield(M_INPUTS, 'eo') ~=1
        error('eo must be provided in inputs, or include RH and specify method for computing eo')
    end
    
    eo = M_INPUTS.eo;
elseif M_OPTIONS.Method_vpr==1
    error('this has not been coded yet')
elseif M_OPTIONS.Method_vpr==2
    % Clausius-Clapeyron e_sat in mb (hPa) from Murray 1967
    sat_vap_pressure = 6.1078.*exp((17.2693882.*(Ta-273.16))./(Ta-35.86));
    eo = M_INPUTS.RH .* sat_vap_pressure;
    
    % Convert eo from hPa to kPa
    eo = eo./10;
end

c_or_s = M_INPUTS.c_or_s;

%%%%%%% at this point, eo should be in kPa and Ta should be in K

%% Method-specific default parameters


if M_OPTIONS.Method_all==1
    % 1 = Brutsaert (1982)
    nparams = 1;
    P(1) = 0.22;
elseif M_OPTIONS.Method_all==2
    % 2 = Iziomon et al. (2003)
    nparams = 5;
    P(1) = 0.35; P(2) = 212; P(3) = 0.50; P(4) = 1489; P(5) = 2;
elseif M_OPTIONS.Method_all==3
    % 3 = Jacobs (1978)
    nparams = 1;
    P(1) = 0.26;
elseif M_OPTIONS.Method_all==4
    % 4 = Keding (1989)
    nparams = 2;
    P(1) = 0.153; P(2) = 2.183;
elseif M_OPTIONS.Method_all==5
    % 5 = Maykut and Church (1973)
    nparams = 2;
    P(1) = 0.22; P(2) = 2.75;
elseif M_OPTIONS.Method_all==6
    % 6 = Sugita and Brutsaert (1993)
    nparams = 2;
    P(1) = 0.0496; P(2) = 2.45;
elseif M_OPTIONS.Method_all==7
    % 7 = Unsworth and Monteith (1975)
    nparams = 2;
    P(1) = 0.84; P(2) = 0.84;
elseif M_OPTIONS.Method_all==8
    % 8 = Kimball et al. (1982)
    nparams = 9;
    P(1) = 11; P(2) = -0.6732; P(3) = 0.6240 .* (10^-2); P(4) = 0.9140 .* (10^-5);
    P(5) = 0.24; P(6) = 2.98 .* (10^-6); P(7) = 3000; P(8) = 1.4; P(9) = 0.4;
elseif M_OPTIONS.Method_all==9
    % 9 = Crawford and Duchon (1999)
    nparams = 0;
elseif M_OPTIONS.Method_all==10
    % 10 = Lhomme et al. (2007)
    nparams = 2;
    P(1) = 1.37; P(2) = 0.34;
end


%%% check for existence of parameters and set default values if no params were specified

if nparams>0
    
    params_spec = zeros(1,nparams); % initialize variable to get size of each param
    
    for j=1:nparams
        
        %%% check if this parameter is a field
        EXcheck = isfield(M_PARAMS, ['P' num2str(j) '_all']);
        
        %%% if not, assign the default value
        if EXcheck==0
            eval(['M_PARAMS.P' num2str(j) '_all = P(j);']);
            Px = P(j);
        else
            eval(['Px = M_PARAMS.P' num2str(j) '_all;']);
        end
        
        
        if isarray(Px) ~=1 && numel(Px)~=1
            error('M_PARAMS.Pn_all variables should be single values or arrays')
        end
        
        params_spec(1,j) = numel(Px);
    end
    
    
    %%% assemble parameter sets
    if nanvar(params_spec)~=0
        error('All parameter sets must have the same number of elements')
    else
        if params_spec(1,1)>1
            
            %%% need to do repeat matrix
            for j=1:nparams
                eval(['Px = M_PARAMS.P' num2str(j) '_all;']);
                
                if size(Px,2)==1 && size(Px,1)>1
                    Px = Px';
                end
                
                eval(['M_PARAMS.P' num2str(j) '_all = repmat(Px, size(To,1), 1);']);
            end
        end
    end
    
end

%% Code

%%% First, calculate clear-sky emissivity for the selected method
e_clr = inc_LW_e_clr(M_INPUTS, M_PARAMS, M_OPTIONS);

%%% double check to make sure we have e_clr within realistic limits [0 1]
e_clr(e_clr<0) = 0;
e_clr(e_clr>1) = 1;

%%% Now, calculate all-sky emissivity
if M_OPTIONS.Method_all==1
    % 1 = Brutsaert (1982)
    e_all = (1+ M_PARAMS.P1_all .*c_or_s).*e_clr;
elseif M_OPTIONS.Method_all==2
    % 2 = Iziomon et al. (2003)
    Mzz = (M_PARAMS.P3_all-M_PARAMS.P1_all)/(M_PARAMS.P4_all-M_PARAMS.P2_all);
    Z = Mzz .*(M_PARAMS.STA_Elev - M_PARAMS.P2_all) + M_PARAMS.P1_all;
    e_all = (1+Z.*c_or_s.^M_PARAMS.P5_all).*e_clr;
elseif M_OPTIONS.Method_all==3
    % 3 = Jacobs (1978)
    e_all = (1+ M_PARAMS.P1_all .*c_or_s).*e_clr;
elseif M_OPTIONS.Method_all==4
    % 4 = Keding (1989)
    e_all = (1+ M_PARAMS.P1_all .*c_or_s.^M_PARAMS.P2_all).*e_clr;
elseif M_OPTIONS.Method_all==5
    % 5 = Maykut and Church (1973)
    e_all = (1+ M_PARAMS.P1_all .*c_or_s.^M_PARAMS.P2_all).*e_clr;
elseif M_OPTIONS.Method_all==6
    % 6 = Sugita and Brutsaert (1993)
    e_all = (1+ M_PARAMS.P1_all .*c_or_s.^M_PARAMS.P2_all).*e_clr;
elseif M_OPTIONS.Method_all==7
    % 7 = Unsworth and Monteith (1975)
    e_all = (1-M_PARAMS.P1_all.*c_or_s).*e_clr + M_PARAMS.P2_all.*c_or_s;
elseif M_OPTIONS.Method_all==8
    % 8 = Kimball et al. (1982)
    Tc = Ta - M_PARAMS.P1_all; %%% no seasonal variation in Tc yet... could add in later
    f8 = M_PARAMS.P2_all + M_PARAMS.P3_all.*Tc - M_PARAMS.P4_all.* Tc.^2;
    e8z = M_PARAMS.P5_all + M_PARAMS.P6_all .* eo.^2 .* exp(M_PARAMS.P7_all./Ta);
    tau8 = 1- e8z.*(M_PARAMS.P8_all-M_PARAMS.P9_all.*e8z);
    
    Lclr = e_clr.*stefan.*Ta.^4;
    
    Ld = Lclr + tau8.*c_or_s.*f8.*stefan.*(Tc.^4);
    
    e_all = Ld./(stefan.*(Tc.^4));  % effective emissivity
    
    %%% make sure e_all is in realistic limits [0 1]
    e_all(e_all<0) = 0;
    e_all(e_all>1) = 1;
    
    %%% recalculate Ld
    LWdwn = e_all.*stefan.*(Tc.^4);        % note it is Tc not Ta for this method
    
elseif M_OPTIONS.Method_all==9
    % 9 = Crawford and Duchon (1999)
    e_all = (1-c_or_s) + c_or_s.*e_clr;
elseif M_OPTIONS.Method_all==10
    % 10 = Lhomme et al. (2007)
    e_all = (M_PARAMS.P1_all - M_PARAMS.P2_all.*c_or_s).*e_clr;
end

if M_OPTIONS.Method_all ~= 8
    %%% make sure e_all is in realistic limits [0 1]
    e_all(e_all<0) = 0;
    e_all(e_all>1) = 1;
    LWdwn = e_all.*stefan.*Ta.^4;
end
