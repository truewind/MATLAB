% Calculates clear sky atmospheric emissivity based on a specified method
%
% RELEASE NOTES
%   Version 1.0 Written by Mark Raleigh (raleigh@ucar.edu), Oct 2013)
%   Version 2.0 Overhauled by Mark Raleigh (Feb 2015) to have structure inputs
%       and to correct errors in the Flerchinger formulas
%
% SYNTAX
%   e_clr = inc_LW_e_clr(M_INPUTS, M_PARAMS, M_OPTIONS)
%
% INPUTS
%   M_INPUTS = structure with the following variables (must have exact name)
%       Ta = air temperature, K or C (will assume C if mean is less than 40 C)
%       eo = vapor pressure (kPa) - optional if RH is provided and option Method_vpr>0
%       RH = fractional relative humidity
%       Qsi = incoming shortwave radiation (W m^-2)
%
%   M_PARAMS = structure with the coefficients used for the method.  If you
%   want to use the defaults, set M_PARAMS.P1_clr = []; or simply exclude the 
%   variable in the structure.  Variables include:
%       STA_Elev = elevation (m)
%       P1_clr = 1xn array of parameter values for PARAMETER 1 (method-specific) for clear sky method
%       P2_clr = 1xn array of parameter values for PARAMETER 2 (method-specific) for clear sky method
%       ...
%       Pn_clr = 1xn array of parameter values for PARAMETER n (method-specific) for clear sky method
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
%
%           %%% additional methods from Juszak and Pellicciotti (2013)
%           14 = Maykut and Church (1973)
%           15 = Konzelmann et al. (1994)
%           16 = Dilley and O'Brien (A) (1998)
%
%           %%% other methods
%           17 = Campbell and Norman (1998) as cited by Walter et al (2005)
%           18 = Long and Turner (2008) - based on Brutsaert (1975)
%           19 = Ohmura (1982) as cited by Howard and Stull 2013
%           20 = Efimova (1961) as cited by Key et al (1996)
% 
% OUTPUTS
%   e_clr = clear-sky emissivity of the atmosphere
%
% NOTES
%   The list of LW methods 1-13 is populated based on Flerchinger et al. (2009).
%   Additional empirical LW models may be added to this list.

function e_clr = inc_LW_e_clr(M_INPUTS, M_PARAMS, M_OPTIONS)

%% constants

stefan = 5.67 * (10^-8);        % Stefan-Boltzmann constant (J/s/m^2/K^4)
T_C2K = 273.15;                 % conversion from C to K
% E_L2W = 41840/86400;            % conversion from langleys/day to W/m2

%% Checks

if M_OPTIONS.Method_clr-floor(M_OPTIONS.Method_clr)~=0 || M_OPTIONS.Method_clr < 1 || M_OPTIONS.Method_clr>20
    error('Invalid M_OPTIONS.Method_clr')
end

if M_OPTIONS.Method_vpr-floor(M_OPTIONS.Method_vpr)~=0 || M_OPTIONS.Method_vpr < 0 || M_OPTIONS.Method_vpr>2
    error('Invalid M_OPTIONS.Method_vpr')
end

% Now check to see if Ta is in Celsius, and if so, convert to K
if nanmean(nanmean(M_INPUTS.Ta)) < 40
    % then assume the input was in Celsius
    M_INPUTS.Ta = M_INPUTS.Ta + T_C2K;
end

% Now check to make sure RH is in fractional
if isfield(M_INPUTS, 'RH')==1
    if nanmean(nanmean(M_INPUTS.RH)) > 1
        disp('Assuming RH was input in percent.  Converting to fractional')
        M_INPUTS.RH = M_INPUTS.RH/100;
    end
end


%% Common variables

Ta = M_INPUTS.Ta;
RH = M_INPUTS.RH;

if M_OPTIONS.Method_vpr==0
    if isfield(M_INPUTS, 'eo') ~=1
        error('eo must be provided in inputs, or include RH and specify method for computing eo')
    end
    
    eo = M_INPUTS.eo;
%     esat = M_INPUTS.esat;
elseif M_OPTIONS.Method_vpr==1
    error('this has not been coded yet')
elseif M_OPTIONS.Method_vpr==2
    % Clausius-Clapeyron e_sat in mb (hPa) from Murray 1967
    sat_vap_pressure = 6.1078.*exp((17.2693882.*(Ta-T_C2K))./(Ta-35.86));
    eo = M_INPUTS.RH .* sat_vap_pressure;
    
    % Convert eo from hPa to kPa
    eo = eo./10;
%     esat = sat_vap_pressure./10;
end



%%%%%%% at this point, eo should be in kPa and Ta should be in K, RH should
%%%%%%% be fractional

% Prata (1996) approximation for precipitable water
w = 465 .* (eo./Ta);            % note: corrected to be 465 instead of 4650 (error in Flerchinger Table 1)

%% Method-specific default parameters


if M_OPTIONS.Method_clr==1
    % 1 = Angstrom (1918)
    nparams = 3;
    P(1) = 0.83; P(2) = 0.18; P(3) = 0.67;  % note: corrected coefficient of P(3) to be 0.67 instead of 0.067 (error in Flerchinger Table 1)
elseif M_OPTIONS.Method_clr==2
    % 2 = Brunt (1932)
    nparams = 2;
    P(1) = 0.52; P(2) = 0.205;
elseif M_OPTIONS.Method_clr==3
    % 3 = Brutsaert (1975)
    nparams = 2;
    P(1) = 1.723; P(2) = (1/7);
elseif M_OPTIONS.Method_clr==4
    % 4 = Garratt (1992)
    nparams = 3;
    P(1) = 0.79; P(2) = 0.17; P(3) = 0.96;
elseif M_OPTIONS.Method_clr==5
    % 5 = Idso and Jackson (1969) (Idso-1)
    nparams = 2;
    P(1) = 0.261; P(2) = 0.00077;
elseif M_OPTIONS.Method_clr==6
    % 6 = Idso (1981) (Idso-2)
    nparams = 3;
    P(1) = 0.70; P(2) = 5.95 .* (10.^-4); P(3) = 1500;
elseif M_OPTIONS.Method_clr==7
    % 7 = Iziomon et al. (2003)
    nparams = 6;
    P(1) = 0.35; P(2) = 100; P(3) = 212;  % low land site (P(3)=elev m)
    P(4) = 0.43; P(5) = 115; P(6) = 1489; % higher elev site (P(6) = elev m)
elseif M_OPTIONS.Method_clr==8
    % 8 = Keding (1989)
    nparams = 3;
    P(1) = 0.92; P(2) = 0.7; P(3) = 1.2;
elseif M_OPTIONS.Method_clr==9
    % 9 = Niemela et al. (2001)
    nparams = 4;
    P(1) = 0.72; P(2) = 0.09; P(3) = 0.2; P(4) = 0.76;
elseif M_OPTIONS.Method_clr==10
    % 10 = Prata (1996)
    nparams = 3;
    P(1) = 1.2; P(2) = 3; P(3) = 0.5;
elseif M_OPTIONS.Method_clr==11
    % 11 = Satterlund (1979)
    nparams = 2;
    P(1) = 1.08; P(2) = 2016;
elseif M_OPTIONS.Method_clr==12
    % 12 = Swinbank (1963)
    nparams = 2;
    P(1) = (5.31 * (10^-13)); P(2) = 6;
elseif M_OPTIONS.Method_clr==13
    % 13 = Dilley and O'Brien (1998)
    nparams = 3;
    P(1) = 59.38; P(2) = 113.7; P(3) = 96.96;
elseif M_OPTIONS.Method_clr==14
    % 14 = Maykut and Church (1973)
    nparams = 1;
    P(1) = 0.7855;
elseif M_OPTIONS.Method_clr==15
    % 15 = Konzelmann et al. (1994)
    nparams = 3;
    P(1) = 0.23; P(2) = 0.484; P(3) = 1.8;
elseif M_OPTIONS.Method_clr==16
    % 16 = Dilley and O'Brien (A) (1998)
    nparams = 3;
    P(1) = 2.232; P(2) = -1.875; P(3) = 0.7356;
elseif M_OPTIONS.Method_clr==17
    % 17 = Campbell and Norman (1998) as cited by Walter et al (2005)
    nparams = 2;
    P(1) = 0.72 - 0.005*T_C2K;      % modified so Ta is in (K) instead of (C) ....should be -0.6458 after the maths...
    P(2) = 0.005;
elseif M_OPTIONS.Method_clr==18
    % 18 = Long and Turner (2008) - based on Brutsaert (1975)
    nparams = 5;
    %%% simplify the k-value determination - use one value for day, one for
    %%% night (where night/day is based on threshold in Qsi data)
    P(1) = 1.18;            % daytime k value, varies diurnally (and likely spatially) - see Fig 2 
    P(2) = 1.28;            % nighttime k value
    P(3) = (1.39e-11 + 3.36e-12 + 1.47e-11 + 4.07e-12)/4;          % "a" coefficient - average of four datasets
    P(4) = (4.8769+5.1938+4.8768+5.1421)/4;                        % "b" coefficient - average of four datasets
    P(5) = (1/7);           % brutsaert exponent
elseif M_OPTIONS.Method_clr==19
    % 19 = Ohmura (1982) as cited by Howard and Stull 2013
    nparams = 2;
    P(1) = 8.733 * 10^-3; P(2) = 0.788;
elseif M_OPTIONS.Method_clr==20
    % 20 = Efimova (1961) as cited by Key et al (1996)
    nparams = 2;
    P(1) = 0.746;
    P(2) = 0.0066*10;  % need to multiply by 10 to account for e in mb isntead of kPa
end

         


%%% check for existence of parameters and set default values if no params were specified
params_spec = zeros(1,nparams); % initialize variable to get size of each param

for j=1:nparams
    
    %%% check if this parameter is a field
    EXcheck = isfield(M_PARAMS, ['P' num2str(j) '_clr']);
    
    %%% if not, assign the default value
    if EXcheck==0
        eval(['M_PARAMS.P' num2str(j) '_clr = P(j);']);
        Px = P(j);
    else
        eval(['Px = M_PARAMS.P' num2str(j) '_clr;']);
    end
    
    
    if isarray(Px) ~=1 && numel(Px)~=1
        error('M_PARAMS.Pn_clr variables should be single values or arrays')
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
            eval(['Px = M_PARAMS.P' num2str(j) '_clr;']);
            
            if size(Px,2)==1 && size(Px,1)>1
                Px = Px';
            end
            
            eval(['M_PARAMS.P' num2str(j) '_clr = repmat(Px, size(To,1), 1);']);
        end
    end
end


%% Code

if M_OPTIONS.Method_clr==1
    % 1 = Angstrom (1918)
    e_clr = (M_PARAMS.P1_clr - M_PARAMS.P2_clr .* 10.^(-1.* M_PARAMS.P3_clr .*eo));
elseif M_OPTIONS.Method_clr==2
    % 2 = Brunt (1932)
    e_clr = (M_PARAMS.P1_clr + M_PARAMS.P2_clr .*sqrt(eo));
elseif M_OPTIONS.Method_clr==3
    % 3 = Brutsaert (1975)
    e_clr = M_PARAMS.P1_clr .* (eo./Ta).^(M_PARAMS.P2_clr);
elseif M_OPTIONS.Method_clr==4
    % 4 = Garratt (1992)
    e_clr = M_PARAMS.P1_clr - M_PARAMS.P2_clr .*exp(-1.* M_PARAMS.P3_clr .*eo);
elseif M_OPTIONS.Method_clr==5
    % 5 = Idso and Jackson (1969) (Idso-1)
    e_clr = 1 - M_PARAMS.P1_clr .*exp(-1.* M_PARAMS.P2_clr .*(Ta-T_C2K).^2);
elseif M_OPTIONS.Method_clr==6
    % 6 = Idso (1981) (Idso-2)
    e_clr = M_PARAMS.P1_clr + (M_PARAMS.P2_clr .* eo .* exp(M_PARAMS.P3_clr./Ta));
elseif M_OPTIONS.Method_clr==7
    % 7 = Iziomon et al. (2003)
    Mxz = (M_PARAMS.P4_clr-M_PARAMS.P1_clr)/(M_PARAMS.P6_clr-M_PARAMS.P3_clr);
    Myz = (M_PARAMS.P5_clr-M_PARAMS.P2_clr)/(M_PARAMS.P6_clr-M_PARAMS.P3_clr);
    X = Mxz .*(M_PARAMS.STA_Elev - M_PARAMS.P3_clr) + M_PARAMS.P1_clr;
    Y = Myz .*(M_PARAMS.STA_Elev - M_PARAMS.P3_clr) + M_PARAMS.P2_clr;
    e_clr = 1 - X.*exp(-Y.*eo./Ta);
elseif M_OPTIONS.Method_clr==8
    % 8 = Keding (1989)
    e_clr = M_PARAMS.P1_clr - M_PARAMS.P2_clr .*(10.^(-1.* M_PARAMS.P3_clr .*eo));
elseif M_OPTIONS.Method_clr==9
    % 9 = Niemela et al. (2001)
    e_clr =         M_PARAMS.P1_clr + M_PARAMS.P2_clr .*(eo- M_PARAMS.P3_clr);
    e_clr(eo<M_PARAMS.P3_clr) = M_PARAMS.P1_clr - M_PARAMS.P4_clr .* (eo(eo<M_PARAMS.P3_clr)-M_PARAMS.P3_clr);
elseif M_OPTIONS.Method_clr==10
    % 10 = Prata (1996)
    e_clr = 1 - (1+w).*exp(-1*(M_PARAMS.P1_clr + M_PARAMS.P2_clr .*w).^M_PARAMS.P3_clr);
elseif M_OPTIONS.Method_clr==11
    % 11 = Satterlund (1979)
    e_clr = M_PARAMS.P1_clr .*(1-exp(-((10.*eo).^(Ta/M_PARAMS.P2_clr))));  % note: corrected so that P2 paramter is in the exponent ^(Ta/P2), error in Flerchinger Table 1
elseif M_OPTIONS.Method_clr==12
    % 12 = Swinbank (1963)
    Lclr = M_PARAMS.P1_clr .* (Ta .^ M_PARAMS.P2_clr);
    e_clr = Lclr./(stefan.*(Ta.^4));    % effective emissivity
elseif M_OPTIONS.Method_clr==13
    % 13 = Dilley and O'Brien (1998)
    Lclr = M_PARAMS.P1_clr + M_PARAMS.P2_clr .*((Ta/T_C2K).^6) + M_PARAMS.P3_clr .*sqrt(w./2.5); % note: corrected to be 2.5 instead of 25 (error in Flerchinger Table 1)
    e_clr = Lclr./(stefan.*(Ta.^4));    % effective emissivity
elseif M_OPTIONS.Method_clr==14
    % 14 = Maykut and Church (1973)
    e_clr = M_PARAMS.P1_clr;
elseif M_OPTIONS.Method_clr==15
    % 15 = Konzelmann et al. (1994)
    e_clr = M_PARAMS.P1_clr + M_PARAMS.P2_clr .* (1000*eo./Ta).^(M_PARAMS.P3_clr);
elseif M_OPTIONS.Method_clr==16
    % 16 = Dilley and O'Brien (A) (1998)
    e_clr = 1 - exp(-1.66.*(M_PARAMS.P1_clr + M_PARAMS.P2_clr .*(Ta/T_C2K) + M_PARAMS.P3_clr .*sqrt(w./2.5)));
elseif M_OPTIONS.Method_clr==17
    % 17 = Campbell and Norman (1998) as cited by Walter et al (2005)
    e_clr = M_PARAMS.P1_clr + M_PARAMS.P2_clr.*Ta;
elseif M_OPTIONS.Method_clr==18
    % 18 = Long and Turner (2008) - based on Brutsaert (1975)
    k_array = Ta*0 + M_PARAMS.P1_clr;                % set all values to daytime k
    if isfield(M_INPUTS, 'Qsi')==1
        k_array(M_INPUTS.Qsi<50) = M_PARAMS.P2_clr;      % set night time values to night k
    end
    Ccoeff = k_array + M_PARAMS.P3_clr.*(RH.*100).^M_PARAMS.P4_clr;
    e_clr = Ccoeff .* (eo.*10./Ta).^M_PARAMS.P5_clr;      % times 10 to convert from kPa to mb
elseif M_OPTIONS.Method_clr==19
    % 19 = Ohmura (1982) as cited by Howard and Stull 2013
    e_clr = M_PARAMS.P1_clr .* Ta.^(M_PARAMS.P2_clr);
elseif M_OPTIONS.Method_clr==20
    % 20 = Efimova (1961) as cited by Key et al (1996)
    e_clr = M_PARAMS.P1_clr + M_PARAMS.P2_clr .*eo;
end

%%% apply limits and recalculate Lclr
e_clr(e_clr<0) = 0;
e_clr(e_clr>1) = 1;


