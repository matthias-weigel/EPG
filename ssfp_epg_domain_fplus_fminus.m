function [F0_vector_out,Xi_F_out,Xi_Z_out] = ssfp_epg_domain_fplus_fminus (N_in,alpha_in,TR_in,T1_in,T2_in)

% [F0_vector,Xi_F,Xi_Z] = ssfp_epg_domain_fplus_fminus (N,alpha,TR,T1,T2)
% 
% Calculates the Extended Phase Graph (EPG) for some variants of gradient echo (GE) / steady state free precession (SSFP) sequences
% The particular type of the (idealized) GE/SSFP is specified in the hard coded settings below
% This code uses the Fourier based EPG domains F+(+k), F-(-k), Z(+k), see the mentioned review paper below
% 
% IN :  N    : Number of Repetitions of RF flip angles (= number of RF excitation pulses for one given flip angle alpha)
%       alpha: Constant flip angle or flip angle array of RF excitation pulses in deg to be repeated N times
%       TR   : Repetition time, same unit as T1,T2
%       T1   : Longitudinal relaxation time, same unit as TR,T2
%       T2   : Transverse relaxation time, same unit as TR,T1
% 
% OUT:  F0_vector: Vector of resulting F0 states ("echo intensities")
%       Xi_F     : State evolution matrix of all transverse F states (immediately after the RF pulse in each case)
%       Xi_Z     : State evolution matrix of all longitudinal Z states (immediately after the RF pulse in each case)
% 
% 
% WRITTEN IN 2014 by MATTHIAS WEIGEL         (epg@matthias-weigel.net)
% Last modified in 01/2015   (Release Version 2.3)
% 
% Current affiliation: Radiological Physics, University of Basel Hospital, Basel, Switzerland
% Past    affiliation: Medical Physics, University Medical Center Freiburg, Freiburg, Germany
% Past^2  affiliation: Biophysics, University of Wuerzburg, Wuerzburg, Germany
% 
% This code resides at "http://epg.matthias-weigel.net"
% The code is based on the depiction and discussion of Extended Phase Graphs in the following publication ("EPG-R"):
% 
% Weigel M. J Magn Reson Imaging 2014; doi: 10.1002/jmri.24619. Extended Phase Graphs: Dephasing, RF Pulses, and Echoes - Pure and Simple. 
% 
% Studying and using this code means to acknowledge Matthias Weigel's months of cursing and weeping ...  ;-)
% ... by citing the above mentioned review paper. Thank you :-) 
% 
% 
% Further comments in regard to the code:
% - For transverse magnetization, this code works in both the F+ and F- domain   --   see EPG-R, particularly the comments in the "software coding section"
%   Therefore, both the F+(0) and F-(0) state are assessed (Omega matrix), however, only the F+(0) state is stored in the state evolution matrix Xi (F-(0) is superfluous)
%   One could see it as the "natural approach in regard to the T-matrix"   --   or not ;-)
% - This code is only roughly optimized for speed or for efficiency: a dedicated code (particularly in C++) would be much faster; and also quite shorter without the "gadgets"
% - This code does not consider any TE effects
% - This code automatically uses an initial alpha/2 pulse for alternating flip angles (you can change that, of course, like everything ...)
% - This code was roughly validated by comparing selected examples with ...
%           ... a "dedicated Bloch simulation" for FISP, FLASH, TRUFI
%           ... the software "ssfp_epg_domain_fplus_alone.m" for FISP, FLASH, TRUFI


% Hard coded settings to modify the type of the simulated SSFP sequence - Typical seetings:
% ---------------------------------------------------------------------------------------------
% FISP / GRASS / FAST    : unbalanced = 1; RF_phase_mode =  0; RF_phase_inc = <does not matter>
% FLASH / RF spoiled GE  : unbalanced = 1; RF_phase_mode =  1; RF_phase_inc = 50;
% TrueFISP, balanced SSFP: unbalanced = 0; RF_phase_mode = -1; RF_phase_inc = <does not matter>
% ---------------------------------------------------------------------------------------------
unbalanced    =  1;  %  1 = Unbalanced readout gradient, i.e., ideal gradient moment spoiling 
                     %  0 = Balanced readout gradient 
                    
RF_phase_mode =  0;  %  0 = All RF pulses have the same phase phi=0deg, i.e., rotation around x-axis
                     % +1 = RF spoiling with increment as specified below, starting with phi=0deg
                     % -1 = Alternating RF flip angles, i.e., alternating rotation around +x and -x
                    
RF_phase_inc  = 50;  % Used linear RF phase increment ONLY if RF_phase_mode = +1  ;  A useful value is 50deg



% Initialization of various parameters regarding EPG matrix operators, output space, ...
% -------------------------------------------------------------------------------------------

% Determine effective alpha and phi arrays and dedicated length parameters: SSFP variant settings
alpha_in = alpha_in/180.0*pi;                       % Convert flip angle(array) into [rad] 
fa       = alpha_in;                                % "alpha" is already a function in Matlab: not good

for pn = 2:N_in                                     % The given flip angle (array) is repeated N_in times
    fa = [fa,alpha_in];   %#ok<AGROW>
end

N   = length(fa);
Nm1 = N-1;                                     
Np1 = N+1;  

switch (RF_phase_mode)
    case -1                                         % Realize alternating RF flip angles 
        fa(2:2:end) = -fa(2:2:end);
        phi (1:N)    = 0;
        if (N>1)                                    % More than one alternating flip angle ?  Use initial fa/2 pulse !
            fa(1)=fa(1)/2;
        end
        
    case +1                                         % Realize RF spoiling with dedicated RF phase angle
        pn  = 1:N;
        phi = (pn-1).*pn/2.0 * RF_phase_inc / 180.0 * pi;        
        
    otherwise                                       % Realize "normal" flip angles with phi=0
        phi(1:N) = 0;
end
    

% Define relaxation operator E elements
E1 = exp(-TR_in/T1_in);                             % Eq.[23] in EPG-R 
E2 = exp(-TR_in/T2_in);                             % Eq.[23] in EPG-R


% Generate state matrices Omega before and after RF: Eq.[26] in EPG-R
Omega_preRF  = zeros(3,Np1);
Omega_postRF = zeros(3,N);


% Generate state evolution matrices Xi_F and Xi_Z: Eqs.[27] in EPG-R
% Here, they will (only) contain all post-RF states == output variables
Xi_F_out = zeros(2*N-1,N);                          % Xi_F with F+ and F- states, Eq.[27a] in EPG-R               
Xi_Z_out = zeros(N,N);                              % Xi_Z with Z states, upper half of Eq.[27b] in EPG-R


% Starting with equilibrium magnetization M0=1 for t<0   -   you could change that !
Omega_preRF(3,1) = 1;    



% Starting the calculation of EPG states over time - "flow of magnetization in the EPG"
% -------------------------------------------------------------------------------------------
for pn = 1:N                                                    % Loop over RF pulse #pn
       
    % Generate T matrix operator elements (RF pulse representation):  Eq.[15] or Eq.[18] in EPG-R
    % Since this generic SSFP example may change RF flip angle and/or phase, the T matrix has to be updated within the RF pulse loop:
    T(1,1) =                            cos(fa(pn)/2)^2;   
    T(1,2) =         exp(+2i*phi(pn)) * sin(fa(pn)/2)^2;
    T(1,3) = -1.0i * exp(+1i*phi(pn)) * sin(fa(pn)  )  ;
    T(2,1) =         exp(-2i*phi(pn)) * sin(fa(pn)/2)^2;
    T(2,2) =                            cos(fa(pn)/2)^2;
    T(2,3) = +1.0i * exp(-1i*phi(pn)) * sin(fa(pn)  )  ;
    T(3,1) = -0.5i * exp(-1i*phi(pn)) * sin(fa(pn)  )  ;
    T(3,2) = +0.5i * exp(+1i*phi(pn)) * sin(fa(pn)  )  ;           
    T(3,3) =                            cos(fa(pn)  )  ;
    
    
    % In the following, further loops over the states' dephasing k will be needed to realize operators T, E, and S
    % --> Note that we deal with integral k units here, see EPG-R
    % Intead of "for loops over k", Matlab's "indexing capability" will be used, which is much faster for large values
    % To improve efficiency, only states are calculated that could be generated at all in the EPG before
    k = 1:pn;                                                   % "k loop index" with limit pn
    
    
    % T matrix operator: RF pulse acting, mixing of F+, F-, and Z states
    % Expand T matrix relations from Eq.[15] or Eq.[18] in EPG-R
    Omega_postRF(:,k) = T * Omega_preRF(:,k);                   % matrix product !  faster than explicit writing out

    
    % Store these post-RF states of the current Omega state matrix in the Xi state evolution matrices
    Xi_F_out(Np1-k,pn)        =      Omega_postRF(1,k);
    Xi_F_out(Nm1+k(2:end),pn) = conj(Omega_postRF(2,k(2:end))); % this is a matter of definition, EPG-R defines Xi_F as a coherent set fully in the F+ domain, you could also omit the 'conj' ...
    Xi_Z_out(Np1-k,pn)        =      Omega_postRF(3,k);

    
    % E matrix operator: Experienced relaxation from the states until the next TR
    % Expand E matrix relations from Eqs.[23] and [24] in EPG-R 
    Omega_preRF(1:2,k)      = E2 * Omega_postRF(1:2,k);         % T2 relaxation for F+ and F-
    Omega_preRF(3,k(2:end)) = E1 * Omega_postRF(3,k(2:end));    % T1 relaxation for Z with k>0
    Omega_preRF(3,1)        = E1 * Omega_postRF(3,1) + 1 - E1;  % T1 recovery for Z with k=0 
    
        
    % S operator: Further dephasing / shifting of F+ and F- states to the next TR
    % With integral k units (see EPG-R):  Delta(k) = +1 (unbalanced, up shifting) , Delta(k) = 0 (balanced, skip S operator)     
    if (unbalanced == 1)      
        Omega_preRF(1,k+1) = Omega_preRF(1,k);                  % dephase F+ (this expression only works with indexing !! no overwriting then) 
        Omega_preRF(2,k)   = Omega_preRF(2,k+1);                % dephase F- 
        Omega_preRF(1,1)   = conj(Omega_preRF(2,1));            % generate conjugate pendant F+0 from F-0, see EPG-R
    end
  
end % of for loop over pn


% Output: "make nice zeros"
% Erase some float point accuracy errors
Xi_F_out(abs(Xi_F_out)<eps*1e3) = 0;
Xi_Z_out(abs(Xi_Z_out)<eps*1e3) = 0;


% Output: define "echoes" separately
F0_vector_out = Xi_F_out(N,:);

end
% bye, bye ... :-P
