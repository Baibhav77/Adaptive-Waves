function [A, T, R] = ATR1D( stack )
%% FUNCTION ATR1D 
% calculates absorptance, transmittance and reflectance for 1D multilayered stack
% -------------------------------------------------------------------------
%% INPUT
% -------------------------------------------------------------------------
% stack - stack dictionary, see 'set_stack.m' for details
% -------------------------------------------------------------------------
%% OUTPUT
% -------------------------------------------------------------------------
% A      - absorptance
% T      - transmittance
% R      - reflectance
% -------------------------------------------------------------------------
%% REFERENCE
% -------------------------------------------------------------------------
% See for details:
% K. Ohta and H. Ishida, 
% Matrix formalism for calculation of electric field intensity of light 
%                                       in stratified multilayered films 
% Appl. Opt. 29(13), 1952–1959 (1990); 10.1364/AO.29.001952
% 
% incoherent case is treated in:
% C. C. Katsidis and D. I. Siapkas, 
% General transfer-matrix method for optical multilayer systems with coherent, 
% partially coherent, and incoherent interference 
% Appl. Opt. 41(19), 3978–3987 (2002); ﻿10.1364/AO.41.003978
% -------------------------------------------------------------------------
%% CHECKING INPUT AND ALLOCATING WORKSPACE
% -------------------------------------------------------------------------
lam = stack{"wavelength"};
d = stack{"thick"};
theta0 = stack{"angle"};
n = stack{"nk"};
%
rs = zeros(numel(lam),numel(d)+1); rp = zeros(numel(lam),numel(d)+1);
ts = zeros(numel(lam),numel(d)+1); tp = zeros(numel(lam),numel(d)+1);
t.s = cell(numel(lam),numel(d)+1);    t.p = cell(numel(lam),numel(d)+1);
ph  = cell(numel(lam),numel(d));
%
R.s = zeros(numel(lam),1); R.p = zeros(numel(lam),1);
T.s = zeros(numel(lam),1); T.p = zeros(numel(lam),1);
I.s = cell(numel(lam),1); I.p = cell(numel(lam),1);
% -------------------------------------------------------------------------
%% CALCULATING PHASE SHIFTS FOR EACH LAYER
% -------------------------------------------------------------------------
theta = zeros( size( n, 1 ), size( n, 2 ) );                      
theta(:,1) = theta0;
% -------------------------------------------------------------------------
for i = 2:size( n, 2 )
    theta(:,i) = asin( n(:,i-1)./n(:,i).*sin( theta(:,i-1) ) );   % Snell's law
end
% -------------------------------------------------------------------------
delt = 2.*pi.*n(:,2:end-1).*...
       bsxfun( @times, bsxfun( @times, d, cos( theta(:,2:end-1) ) ), 1./lam );    
edp = exp(  1i.*delt );
edm = 1./edp;
% -------------------------------------------------------------------------
%% CALCULATING FRESNEL COEFFICIENTS FOR EACH BOUNDARY
% -------------------------------------------------------------------------
for i = 1:numel(d)+1
% reflection
% -------------------------------------------------------------------------
    rs(:,i) = ( n(:,i).*cos( theta(:,i) )   - n(:,i+1).*cos( theta(:,i+1) ) )./...
              ( n(:,i).*cos( theta(:,i) )   + n(:,i+1).*cos( theta(:,i+1) ) );
    rp(:,i) = ( n(:,i).*cos( theta(:,i+1) ) - n(:,i+1).*cos( theta(:,i) ) )./...
                 ( n(:,i).*cos( theta(:,i+1) ) + n(:,i+1).*cos( theta(:,i) ) );
% transmission
% -------------------------------------------------------------------------
    ts(:,i) =            2*n(:,i).*cos( theta(:,i) )./...
                ( n(:,i).*cos( theta(:,i) ) + n(:,i+1).*cos( theta(:,i+1) ) );
    tp(:,i) =            2*n(:,i).*cos( theta(:,i) )./...
                ( n(:,i).*cos( theta(:,i+1) ) + n(:,i+1).*cos( theta(:,i) ) );
end
% -------------------------------------------------------------------------
%% CALCULATING TRANSFER MATRIX ELEMENTS FOR EACH LAYER
% -------------------------------------------------------------------------
for i = 1:numel(lam)
    t.s{i,1} = [ 1, rs(i,1); rs(i,1), 1 ] ./ ts(i,1);
    t.p{i,1} = [ 1, rp(i,1); rp(i,1), 1 ] ./ tp(i,1);
    for j = 2:numel(d)+1
        ph{i,j-1} = [ edm(i,j-1), 0; 0, edp(i,j-1) ];
        t.s{i,j} = [ 1, rs(i,j); rs(i,j), 1 ] ./ ts(i,j);
        t.p{i,j} = [ 1, rp(i,j); rp(i,j), 1 ] ./ tp(i,j);
    end
end
% -------------------------------------------------------------------------
%% GETTING S & P TRANSMITTANCE & REFLECTANCE
% -------------------------------------------------------------------------
% setting up useful coefficients
% -------------------------------------------------------------------------
cfTs = real(n(:,end).*cos(theta(:,end))./n(:,1)./cos(theta0));
cfTp = real(conj(n(:,end)).*cos(theta(:,end))./conj(n(:,1))./cos(theta0));
% -------------------------------------------------------------------------
nincoh = stack{"nincoh"};
for i = 1:numel(lam)
    I.s{i} = eye(2); I.p{i} = eye(2);                                    % initiate the values of intensity with unity
    j = 1;                                                               % starting with the first layer
    while j <= numel(d)                                                  % looping through the whole stack
        tmp.s = t.s{i,j}; tmp.p = t.p{i,j};                              % setting up transmittance coefficients for the first layer in a coherent batch
        while nincoh(j)                                                  % while coherence is preserved, we use conventional T-matrix and track phases
            tmp.s = tmp.s*ph{i,j}*t.s{i,j+1};                            % conventional T-matrix s-polarized
            tmp.p = tmp.p*ph{i,j}*t.p{i,j+1};                            % conventional T-matrix p-polarized
            j = j + 1;                                                   % next layer
            if j > numel(d)                                              % if we reached the end of the stack
                I.s{i} = I.s{i}*abs(tmp.s).^2;                           % final intensity s-pol
                I.p{i} = I.p{i}*abs(tmp.p).^2;                           % final intensity p-pol
                break;                                                   % breaking the loop
            end                                                          % if no incoherence involved, this ends up as a conventional T-matrix
        end
        if j <= numel(d) && ~nincoh(j)                                   % if there is an incoherent layer, we continue
            I.s{i} = I.s{i}*abs(tmp.s).^2*abs(ph{i,j}).^2;               % independent on the location of incoherent layer,
            I.p{i} = I.p{i}*abs(tmp.p).^2*abs(ph{i,j}).^2;               % we calculate intensity and square of the phase
            j = j + 1;                                                   % next layer
            if j > numel(d)                                              % if we reached the exit medium,
                I.s{i} = I.s{i}*abs(t.s{i,j}).^2;                        % get transmittance
                I.p{i} = I.p{i}*abs(t.p{i,j}).^2;                          
                break;                                                   % and exit the loop
            elseif ~nincoh(j)                                            % if it's not the end and yet another incoherent layer
                I.s{i} = I.s{i}*abs(t.s{i,j}).^2*abs(ph{i,j}).^2;        % propagate intensity and square of the phase
                I.p{i} = I.p{i}*abs(t.p{i,j}).^2*abs(ph{i,j}).^2;
                j = j + 1;
                    if j > numel(d)                                              % if we reached the exit medium,
                        I.s{i} = I.s{i}*abs(t.s{i,j}).^2;                        % get transmittance
                        I.p{i} = I.p{i}*abs(t.p{i,j}).^2;                          
                        break;% next layer
                    end
            end
        end
    end
    R.s(i) = I.s{i}(2,1)/I.s{i}(1,1);                                    % collect reflectance
    R.p(i) = I.p{i}(2,1)/I.p{i}(1,1);
    T.s(i) = cfTs(i)./I.s{i}(1,1);                                       % and transmittance
    T.p(i) = cfTp(i)./I.p{i}(1,1);
end
% -------------------------------------------------------------------------
%% CONSTRUCTING OUTPUT
% -------------------------------------------------------------------------
T.sp = (T.s + T.p)/2; R.sp = (R.s + R.p)/2;
A.s = 1 - T.s - R.s; A.p = 1 - T.p - R.p; A.sp = 1 - T.sp - R.sp;
% -------------------------------------------------------------------------
end
% -------------------------------------------------------------------------

%%% BELOW IS TO BE INCORPORATED
% -------------------------------------------------------------------------
%% CALCULATING PRODUCT OF TRANSFER MATRICES
% -------------------------------------------------------------------------
% for i = 1:numel(lam)
% %     S.s{1,i} = t.s{1,i}; S.p{1,i} = t.p{1,i};
% %     for j = 2:numel(d)+1
% %         S.s{j,i} = S.s{j-1,i}*ph{j-1,i}*t.s{j,i};
% %         S.p{j,i} = S.p{j-1,i}*ph{j-1,i}*t.p{j,i};
% %     end
% 
%     for j = 1:numel(d)+1
%         if j == 1 
%             S.s{j,i} = t.s{1,i}; S.p{j,i} = t.p{1,i};
%         else
%             S.s{j,i} = eye(2); S.p{j,i} = eye(2);
%         end
% 
%         for jj = j:numel(d)
%             S.s{j,i} = S.s{j,i}*ph{jj,i}*t.s{jj+1,i};
%             S.p{j,i} = S.p{j,i}*ph{jj,i}*t.p{jj+1,i};
%         end
%     end
% 
%     for j = 1:numel(d)+2
%         if j ==numel(d)+2
%             Sfwd.s(j,i) =   real(n(j,i))*abs(1/S.s{1,i}(1,1)).^2*...
%                 cos(theta(j,i))/cos(theta0);
%             Sbck.s(j,i) =   0;
%             Sfwd.p(j,i) = 1/real(n(j,i))*abs(1/S.p{1,i}(1,1)).^2*...
%                 cos(theta(j,i))/cos(theta0);
%             Sbck.p(j,i) = 0;
%         else
%         Sfwd.s(j,i) =   real(n(j,i))*abs(S.s{j,i}(1,1)/S.s{1,i}(1,1)).^2*...
%             cos(theta(j,i))/cos(theta0);
%         Sbck.s(j,i) =   real(n(j,i))*abs(S.s{j,i}(2,1)/S.s{1,i}(1,1)).^2*...
%             cos(theta(j,i))/cos(theta0);
%         Sfwd.p(j,i) = 1/real(n(j,i))*abs(S.p{j,i}(1,1)/S.p{1,i}(1,1)).^2*...
%             cos(theta(j,i))/cos(theta0);
%         Sbck.p(j,i) = 1/real(n(j,i))*abs(S.p{j,i}(2,1)/S.p{1,i}(1,1)).^2*...
%             cos(theta(j,i))/cos(theta0);
%         end
%     end
% 
%     for j = 1:numel(d)
%         Aabs.s(j,i) = (Sfwd.s(j,i) - Sbck.s(j,i)) - (Sfwd.s(j+2,i) - Sbck.s(j+2,i));
%         Aabs.p(j,i) = (Sfwd.p(j,i) - Sbck.p(j,i)) - (Sfwd.p(j+2,i) - Sbck.p(j+2,i));
%         Aabs.sp(j,i) = (Aabs.s(j,i) + Aabs.p(j,i))/2;
%     end
% 
%     % total reflection and transmission for fields, so squared
%     R.s(i) = abs(S.s{1,i}(2,1)/S.s{1,i}(1,1)).^2;
%     R.p(i) = abs(S.p{1,i}(2,1)/S.p{1,i}(1,1)).^2;
%     T.s(i) = Sfwd.s(end,i);
%     T.p(i) = Sfwd.p(end,i);
% end