function results = ReelleMatKonst
%USAGE: 
%Tatt direkte fra Lohne's appendiks C og tilpasset egne data

DD = lDisc();
results(size(DD,2)) = struct();
%==============
%   INNDATA
%==============

f1_1 = [DD.f1_R1];
f1_2 = [DD.f1_R2];
f2_1 = [DD.f2_R1];
f1T = [DD.f1_TE1];
f2T = [DD.f2_TE1];
t = [DD.t];
a = [DD.a];
rho = [DD.rho];
DT = 2*a/t;
C0_T = [DD.B_omega];
tan_d = 1./[DD.Q_e];
e31 = 1.11; %Fra Ferroperm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Iterativt finne sigmaP og vP %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n1_1 = 2.12; %Dette er kun en initialverdi i think
n1_2 = n1_1.*f1_2./f1_1;
sigma1_p = 1 - n1_1.*besselj(0,n1_1)./besselj(1,n1_1);
sigma2_p = 1 - n1_2.*besselj(0,n1_2)./besselj(1,n1_2);
x = 0.275; %initialverdi for hva n1_1 endres med

while abs(sigma1_p - sigma2_p) > 1E-8
    if sigma1_p < sigma2_p
        n1_1 = n1_1 - x;
    else
        n1_1 = n1_1 + x;
    end
    
    n1_2 = n1_1.*f1_2./f1_1;
    sigma1_p = 1 - n1_1.*besselj(0,n1_1)./besselj(1,n1_1);
    sigma2_p = 1 - n1_2.*besselj(0,n1_2)./besselj(1,n1_2);
    x = x/2;
end
sigma_p = sigma1_p;
v_p = (2*pi.*f1_1.*a)./n1_1;
n2_1 = (2*pi.*f2_1.*a)./v_p;
Ono_n2_1 = (n2_1.*besselj(0,n2_1))./besselj(1,n2_1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Utregninger for radiell mode %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kp = -sqrt((Ono_n2_1 + sigma_p - 1)./(Ono_n2_1 - 2));
k31 = kp.*sqrt(0.5.*(1 - sigma_p));
eps33_T = (C0_T.*t)./(pi.*a.^2);
eps33_p = eps33_T.*(1 - kp.^2);
s11_E = n1_1.^2./(rho.*(2*pi.*f1_1.*a).^2.*(1 - sigma_p.^2));
s12_E = - (sigma_p.*s11_E);
d31 = k31.*sqrt(eps33_T.*s11_E);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Utregninger for tykkels mode %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v_D = 2.*t.*f2T; %benytter kun tykkelses til keramikken
c33_D = (v_D).^2.*rho;
kt_2 = (( pi.*f1T)./(2.*f2T)).*tan(pi.*(f2T - f1T)./(2.*f2T));
kt = sqrt(kt_2);
c33_E = c33_D.*(1 - kt_2); %regnes ut på nytt under
eps33_S = eps33_T.*(1 - kt_2).*(1 - kp.^2);
e33 = sqrt(kt_2.*(eps33_S.*c33_D));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Materialkonstanter til FEMP  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s13_E = (d31 - e31.*(s11_E + s12_E) )./(e33);
s33_E = (eps33_T - eps33_S - 2.*(e31.^2).*(s11_E + s12_E) - 4.*(e31.*e33.*s13_E) )./(e33.^2);
d33 = 2.*e31.*s13_E + e33.*s33_E;
c11_E = (s11_E.*s33_E - (s13_E.^2) )./( (s11_E - s12_E).* (s33_E.*(s11_E + s12_E) - 2.*(s13_E.^2) ) );
c12_E = ((-s12_E).*s33_E + (s13_E.^2) )./( (s11_E - s12_E).*(s33_E.*(s11_E + s12_E) - 2.*(s13_E.^2)));
c13_E = (-s13_E)./(s33_E.*(s11_E + s12_E) - 2.*(s13_E.^2));
c33_E = (s11_E + s12_E)./(s33_E.*(s11_E + s12_E) - 2.*(s13_E.^2) );
c66_E = (c11_E - c12_E)./2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Lagrer materialkonstanter    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c11_E = num2cell(c11_E);
[results.c11_E ] = deal(c11_E {:});
c12_E = num2cell(c12_E);
[results.c12_E ] = deal(c12_E {:});
c13_E = num2cell(c13_E);
[results.c13_E ] = deal(c13_E {:});
c33_E = num2cell(c33_E);
[results.c33_E ] = deal( c33_E{:});
c66_E = num2cell(c66_E);
[results.c66_E ] = deal(c66_E {:});

e33 = num2cell(e33);
[results.e33 ] = deal(e33 {:});
eps33_S = num2cell(eps33_S);
[results.eps33_S ] = deal(eps33_S {:});

d33 = num2cell(d33);
[results.d33 ] = deal(d33 {:});
kt = num2cell(kt);
[results.kt ] = deal(kt {:});
eps33_p = num2cell(eps33_p);
[results.eps33_p ] = deal(eps33_p {:});