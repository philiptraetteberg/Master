function results = KomplekseMatKonst()
%USAGE: Taken from Appendiks C Lohne

DD = lDisc();
results(size(DD,2)) = struct();
%==============
%   INNDATA
%==============

t = [DD.t];
a = [DD.a];
%d = 2.*a;
rho = [DD.rho];
C0_T = [DD.B_omega]; %relatert til B / omega, finn ut hvordan
Qe = [DD.Q_e];
f1_1 = [DD.f1_R1];
BW_f1_1 = [DD.BW_f1_R1];
f1_2 = [DD.f1_R2];
BW_f1_2 = [DD.BW_f1_R2];
f2_1 = [DD.f2_R1];
BW_f2_1 = [DD.BW_f2_R1];
f1T = [DD.f1_TE1];
BW_f1T = [DD.BW_f1_TE1];
f2T = [DD.f2_TE1];
BW_f2T = [DD.BW_f2_TE1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%HER ANTAS E31%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e31 = 1.11*(1-1i/166);
%c44 = 'a';
%e15 = 'b';
%eps11_S = 'c';

%utregning av radielle
KOMP_f1_1 = f1_1.*(1-1i.*((BW_f1_1)./f1_1)).^(-0.5);
KOMP_f1_2 = f1_2.*(1-1i.*((BW_f1_2)./f1_2)).^(-0.5);
KOMP_f2_1 = f2_1.*(1-1i.*((BW_f2_1)./f2_1)).^(-0.5);

rs = KOMP_f1_2./KOMP_f1_1;

%Gjør rs reell hvis image delen er liten

for i = 1:size(rs)
    if abs(imag(rs(i))) <= (2*sqrt(2)*1e-3)/(sqrt((real(KOMP_f1_1)/imag(KOMP_f1_1)) ))
        rs(i) = real(rs(i));
    end    
end

n1_1 = 11.2924 - 7.63859.*rs + 2.13559.*(rs.^2) - 0.215782.*(rs.^3);
sigma_p = 97.527023 - 126.91730.*rs + 63.400384.*(rs.^2) - 14.34044.*(rs.^3) + 1.2312109.*(rs.^4);

c11_p = (2.*pi.*KOMP_f1_1.*a).^2.*rho./(n1_1.^2);
s11_E = 1./(c11_p.*(1-sigma_p.^2));
s12_E = - sigma_p.*s11_E;

n2_1 = (KOMP_f2_1.*n1_1)./KOMP_f1_1;
Ono_n2_1 = n2_1.*besselj(0,n2_1)./besselj(1,n2_1);
k_p2 = ((1- sigma_p - Ono_n2_1)./2); %this is (k^p)^2
eps33_T = (( C0_T.*t)./(pi.*a.^2)).*(1- (1i./Qe));
eps33_p = eps33_T./(1+(2.*(k_p2))./(1+sigma_p) );
kp_2 = 1 - (eps33_p./eps33_T); %this is k_p^2
k31 = -sqrt(kp_2.*0.5.*(1-sigma_p)); %negative sign here implies root of kp_2 is negative.
d31 = k31.*sqrt(eps33_T.*s11_E); %the thing above is done such that d31 is negative, which coincides with ferroperms value


%tykkelsesting
KOMP_f1_T = f1T.*(1-1i.*(BW_f1T./f1T)).^(-0.5);
KOMP_f2_T = f2T.*(1-1i.*(BW_f2T./f2T)).^(-0.5);

kt = sqrt( (pi.*KOMP_f1_T)./(2.*KOMP_f2_T).*tan( (pi.*(KOMP_f2_T - KOMP_f1_T) )./(2.*KOMP_f2_T) ));
c33_D = 4.*rho.*(t.*KOMP_f1_T).^2;
c33_E = c33_D.*(1-kt.^2); %also calculated below, gives the same results
eps33_S = eps33_T.*(1-kt.^2).*(1-kp_2);
e33 = sqrt(kt.^2.*(eps33_S.*c33_D) );

%mat konst for femp
s13_E = (d31 - e31.*(s11_E + s12_E) ) ./e33;
s33_E = (eps33_T - eps33_S - 2.*(e31.^2).*(s11_E + s12_E) - 4.*(e31.*e33.*s13_E))./(e33.^2);
d33 = 2.*e31.*s13_E + e33.*s33_E;
c11_E = (s11_E.*s33_E - (s13_E.^2))./((s11_E - s12_E).*(s33_E.*(s11_E + s12_E) - 2.*(s13_E.^2) ));
c12_E = (( - s12_E).*s33_E + (s13_E.^2))./((s11_E - s12_E).*(s33_E.*(s11_E + s12_E) - 2.*(s13_E.^2)));
c13_E = (-s13_E)./(s33_E.*(s11_E+s12_E)-2.*(s13_E.^2) );
c33_E = (s11_E + s12_E)./(s33_E.*(s11_E + s12_E) - 2.*(s13_E.^2) );
c66_E = (c11_E - c12_E)./2;

%Q-verider
Qeps33_S = real(eps33_S)./imag(eps33_S);
%Qeps11_S = real(eps11_S)./imag(eps11_S);
Qe31 = real(e31)./imag(e31);
Qe33 = real(e33)./imag(e33);
%Qe15 = real(e15)./imag(e15);
Qc11_E = real(c11_E)./imag(c11_E);
Qc12_E = real(c12_E)./imag(c12_E);
Qc13_E = real(c13_E)./imag(c13_E);
Qc33_E = real(c33_E)./imag(c33_E);
%Qc44_E = real(c44_E)./imag(c44_E);
Qc66_E = real(c66_E)./imag(c66_E);

d31 = num2cell(d31);
[results.d31] = deal(d31{:});
s13_E = num2cell(s13_E);
[results.s13_E ] = deal(s13_E {:});
s33_E = num2cell(s33_E);
[results.s33_E ] = deal(s33_E {:});
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

%Q-verider
Qeps33_S = num2cell(Qeps33_S);
[results.Qeps33_S ] = deal(Qeps33_S {:});
%Qeps11_S = num2cell(Qeps11_S);
%[results.Qeps11_S ] = deal(Qeps11_S {:});
Qe31 = num2cell(Qe31);
[results.Qe31 ] = deal(Qe31 {:});
Qe33 = num2cell(Qe33);
[results.Qe33 ] = deal(Qe33 {:});
%Qe15 = num2cell(Qe15);
%[results.Qe15 ] = deal(Qe15 {:});
Qc11_E = num2cell(Qc11_E);
[results.Qc11_E ] = deal(Qc11_E {:});
Qc12_E = num2cell(Qc12_E);
[results.Qc12_E ] = deal(Qc12_E {:});
Qc13_E = num2cell(Qc13_E);
[results.Qc13_E ] = deal(Qc13_E {:});
Qc33_E = num2cell(Qc33_E);
[results.Qc33_E ] = deal(Qc33_E {:});
Qc66_E = num2cell(Qc66_E);
[results.Qc66_E ] = deal(Qc66_E {:});


