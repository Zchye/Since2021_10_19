function CHParam = ConfigChannel(param, LinkDir, siteIdx, RNTI)
    %LinkDir    -   0 for DL, 1 for UL
    
    %Select proper part of the Table 7.5-6 according to the scenario
    FFTabs = param.FastFadingTabs;
    switch param.Scenario
        case 'UMa'
            Tab = FFTabs.UMa;
        case 'UMi'
            Tab = FFTabs.UMi;
        case 'RMa'
            Tab = FFTabs.RMa;
        otherwise 
            error('Scenarios other than UMa, UMi, RMa are yet to be implemented.');
    end
    
    %Generate channel parameters
    %step 1
    
    %MXC_1 
    if LinkDir == 0 % DL
        TxPos = param.GNBPositions(siteIdx,:);
        RxPos = param.UEPositions{siteIdx}(RNTI,:);
    else % UL
        TxPos = param.UEPositions{siteIdx}(RNTI,:);
        RxPos = param.GNBPositions(siteIdx,:);
    end
    LOSAngles = getLOSAngles(TxPos,RxPos);
    bearing = param.gNBBearing(siteIdx);
    downtilt = param.AntennaDowntilt;
    slant = param.AntennaSlant;
    %step2
    %gNBPos = TxPos;
    %UEPos = RxPos;
    if LinkDir == 0 % DL
        gNBPos = TxPos;
        UEPos = RxPos;
    else % UL
        gNBPos = RxPos;
        UEPos = TxPos;
    end
    %MXC_1 
    
    LOS = assignLOS(param,gNBPos,UEPos);
    %Convert LOS condition from logical to character array
    %If LOS=1, PropCond='LOS'; else PropCond='NLOS'.
    PropCond = strrep(char(LOS*double('LOS ')+(1-LOS)*double('NLOS')),' ','');
    %step 3
    [Pathloss, sigma_SF] = calculatePathloss(param,gNBPos,UEPos,LOS,LinkDir);
    %step 4
    LSPs = genLSPs(param,gNBPos,UEPos,LOS,LinkDir);
    %step 5
    [tau_n,tau_nForStep6] = genTau_n(param,LSPs,LOS);
    %step 6
    [P_n,P_per_ray,RmIdx] = genP_n(param,LOS,LSPs,tau_nForStep6);
    tau_n(RmIdx) = [];  %Remove the delays corresponding to the removed powers
    %step 7
    ADAngles = genADAngles(param, gNBPos, UEPos, LOS, LSPs, P_n, LOSAngles,LinkDir);
    %step 9
    XPR = genXPR(param, LOS);
    %step 11 c_DS
    c_DS = readTab7Dot5Dash6(Tab,'c_DS',PropCond,param.DLCarrierFreq/1e9);
    
    %Construct CHParam
    CHParam.LOSAngles = LOSAngles;
    CHParam.bearing = bearing;
    CHParam.downtilt = downtilt;
    CHParam.slant = slant;
    CHParam.LOS = LOS;
    CHParam.Pathloss = Pathloss;
    CHParam.sigma_SF = sigma_SF;
    CHParam.LSPs = LSPs;
    CHParam.tau_n = tau_n;
    CHParam.P_n = P_n;
    CHParam.P_per_ray = P_per_ray;
    CHParam.ADAngles = ADAngles;
    CHParam.XPR = XPR;
    CHParam.c_DS = c_DS;  

end

%% step 1 c)
function LOSAngles=getLOSAngles(TxPos, RxPos)
    %getLOSAngles calculates the LOS AOD, LOS ZOD, LOS AOA, and LOS ZOA
    %of the directed link between a transmitter and a receiver
    %TxPos - the Cartesian coordinates of the transmitter. Must be a 3-D
    %column vector.
    %RxPos - the Cartesian coordinates of the receiver. Must be a 3-D
    %column vector.
    %LOSAngles - a structure with fields AOD, ZOD, AOA, ZOA. The values are
    %in degrees.
    
    %get the 3-D vector that is from the transmitter to the receiver
    v=RxPos-TxPos;
    
    %transform the Cartesion coordinates of v to spherical coordinates
    [az,el,~]=cart2sph(v(1),v(2),v(3));
    
    %From radian to degree. Note that the zenith angle is complementary to
    %the elevation angle, and that the angle or arrival is supplementary to
    %the angle of departure in the global coordinate system.
    AOD=rad2deg(az);
    ZOD=90-rad2deg(el);%zenith angle is complementary to elevation angle
    AOA=AOD+180;
    ZOA=ZOD+180;
    
    %Confine the angles to be in [0,360)
    LOSAngles.AOD=mod(AOD,360);
    LOSAngles.ZOD=mod(ZOD,360);
    LOSAngles.AOA=mod(AOA,360);
    LOSAngles.ZOA=mod(ZOA,360);
end

%% step 2
function LOS=assignLOS(param,gNBPos,UEPos)
    %assignLOS assigns propagation conditions (LOS/NLOS) for each BS-UT
    %link according to 3GPP TR 38.901 Table 7.4.2-1
    %param - a structure containing the field Scenario
    %gNBPos - the Cartesian coordinates of a gNB, must be a 3-D column
    %vector
    %UEPos - the Cartesion coordinates of a UE, must be a 3-D column vector
    %LOS - 1 for LOS, 0 for NLOS.
    
    %get the parameters in 3GPP TR 38.901 Figure7.4.1-1
    d_2D=norm(gNBPos(1:2)-UEPos(1:2));
    h_UT=UEPos(3);
    
    %select the function to calculate LOS probability
    switch param.Scenario
        case 'RMa'
            h=@RMaProb;
        case 'UMi'
            h=@UMiProb;
        case 'UMa'
            h=@UMaProb;
        otherwise
            error('Scenarios other than RMa, UMi, UMa are yet to be implemented.');
    end
    
    %calculate LOS probability
    Prob=h(d_2D,h_UT);
    
    %generate Bernoulli random variable with the probability Prob
    LOS=binornd(1,Prob);

end

function Prob=RMaProb(d_2D,~)
    %see 3GPP TR 38.901 Table 7.4.2-1
    if d_2D<=10
        Prob=1;
    else
        Prob=exp(-(d_2D-10)/1000);
    end
end

function Prob=UMiProb(d_2D,~)
    %see 3GPP TR 38.901 Table 7.4.2-1
    if d_2D<=18
        Prob=1;
    else
        Prob=18/d_2D+exp(-d_2D/36)*(1-18/d_2D);
    end
end

function Prob=UMaProb(d_2D,h_UT)
    %see 3GPP TR 38.901 Table 7.4.2-1
    if h_UT<=13
        C=0;
    elseif 13<h_UT && h_UT<=23
        C=((h_UT-13)/10)^1.5;
    else
         error('The height of a UE must not exceed 23m.');
    end

    if d_2D<=18
        Prob=1;
    else
        Prob=(18/d_2D+exp(-d_2D/63)*(1-18/d_2D))*(1+C*5/4*(d_2D/100)^3*exp(-d_2D/150));
    end
end

%% step 3
function [Pathloss, sigma_SF]=calculatePathloss(param,gNBPos,UEPos,LOS,LinkDir)
    %Calculate pathloss with formulas in Table 7.4.1-1 for each BS-UT link to be modelled
    %param - a structure containing the field Scenario
    %gNBPos - the Cartesian coordinates of a gNB, must be a 3-D column
    %vector
    %UEPos - the Cartesion coordinates of a UE, must be a 3-D column vector
    %LOS - 1 for LOS, 0 for NLOS.
    %get the parameters in 3GPP TR 38.901 Figure7.4.1-1
    
    %get distance and heights
    d_2D=norm(gNBPos(1:2)-UEPos(1:2));
    d_3D=norm(gNBPos(1:3)-UEPos(1:3));
    h_BS = gNBPos(3);
    h_UT=UEPos(3);
    
    %constants
    if ~LinkDir %For DL
        f_c = param.DLCarrierFreq/1e9; %f_c is carrier frequency in GHz
    else %For UL
        f_c = param.ULCarrierFreq/1e9; %f_c is carrier frequency in GHz
    end
    
    
    %select the function to calculate LOS probability
    switch param.Scenario
        case 'UMi'
            h=@UMiPathloss;
        case 'UMa'
            h=@UMaPathloss;
        case 'RMa'
            h=@RMaPathloss;
        otherwise
            error('Scenarios other than RMa, UMi, UMa are yet to be implemented.');
    end
    
    %calculate Pathloss

    [Pathloss, sigma_SF]=h(d_2D,d_3D,h_BS,h_UT,f_c,LOS);
    
    

end

function [Pathloss, sigma_SF]=UMiPathloss(d_2D,d_3D,h_BS,h_UT,f_c,LOS)
    %see 3GPP TR 38.901 Table 7.4.1-1
    
    %generate coefficients
    c = 300000000;
    h_e = 1;
    h_bs_prime = h_BS - h_e;
    h_ut_prime = h_UT - h_e;
    d_bp_prime = (4 * h_bs_prime * h_ut_prime * f_c * 1000000000) / c;
    
    %calculate LOS pathloss
    if (d_2D >= 10) && (d_2D <= d_bp_prime)
        PL_umi_los = 32.4 + 21*log10(d_3D) + 20*log10(f_c);
    elseif (d_2D >= d_bp_prime) && (d_2D <=5000)
        PL_umi_los = 32.4 + 40*log10(d_3D) + 20*log10(f_c) - 9.5*log10((d_bp_prime)^2+(h_BS-h_UT)^2);
    else
%         error('d_2D out of bound.');
        PL_umi_los = 32.4 + 40*log10(d_3D) + 20*log10(f_c) - 9.5*log10((d_bp_prime)^2+(h_BS-h_UT)^2);
    end 
    
    if LOS
        Pathloss = PL_umi_los;
        sigma_SF = 4;
    else
        %calculate NLOS pathloss
        PL_umi_nlos_prime = 35.3*log10(d_3D)+22.4+21.3*log10(f_c)-0.3*(h_UT-1.5);
        Pathloss = max(PL_umi_los, PL_umi_nlos_prime);
        sigma_SF = 7.82;
    end
end

function [Pathloss, sigma_SF]=UMaPathloss(d_2D,d_3D,h_BS,h_UT,f_c,LOS)
    %see 3GPP TR 38.901 Table 7.4.1-1
    
    %generate coefficients
    c = 300000000;
    h_e = 1;
    h_bs_prime = h_BS - h_e;
    h_ut_prime = h_UT - h_e;
    d_bp_prime = (4 * h_bs_prime * h_ut_prime * f_c * 1000000000) / c;
    
    %calculate LOS pathloss
    if (d_2D >= 10) && (d_2D <= d_bp_prime)
        PL_uma_los = 28.0 + 22*log10(d_3D) + 20*log10(f_c);
    elseif (d_2D >= d_bp_prime) && (d_2D <=5000)
        PL_uma_los = 28.0 + 40*log10(d_3D) + 20*log10(f_c) - 9*log10((d_bp_prime)^2+(h_BS-h_UT)^2);
    else
%         error('d_2D out of bound.');
        PL_uma_los = 28.0 + 40*log10(d_3D) + 20*log10(f_c) - 9*log10((d_bp_prime)^2+(h_BS-h_UT)^2);
    end 
    
    if LOS
        Pathloss = PL_uma_los;
        sigma_SF = 4;
    else
        %calculate NLOS pathloss
        PL_uma_nlos_prime = 13.54 + 39.08*log10(d_3D)+20*log10(f_c)-0.6*(h_UT-1.5);
        Pathloss = max(PL_uma_los, PL_uma_nlos_prime);
        sigma_SF = 6;
    end
end

function [Pathloss, sigma_SF]=RMaPathloss(d_2D,d_3D,h_BS,h_UT,f_c,LOS)
    %see 3GPP TR 38.901 Table 7.4.1-1
    
    %generate coefficients
    c = 300000000;
    h = 5;
    w = 20;
    d_bp = (2 * pi * h_BS * h_UT * f_c * 1000000000) / c;
    
    %calculate LOS pathloss
    if (d_2D >= 10) && (d_2D <= d_bp)
        PL_rma_los = 20*log10(40*pi*d_3D*f_c/3)+min(0.03*(h^1.72),10)*log10(d_3D)-min(0.044*(h^1.72),14.77)+0.002*log10(h)*d_3D;
        sigma = 4;
    elseif (d_2D >= d_bp) && (d_2D <=10000)
        PL_rma_los = 20*log10(40*pi*d_bp*f_c/3)+min(0.03*(h^1.72),10)*log10(d_bp)-min(0.044*(h^1.72),14.77)+0.002*log10(h)*d_bp+40*log10(d_3D/d_bp);
        sigma = 6;
    else
%         error('d_2D out of bound.');
        PL_rma_los = 20*log10(40*pi*d_bp*f_c/3)+min(0.03*(h^1.72),10)*log10(d_bp)-min(0.044*(h^1.72),14.77)+0.002*log10(h)*d_bp+40*log10(d_3D/d_bp);
        sigma = 6;
    end 
    
    if LOS
        Pathloss = PL_rma_los;
        sigma_SF = sigma;
    else
        %calculate NLOS pathloss
        PL_rma_nlos_prime = 161.04-7.1*log10(w)+7.5*log10(h)-(24.34-((h/h_BS)^2))*log10(h_BS)+(43.42-3.1*log10(h_BS))*(log10(d_3D)-3)+20*log10(f_c)-(3.2*((log10(11.75*h_UT))^2)-4.97);
        Pathloss = max(PL_rma_los, PL_rma_nlos_prime);
        sigma_SF = 8;
    end
end

%% step 4
function LSPs=genLSPs(param,gNBPos,UEPos,LOS,LinkDir)
    %Generate large scale parameters. See 3GPP TR 38.901 section 7.5 step 4
    %LinkDir - link direction. 0 for DL, 1 for UL
    %In ns-3, the LSPs are generated by three-gpp-channel-model.cc. They
    %use a MATLAB file to generate the covariance matrix which can be found
    %in 
    %https://github.com/nyuwireless-unipd/ns3-mmwave/blob/master/src/mmwave/model/BeamFormingMatrix/SqrtMatrix.m
    %However, the values used there do not agree with 3GPP TR 38.901 Table
    %7.5-6
    %Here we substitute the matrix with the nearest positive semidefinite
    %matrix.
    d_2D = norm(gNBPos(1:2)-UEPos(1:2));
    h_UT = UEPos(3);
    h_BS = gNBPos(3);
    if ~LinkDir %For DL
        fc = param.DLCarrierFreq/1e9; %fc is carrier frequency in GHz
    else %For UL
        fc = param.ULCarrierFreq/1e9; %fc is carrier frequency in GHz
    end
    
    %Select tables/table parts according to the specified scenario
    switch param.Scenario
        case 'UMi'
            T=param.FastFadingTabs.UMi;
            W=param.FastFadingTabs.ZSDZODOffsetUMi;
        case 'UMa'
            T=param.FastFadingTabs.UMa;
            W=param.FastFadingTabs.ZSDZODOffsetUMa;
        case 'RMa'
            T=param.FastFadingTabs.RMa;
            W=param.FastFadingTabs.ZSDZODOffsetRMa;
        otherwise
            error('Scenarios other than UMi, UMa, RMa are yet to be implemented.');
    end
    
    %Convert LOS condition from logical to character array
    %If LOS=1, PropCond='LOS'; else PropCond='NLOS'.
    PropCond=strrep(char(LOS*double('LOS ')+(1-LOS)*double('NLOS')),' ','');
    
    %Construct covariance table (CovT)
    if LOS
        parList={'SF','K','DS','ASD','ASA','ZSD','ZSA'};%The order of parameters follows 3GPP TR 38.901 section 7.5 step 4
    else %NLOS, no Riccian K factor
        parList={'SF','DS','ASD','ASA','ZSD','ZSA'};
    end
    Npar=length(parList);%the number of parameters
    varTypes=cell(1,Npar);
    varTypes(:)={'double'};
    CovT=table('Size',[Npar,Npar],'VariableTypes',varTypes);%Initialize CovT
    CovT.Properties.RowNames=parList;
    CovT.Properties.VariableNames=parList;
    
    %Select data from Table 7.5-6 to populate CovT
    mu=zeros(Npar,1);%Preallocate means
    for ii=1:Npar
        %fill variances and means
        rn=parList{ii};%row name
        U=T;%Table to read, all LSPs are in table T except ZSD
        if isequal(rn,'ZSD')
            U=W;%ZSD is in table W
            pn='sigma_lg_ZSD';
            pnmu='mu_lg_ZSD';
            switch param.Scenario
                case 'UMa'
                    %Read sigma and comput variance
                    CovT{rn,rn} = (U{pn,PropCond}{1})^2;
                    %Read mu
                    mu(ii) = U{pnmu,PropCond}{1}(d_2D,h_UT);
                case 'UMi'
                    %Read sigma and comput variance
                    CovT{rn,rn} = (U{pn,PropCond}{1})^2;
                    %Read mu
                    mu(ii) = U{pnmu,PropCond}{1}(d_2D,h_UT,h_BS);
                case 'RMa'
                    %Read sigma and comput variance
                    CovT{rn,rn} = (U{pn,PropCond}{1})^2;
                    %Read mu
                    mu(ii) = U{pnmu,PropCond}{1}(d_2D,h_UT);
                otherwise 
                    error('Scenarios other than UMi, UMa, RMa are yet to be implemented.');
            end
        elseif isequal(rn,'SF') || isequal(rn,'K')
            pn=['sigma_',rn];%parameter name
            pnmu=['mu_',rn];
            %Read sigma and comput variance
            CovT{rn,rn}=(readTab7Dot5Dash6(U,pn,PropCond,fc))^2;
            %Read mu
            mu(ii)=readTab7Dot5Dash6(U,pnmu,PropCond,fc);
        else
            pn=['sigma_lg_',rn];
            pnmu=['mu_lg_',rn];
            %Read sigma and comput variance
            CovT{rn,rn}=(readTab7Dot5Dash6(U,pn,PropCond,fc))^2;
            %Read mu
            mu(ii)=readTab7Dot5Dash6(U,pnmu,PropCond,fc);
        end
        %fill covariances
        for jj=1:Npar
            cn=parList{jj};%column name
            pn=[rn,'vs',cn];
            %Copy cross-correlation to covariance
            CovT{rn,cn}=readTab7Dot5Dash6(T,pn,PropCond,fc);
        end
    end
   
    %Read sigma_SF and PL
    [Pathloss, sigma_SF]=calculatePathloss(param,gNBPos,UEPos,LOS);
    CovT{'SF','SF'}=sigma_SF^2;
    mu(1)=Pathloss;
    
    %Create covariance matrix (CovM) from CovT. Note that now CovM is lower
    %triangular
    CovM=table2array(CovT);
    
    %Make CovM symmetric by flipping the lower part to the upper part
    CovM=CovM+triu(CovM',1);
     
    %Cholesky decomposition. L*L'=CovM
    L=cholDcmp(CovM,'lower');
    
    %Generate Npar iid N(0,1) Gaussian random variables
    x=randn(Npar,1);
    
    %Generate Gaussian random vector with covariance matrix CovM
    y=L*x;
    s=y+mu;
    
    %Generate LSPs with cross-correlations according to 3GPP TR 38.901
    %Table 7.5-6
    for ii=1:Npar
        switch parList{ii}
            case 'SF'
                LSPs.SF=s(ii);  %in dB
            case 'K'
                LSPs.K=s(ii);   %in dB
            case 'DS'
                LSPs.DS=10^s(ii);
            case 'ASD'
                LSPs.ASD=min(10^s(ii),104);
            case 'ASA'
                LSPs.ASA=min(10^s(ii),104);
            case 'ZSD'
                LSPs.ZSD=min(10^s(ii),52);
            case 'ZSA'
                LSPs.ZSA=min(10^s(ii),52);
        end
    end
end

function y=readTab7Dot5Dash6(T,row,col,fc)
    %Read a value from table 'T', specified by row name 'row' and column name
    %'col'
    %The data type in the tables could be double or function handle
    %If the data is a function handle, the function should be evaluated
    %at fc
    
    %Check if the row exist. If not, read 0.
    if isempty(find(strcmp(T.Properties.RowNames,row),1))
        y=0;
        return
    end
    
    %If the data is a function handle, the function should be evaluated
    %at fc
    temp=T{row,col}{1};
    if isa(temp,'function_handle')
        y=temp(fc);
    elseif isa(temp,'double')
        y=temp;
    else
        y=0;
    end
end

%% step 5
function [tau_n,tau_nForStep6] = genTau_n(param,LSPs,LOS)
    %Generate cluster delays tau_n. See 3GPP TR 38.901 section 7.5 step 5.
    
    %Convert LOS condition from logical to character array
    %If LOS=1, PropCond='LOS'; else PropCond='NLOS'.
    PropCond = strrep(char(LOS*double('LOS ')+(1-LOS)*double('NLOS')),' ','');
    
    FastFadingTabs = param.FastFadingTabs;
    
    %Select table according to the scenario
    switch param.Scenario
        case 'RMa'
            Tab = FastFadingTabs.RMa;
        case 'UMa'
            Tab = FastFadingTabs.UMa;
        case 'UMi'
            Tab = FastFadingTabs.UMi;
        otherwise
            error('Scenarios other than RMa, UMa, UMi are yet to be implemented.');
    end
    
    %Read parameters from the table selected
    N_cluster = Tab{'N_cluster',PropCond}{1};
    r_tau = Tab{'r_tau',PropCond}{1};
    
    %Delay spread
    DS = LSPs.DS;
    
    %Generate tau_n
    %MXC_1 quick fix 
    %MXC_2 made smaller
    %X = rand(N_cluster,1);
    X = 0.4 * rand(N_cluster,1) + 0.6;
    %MXC_1
    tau_prime = -r_tau*DS*log(X);
    tau_n = sort(tau_prime-min(tau_prime));
    
    %Generate tau_n to be used for generating cluster powers(step 6)
    tau_nForStep6=tau_n;
    
    %In case of LOS, scale the delays
    if LOS
        K = LSPs.K;
        C_tau = 0.7705-0.0433*K+0.0002*K^2+0.000017*K^3;
        tau_n = tau_n/C_tau;
    end
end

%% step 6
function [P_n,P_per_ray,RmIdx] = genP_n(param,LOS,LSPs,tau_n)
    %Generate cluster powers P_n. See 3GPP TR 38.901 section 7.5 step 6.

    %Convert LOS condition from logical to character array
    %If LOS=1, PropCond='LOS'; else PropCond='NLOS'.
    PropCond = strrep(char(LOS*double('LOS ')+(1-LOS)*double('NLOS')),' ','');

    %Select table according to the scenario
    switch param.Scenario
        case 'UMi'
            Table=param.FastFadingTabs.UMi;
        case 'UMa'
            Table=param.FastFadingTabs.UMa;
        case 'RMa'
            Table=param.FastFadingTabs.RMa;
        otherwise
            error('Scenarios other than UMi, UMa, RMa are yet to be implemented.');
    end

    %Read the per cluster shadowing std [dB]
    zeta_pcshadow = Table{'zeta_pcshadow',PropCond}{1};
    %Read number of clusters
    N_cluster = Table{'N_cluster',PropCond}{1};
    %Read number of rays per cluster
    M_rays_per_cluster = Table{'M_raypc',PropCond}{1};
    %Read delay scaling parameter
    r_tau = Table{'r_tau',PropCond}{1};

    %Delay spread
    DS = LSPs.DS;

    %Ricean K-factor converted to linear scale
    %TR38.901 uses 10log10, MATLAB db2mag uses 20log10
    if LOS
    K_R = 10 ^ (LSPs.K / 10); 
    end

    %per cluster shadowing term in [dB]
    Z_n = normrnd(0,zeta_pcshadow,[N_cluster,1]);

    %cluster powers
    P_prime= exp(-tau_n * ((r_tau-1) / (r_tau * DS))).*(10.^(-Z_n/10));
    %i think this is right, double check for ,* and .^

    %normalize the cluster powers so that the sum power of all cluster powers is equal to one
    P_n = P_prime / sum(P_prime);

    if LOS
        P_1_LOS = K_R / (K_R + 1);
        P_n = P_n * (1 / (K_R + 1));
        P_n(1) = P_n(1) + P_1_LOS;
    end

    %remove clusters with less than -25 dB power compared to maximum cluster power
    %
    %MXC_1 quick fix
    %this causes error on line 1871 in nrCDLChannel, remove for now, fix later
    %
    %
    %
    
    %maximum cluster power
    P_max = max(P_n);

    %threshold power 
    P_threshold = P_max * (10 ^(-2.5));

    %find the indices of and remove clusters with less than -25 dB power compared to maximum cluster power
    RmIdx = find(P_n < P_threshold);
    P_n(P_n < P_threshold) = [];
    
    %power per ray
    P_per_ray=P_n / M_rays_per_cluster;
    
end

%% step 7
function ADAngles = genADAngles(param, gNBPos,UEPos, LOS, LSPs, P, LOSAngles, LinkDir)
    %Generate arrival angles and departure angles for both azimuth and
    %elevation
    
    if ~LinkDir %For DL
        fc = param.DLCarrierFreq/1e9; %Carrier frequency in GHz
    else
        fc = param.ULCarrierFreq/1e9;
    end
    InputStruct.fc = fc;
    InputStruct.gNBPos = gNBPos;
    InputStruct.UEPos = UEPos;
    
    %Construct InputStruct    
    InputStruct.ASA = LSPs.ASA;
    InputStruct.ASD = LSPs.ASD;
    InputStruct.ZSA = LSPs.ZSA;
    InputStruct.ZSD = LSPs.ZSD;
    if LOS
    InputStruct.K = LSPs.K;
    end
    InputStruct.P = P;
    InputStruct.LOS = LOS;
    InputStruct.param = param;
    InputStruct.AOD = LOSAngles.AOD;
    InputStruct.ZOD = LOSAngles.ZOD;
    InputStruct.AOA = LOSAngles.AOA;
    InputStruct.ZOA = LOSAngles.ZOA;
    
    %Convert LOS condition from logical to character array
    %If LOS=1, PropCond='LOS'; else PropCond='NLOS'.
    InputStruct.PropCond = strrep(char(LOS*double('LOS ')+(1-LOS)*double('NLOS')),' ','');

    %Select proper part of Table 7.5-6 according to the scenario
    FastFadingTabs = param.FastFadingTabs;
    switch param.Scenario
        case 'RMa'
            Tab = FastFadingTabs.RMa;
        case 'UMa'
            Tab = FastFadingTabs.UMa;
        case 'UMi'
            Tab = FastFadingTabs.UMi;
        otherwise
            error('Scenarios other than RMa, UMa, UMi are yet to be implemented.');
    end
    InputStruct.Tab = Tab;
    
    %Read parameters from the table selected
    InputStruct.N_cluster = Tab{'N_cluster',InputStruct.PropCond}{1};
    
    %Generate AOAs
    InputStruct.Spread = 'ASA';
    OutputStruct = genAnglesCodeBlock(InputStruct);
    ADAngles.phi_n_AOA = OutputStruct.phi_n_AOA;
    ADAngles.phi_n_m_AOA = OutputStruct.phi_n_m_AOA;
    
    %Generate AODs
    InputStruct.Spread = 'ASD';
    OutputStruct = genAnglesCodeBlock(InputStruct);
    ADAngles.phi_n_AOD = OutputStruct.phi_n_AOD;
    ADAngles.phi_n_m_AOD = OutputStruct.phi_n_m_AOD;
    
    %Generate ZOAs
    InputStruct.Spread = 'ZSA';
    OutputStruct = genAnglesCodeBlock(InputStruct);
    ADAngles.theta_n_ZOA = OutputStruct.theta_n_ZOA;
    ADAngles.theta_n_m_ZOA = OutputStruct.theta_n_m_ZOA;
    
    %Generate ZODs
    InputStruct.Spread = 'ZSD';
    OutputStruct = genAnglesCodeBlock(InputStruct);
    ADAngles.theta_n_ZOD = OutputStruct.theta_n_ZOD;
    ADAngles.theta_n_m_ZOD = OutputStruct.theta_n_m_ZOD;
       
end

function OutputStruct = genAnglesCodeBlock(InputStruct)
    %Compute the cluster angle and the angles of the rays within the
    %cluster for one of the direction -- AOA, AOD, ZOA, ZOD
    
    if InputStruct.LOS
    K = InputStruct.K;
    end
    P = InputStruct.P;
    N = length(P);
    X = 2*(randi(2,1,N)-1.5); %X is a row vector of length N with elements being 1 and -1
    N_cluster = InputStruct.N_cluster;
    alpha_m = [0.0447;
                      0.1413;
                      0.2492;
                      0.3715;
                      0.5129;
                      0.6797;
                      0.8844;
                      1.1481;
                      1.5195;
                      2.1551]; %Table 7.5-3. The plus-minus signs are assigned in the following few lines
    temp = zeros(20,1);
    for m = 1:20
        temp(m) = alpha_m(ceil(m/2))*(-1)^(m+1);    %Assign plus and minus signs
    end
    alpha_m = temp;
    
    if isequal(InputStruct.Spread(1),'A') %Azimuth angles
        %Table 7.5-2
        num_clusters = [4,5,8,10,11,12,14,15,16,19,20,25];
        C = [0.779, 0.86, 1.018, 1.09, 1.123, 1.146, 1.190, 1.211, 1.226, 1.273, 1.289, 1.358];
        %Store angular spreads and LOS angles in sprd and LOSA respectively
        switch InputStruct.Spread
            case 'ASA'
                sprd = InputStruct.ASA;
                LOSA = InputStruct.AOA; %phi_n_AOX
            case 'ASD'
                sprd = InputStruct.ASD;
                LOSA = InputStruct.AOD;
        end
        
        %Generate C_phi. Equation 7.5-10
        C_phi_NLOS = C(find(num_clusters == N_cluster,1));
        if isempty(C_phi_NLOS)  %Validate the number of clusters. If empty, N_cluster is invalid
            error('Invalid number of clusters. Please refer to 3GPP TR 38.901 Table 7.5-2');
        end
        if InputStruct.LOS
            C_phi = C_phi_NLOS*(1.1035-0.028*K-0.002*K^2+0.0001*K^3);
        else
            C_phi = C_phi_NLOS;
        end
        
        %Equantion 7.5-9
        phi_n_prime = 2*sprd/1.4*sqrt(-log(P/max(P)))/C_phi;
        
        %Euqation 7.5-11 and 7.5-12
        Y = normrnd(0,sprd/7,1,N);   %Y is a row vector of length N with elements being R.V.s ~ N(0,(sprd/7)^2)
        if find(size(phi_n_prime) ~= 1) == 1 %If phi_n_prime is a column vector, make X and Y column vectors as well.
            X = X';
            Y = Y';
        end
        if ~InputStruct.LOS %For NLOS
            phi_n=X.*phi_n_prime+Y+LOSA;    %Equation 7.5-11
        else %For LOS
            phi_n=(X.*phi_n_prime+Y)-(X(1)*phi_n_prime(1)+Y(1)-LOSA);   %Equation 7.5-12
        end
        
        %Equation 7.5-13
        phi_n_m = zeros(N,20);
        switch InputStruct.Spread
            case 'ASA'
                c = InputStruct.Tab{'c_ASA',InputStruct.PropCond}{1}; %Read c_ASA from the table
            case 'ASD'
                c = InputStruct.Tab{'c_ASD',InputStruct.PropCond}{1}; %Read c_ASD from the table
        end
        for n = 1:N
            for m = 1:m
                phi_n_m(n,m) = phi_n(n) + c*alpha_m(m); %Equation 7.5-13
            end
        end
        
        %Construct output
        switch InputStruct.Spread
            case 'ASA'
                OutputStruct.phi_n_AOA = phi_n;
                OutputStruct.phi_n_m_AOA = phi_n_m;
            case 'ASD'
                OutputStruct.phi_n_AOD = phi_n;
                OutputStruct.phi_n_m_AOD = phi_n_m;
        end
        
    else   %Zenith Angles
        num_clusters = [8, 10, 11, 12, 15, 19, 20, 25];
        C = [0.889, 0.957, 1.031, 1.104, 1.1088, 1.184, 1.178, 1.282];
        
        C_theta_NLOS = C(find(num_clusters == N_cluster,1));
        if isempty(C_theta_NLOS)  %Validate the number of clusters. If empty, N_cluster is invalid
            error('Invalid number of clusters. Please refer to 3GPP TR 38.901 Table 7.5-4');
        end
        if InputStruct.LOS %Equation 7.5-15
            C_theta = C_theta_NLOS*(1.3086+0.0339*K-0.0077*K^2+0.0002*K^3);
        else
            C_theta = C_theta_NLOS;
        end
        
        %Function handle for equation 7.5-14
        h_theta_n_prime = @(sprd) -sprd*log(P/max(P))/C_theta;
        
        if isequal(InputStruct.Spread,'ZSA')%For ZSA
            ZSA = InputStruct.ZSA;
            ZOA = InputStruct.ZOA;
            
            %Equation 7.5-14
            theta_n_prime = h_theta_n_prime(ZSA);
            
            %Generate zenith angles
            Y = normrnd(0,ZSA/7,1,N);   %Y is a row vector of length N with elements being R.V.s ~ N(0,(ZSA/7)^2)
            if find(size(theta_n_prime) ~= 1) == 1 %If theta_n_prime is a column vector, make X and Y column vectors as well.
                X = X';
                Y = Y';
            end
            if ~InputStruct.LOS %For NLOS
                theta_n = X.*theta_n_prime+Y+ZOA;   %Equation 7.5-16
            else %for LOS
                theta_n = (X.*theta_n_prime+Y)-(X(1)*theta_n_prime(1)+Y(1)-ZOA);    %Equation 7.5-17
            end
            
            %Read c_ZSA from the table
            c_ZSA = InputStruct.Tab{'c_ZSA',InputStruct.PropCond}{1};
            
            %Equation 7.5-18
            theta_n_m=zeros(N,20);
            for n = 1:N
                for m = 1:20
                    theta_n_m(n,m) = theta_n(n) + c_ZSA*alpha_m(m);
                end
            end
            theta_n_m=mod(theta_n_m,360);   %Wrap theta_n_m into [0,360]
            if theta_n_m>=180
                theta_n_m = 360 - theta_n_m;
            end
            
            %Construct output
            OutputStruct.theta_n_ZOA = theta_n;
            OutputStruct.theta_n_m_ZOA = theta_n_m;
            
        else %For ZSD
            %Stote LOS angles
            ZSD = InputStruct.ZSD;
            ZOD = InputStruct.ZOD;
            
            %Equation 7.5-14
            theta_n_prime = h_theta_n_prime(ZSD);
            
            
            Y = normrnd(0,ZSD/7,1,N);   %Y is a row vector of length N with elements being R.V.s ~ N(0,(ZSD/7)^2)
            if find(size(theta_n_prime) ~= 1) == 1 %If theta_n_prime is a column vector, make X and Y column vectors as well.
                X = X';
                Y = Y';
            end
            
            gNBPos = InputStruct.gNBPos;
            UEPos = InputStruct.UEPos;
            LOS = InputStruct.LOS;
            param = InputStruct.param;
            fc = InputStruct.fc;
            
            if ~InputStruct.LOS %For NLOS
                rowname = 'mu_offset_ZOD';
                mu_offset_ZOD = readZODTab(param,rowname,gNBPos,UEPos,LOS,fc);
                theta_n = X.*theta_n_prime + Y + ZOD + mu_offset_ZOD;   %Equation 7.5-19
            else %For LOS
                theta_n = (X.*theta_n_prime+Y)-(X(1)*theta_n_prime(1)+Y(1)-ZOD);    %Equation 7.5-17
            end
            
            %Equation 7.5-20
            rowname = 'mu_lg_ZSD';
            mu_lg_ZSD = readZODTab(param,rowname,gNBPos,UEPos,LOS,fc);
            theta_n_m = zeros(N,20);
            for n = 1:N
                for m = 1:20
                    theta_n_m(n,m) = theta_n(n) + (3/8)*(10^mu_lg_ZSD)*alpha_m(m);
                end
            end
            
            %Construct output
            OutputStruct.theta_n_ZOD = theta_n;
            OutputStruct.theta_n_m_ZOD = theta_n_m;
        end
    end
                           
end

function y = readZODTab(param,rowname,gNBPos,UEPos,LOS,fc)
    %Read data from 3GPP TR 38.901 Table 7.5-7, 7.5-8, 7.5-9
    %fc is in GHz

    d_2D = norm(gNBPos(1:2)-UEPos(1:2));
    h_UT = UEPos(3);
    h_BS = gNBPos(3);
    %Convert LOS condition from logical to character array
    %If LOS=1, PropCond='LOS'; else PropCond='NLOS'.
    PropCond = strrep(char(LOS*double('LOS ')+(1-LOS)*double('NLOS')),' ','');
    
    switch param.Scenario
        case 'UMa'
            Tab = param.FastFadingTabs.ZSDZODOffsetUMa;
            switch rowname
                case 'mu_lg_ZSD'
                    y = Tab{rowname,PropCond}{1}(d_2D,h_UT);
                case 'sigma_lg_ZSD'
                    y = Tab{rowname,PropCond}{1};
                case 'mu_offset_ZOD'
                    if LOS
                        y = 0;
                    else
                        y = Tab{rowname,PropCond}{1}(fc,d_2D,h_UT);
                    end
            end
        case 'UMi'
            Tab = param.FastFadingTabs.ZSDZODOffsetUMi;
            switch rowname
                case 'mu_lg_ZSD'
                    y = Tab{rowname,PropCond}{1}(d_2D,h_UT,h_BS);
                case 'sigma_lg_ZSD'
                    y = Tab{rowname,PropCond}{1};
                case 'mu_offset_ZOD'
                    if LOS
                        y = 0;
                    else
                        y = Tab{rowname,PropCond}{1}(d_2D);
                    end
            end
        case 'RMa'
            Tab = param.FastFadingTabs.ZSDZODOffsetRMa;
            switch rowname
                case 'mu_lg_ZSD'
                    y = Tab{rowname,PropCond}{1}(d_2D,h_UT);
                case 'sigma_lg_ZSD'
                    y = Tab{rowname,PropCond}{1};
                case 'mu_offset_ZOD'
                    if LOS
                        y = 0;
                    else
                        y = Tab{rowname,PropCond}{1}(d_2D);
                    end
            end
    end
    
end

%% step 8
%Already implemented in nrCDLChannel

%% step 9
function XPR = genXPR(param, LOS)
    %Generate the cross polarization power ratios (XPR)
    
    Tabs = param.FastFadingTabs;
    
    %Select proper table part from Table 7.5-6 according to the scenario
    switch param.Scenario
        case 'UMa'
            Tab = Tabs.UMa;
        case 'UMi'
            Tab = Tabs.UMi;
        case 'RMa'
            Tab = Tabs.RMa;
        otherwise 
            error('Scenarios other than UMa, UMi, RMa are yet to be implemented.');
    end
    
    %Convert LOS condition from logical to character array
    %If LOS=1, PropCond='LOS'; else PropCond='NLOS'.
    PropCond = strrep(char(LOS*double('LOS ')+(1-LOS)*double('NLOS')),' ','');
    
    %Assign the mean (in dB) to XPR.
    %This is NOT compliant with section 7.5 step 9 since it requires XPR
    %drawn independently for each ray and each cluster but nrCDLChannel
    %requires a numeric scalar for the porperty XPR.
    XPR = Tab{'mu_XPR',PropCond}{1};
end

%% Additional functions
function B = nearestPosSemDef(A)
    %Find the nearest positive semidefinite matrix B from a given matrix A
    %in the Frobenius norm
    %Reference: Nicholas J. Higham, Copmuting a Nearest Symmetric Positive
    %Semidefinite Matrix, Linear Algebra and its Applications, Volume 103, 1988, Pages 103-118, ISSN 0024-3795, https://doi.org/10.1016/0024-3795(88)90223-6.
    
    A = (A+A')/2; %Make A symmetric
    [V,D] = eig(A); %Eigenvalue decomposition
    n = size(D,1);
    e = 1e-9;
    %Make the eigenvalues of A nonnegative
    for ii = 1:n
        if D(ii,ii) <0
            D(ii,ii) = e;
        end
    end
    %Create B from the reconstructed D
    B = V*D/V;
    B = (B+B')/2; %To compensate for rounding errors
end

function C = cholDcmp(varargin)
    %Implement Cholesky decomposition
    A = varargin{1};
    while ~all(eig(A)>=0)
        A = nearestPosSemDef(A);
    end
    if nargin == 1
        C = chol(A);
    else
        C = chol(A,varargin{2});
    end
            
end
