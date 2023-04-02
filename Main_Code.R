# Loading R packages
library(mrgsolve)
library(dplyr)
library(zoo)  # the package was loaded for the function na.locf 
library(ggplot2)
library(patchwork)
library(FME)
library(Metrics)
library(tidyr)

# Define a mrgsolve-based pbpk
NanoPBPK.code <- '

$PARAM @annotated

// Cardiac Output and Blood flow
QCC            : 16.5   :L/h/kg^0.75,                    Cardio output,                  (Brown, 1997)
QLC            : 0.02   :unitless,                       Fraction blood flow to liver,   (Brown, 1997, Table 23)
QLuC           : 1      :unitless,                       Fraction blood flow to lung,    (Brown, 1997, Table 23)
QKC            : 0.091  :unitless,                       Fraction blood flow to kidney,  (Brown, 1997, Table 23)
QBrC           : 0.033  :unitless,                       Fraction blood flow to brain,   (Brown, 1997, Table 23)
QSC            : 0.011  :unitless,                       Fraction blood flow to spleen,  (Lin, 2008; Davies and Morries, 1993)
QMC            : 0.159  :unitless,                       Fraction blood flow to muschle, (Brown, 1997, Table 23)
QTC            : 0.033  :unitless,                       Fraction blood flow to tumor,   

// Tissue Volume
BW             : 0.02   :kg,                             Body weight                    
VLC            : 0.055  :unitless,                       Fraction liver tissue,          (Brown, 1997, Table 21)
VLuC           : 0.007  :unitless,                       Fraction lung tissue,           (Brown, 1997, Table 21)
VKC            : 0.017  :unitless,                       Fraction kidney tissue,         (Brown, 1997, Table 21)
VBrC           : 0.017  :unitless,                       Fraction brain tissue,          (Brown, 1997, Table 21)
VSC            : 0.005  :unitless,                       Fraction spleen tissue,         (Lin, 2008; Davies and Morries, 1993)
VBldC          : 0.06   :unitless,                       Fraction blood,                 (Chen, 2015)
VPlasC         : 0.0355 :unitless,                       Fraction plasma,                (Davies and Morris, 1993
VMC            : 0.384  :unitless,                       Fraction muscle tissue,         (Brown, 1997, Table 21)
VTC            : 0.04   :unitless,                       Fraction tumor,                 (Sykes et al., 2014; Wilhelm et al., 2016)


//Blood volume fraction in organs and tissues
BVL            : 0.31   :unitless,                       Liver,                          (Brown, 1997, Table 30)
BVBr           : 0.03   :unitless,                       Brain,                          (Brown, 1997, Table 30)
BVK            : 0.24   :unitless,                       Kidney,                         (Brown, 1997, Table 30)
BVS            : 0.17   :unitless,                       Spleen,                         (Brown, 1997, Table 30)
BVLu           : 0.5    :unitless,                       lungs,                          (Brown, 1997, Table 30)
BVM            : 0.04   :unitless,                       muscle,                         (Brown, 1997, Table 30)
BVR            : 0.04   :unitless,                       Rest of body (assumed to equal to muscle), (Brown, 1997, Table 30)
BVT            : 0.03   :unitless,                       tumor,                          fitted

//Partition coefficients(PC, tissue:plasma)
PL             : 0.08   :unitless,                       liver,                          (Lin, 2016)
PK             : 0.15   :unitless,                       kidney,                         (Lin, 2016)
PBr            : 0.15   :unitless,                       brain,                          (Lin, 2016)
PS             : 0.15   :unitless,                       spleen,                         (Lin, 2016)
PLu            : 0.15   :unitless,                       lungs,                          (Lin, 2016)
PH             : 0.15   :unitless,                       heart,                          (Lin, 2016)
PM             : 0.15   :unitless,                       muscle,                         (Lin, 2016)
PR             : 0.15   :unitless,                       rest of body,                   (Lin, 2016)
PT             : 0.15   :unitless,                       tumor,                          fitted

//Membrane-limited permeability coefficient constants
PALC           : 0.001  :unitless,                       liver,                          (Lin, 2016)
PABrC          : 0.000001 :unitless,                     brain,                          (Lin, 2016)
PAKC           : 0.01   :unitless,                       kidney,                         (Lin, 2016)
PASC           : 0.15   :unitless,                       spleen,                         (Lin, 2016)
PALuC          : 0.001  :unitless,                       lung,                           (Lin, 2016)
PAMC           : 0.00005:unitless,                       muscle,                         (Lin, 2016)
PARC           : 0.00005:unitless,                       rest of body,                   (Lin, 2016)
PATC           : 0.001  :unitless,                       tumor,                          fitted

//Endocytic parameters; RES represent endocytic/phagocytic cells
KLRESrelease   : 0.0015   :1/h,                liver,                          Release rate constant of phyagocytic cells 
KLRESmax       : 0.3      :1/h,                liver,                          Maximmum uptake rate constant of phyagocytic cells
KLRES50        : 48       :h,                  liver,                          Time reaching half maximum uptake rate
KLRESn         : 5        :unitless,           liver,                          Hill coefficient
ALREScap       : 100      :ug/g tissue,        liver,                          uptake capacity per tissue weight

KSRESrelease   : 0.001    :1/h,                spleen,                         Release rate constant of phyagocytic cells 
KSRESmax       : 5        :1/h,                spleen,                         Maximmum uptake rate constant of phyagocytic cells
KSRES50        : 36       :h,                  spleen,                         Time reaching half maximum uptake rate
KSRESn         : 5        :unitless,           spleen,                         Hill coefficient
ASREScap       : 200      :ug/g tissue,        spleen,                         uptake capacity per tissue weight

KKRESrelease   : 0.001    :1/h,                kidney,                         Release rate constant of phyagocytic cells 
KKRESmax       : 0.12     :1/h,                kidney,                         Maximmum uptake rate constant of phyagocytic cells
KKRES50        : 48       :h,                  kidney,                         Time reaching half maximum uptake rate
KKRESn         : 5        :unitless,           kidney,                         Hill coefficient
AKREScap       : 15       :ug/g tissue,        kidney,                         uptake capacity per tissue weight

KLuRESrelease  : 0.003    :1/h,               lung,                           Release rate constant of phyagocytic cells 
KLuRESmax      : 0.085    :1/h,               lung,                           Maximmum uptake rate constant of phyagocytic cells
KLuRES50       : 48       :h,                 lung,                           Time reaching half maximum uptake rate
KLuRESn        : 5        :unitless,          lung,                           Hill coefficient
ALuREScap      : 100      :ug/g tissue,       lung,                           uptake capacity per tissue weight

KMRESrelease   : 0.005    :1/h,                muscle,                         Release rate constant of phyagocytic cells 
KMRESmax       : 0.4      :1/h,                muscle,                         Maximmum uptake rate constant of phyagocytic cells
KMRES50        : 48       :h,                  muscle,                         Time reaching half maximum uptake rate
KMRESn         : 5        :unitless,           muscle,                         Hill coefficient
AMREScap       : 15       :ug/g tissue,        muscle,                         uptake capacity per tissue weight

KRRESrelease   : 0.005    :1/h,                rest of body,                   Release rate constant of phyagocytic cells 
KRRESmax       : 0.4      :1/h,                rest of body,                   Maximmum uptake rate constant of phyagocytic cells
KRRES50        : 48       :h,                  rest of body,                   Time reaching half maximum uptake rate
KRRESn         : 5        :unitless,           rest of body,                   Hill coefficient
ARREScap       : 15       :ug/g tissue,        rest of body,                   uptake capacity per tissue weight

KTRESrelease   : 0.005    :1/h,                tumor, Release rate constant of phyagocytic cells 
KTRESmax       : 0.4      :1/h,                tumor, Maximmum uptake rate constant of phyagocytic cells
KTRES50        : 24       :h,                  tumor, Time reaching half maximum uptake rate
KTRESn         : 5        :unitless,           tumor, Hill coefficient
ATREScap       : 1        :ug/g tissue,         tumor, uptake capacity per tissue weight

//Excretion parameters
KbileC         :  0.00003  :L/hr/kg^0.75,       Bile clearance
KurineC        :  0.000003 :L/hr/kg^0.75,       Urine clearance

$MAIN
double QRC     = 1 - (QLC + QKC  + QSC + QTC + QBrC + QMC);                          //Fraction of blood flow to rest of body
double VRC     = 1 - (VLC + VLuC + VKC  + VSC + VBrC + VMC + VTC + VPlasC);                          //Tissue volume of rest of body

double QC      = QCC * pow(BW, 0.75);                                                    //L/h, Cardiac output (adjusted for plasma)
double QL      = QC * QLC;                                                           //L/h, Blood flow to liver
double QBr     = QC * QBrC;                                                         //L/h, Blood flow to brain
double QK      = QC * QKC;                                                           //L/h, Blood flow to kidney
double QS      = QC * QSC;                                                          //L/h, Blood flow to spleen
double QM      = QC * QMC;                                                           //L/h, Blood flow to muscle
double QT      = QC * QTC;                                                           //L/h, Blood flow to tumor
double QR      = QC * QRC;                                                           //L/h, Blood flow to the rest of body

//Tissue volumes
double VL      = BW * VLC;                                                           //L, Liver volume
double VBr     = BW * VBrC;                                                          //L, Brain volume
double VK      = BW * VKC;                                                           //L, Kidney volume
double VM      = BW * VMC;
double VS      = BW * VSC;                                                          //L, spleen volume
double VLu     = BW * VLuC;                                                          //L, lung volume
double VR      = BW * VRC;                                                           //L, Rest of body
double VT      = BW * VTC;                                                           //L, Tumor
double VBlood  = BW * VBldC;                                                      //L, Blood
double VPlasma = BW * VPlasC;                                                    //L, Plasma

double VLb     = VL * BVL; 							                                            //Weight/volume of capillary blood in liver compartment
double VLt     = VL - VLb; 							                                            //Weight/volume of tissue in liver compartment
double VBrb    = VBr * BVBr; 						                                            //Weight/volume of capillary blood in brain compartment
double VBrt    = VBr - VBrb; 						                                            //Weight/volume of tissue in brain compartment
double VKb     = VK * BVK; 							                                            //Weight/volume of capillary blood in kidney compartment
double VKt     = VK - VKb; 							                                            //Weight/volume of tissue in kidney compartment
double VSb     = VS * BVS; 							                                          //Weight/volume of capillary blood in spleen compartment
double VSt     = VS - VSb; 							                                          //Weight/volume of tissue in spleen compartment
double VLub    = VLu * BVLu; 						                                            //Weight/volume of capillary blood in lung compartment
double VLut    = VLu - VLub; 						                                            //Weight/volume of tissue in lung compartment
double VMb     = VM * BVM; 						                                              //Weight/volume of capillary blood in muscle compartment
double VMt     = VM - VMb; 							                                            //Weight/volume of tissue in muscle compartment
double VRb     = VR * BVR; 							                                            //Weight/volume of capillary blood in rest of body compartment
double VRt     = VR - VRb; 							                                            //Weight/volume of tissue in rest of body compartment
double VTb     = VT * BVT; 							                                            //Weight/volume of capillary blood in tumor compartment
double VTt     = VT - VTb; 							                                            //Weight/volume of tissue in tumor compartment

//Permeability coefficient-surface area cross-product (L/h)
double PAL     = PALC  * QL; 						                                            //L/h, Liver
double PABr    = PABrC * QBr; 						                                          //L/h, Brain
double PAK     = PAKC  * QK; 						                                            //L/h, Kidneys
double PAS     = PASC  * QS; 						                                          //L/h, Spleen
double PALu    = PALuC * QC;  						                                          //L/h, Lungs
double PAM     = PAMC  * QM; 						                                            //L/h, Muscle
double PAR     = PARC  * QR; 						                                            //L/h, Rest of body
double PAT     = PATC  * QT; 						                                            //L/h, Tumor

//Endocytosis rate (1/h)
double KLRESUP = ((KLRESmax)*pow(TIME, KLRESn))/(pow(KLRES50, KLRESn) + pow(TIME, KLRESn)); 	    
double KSRESUP = ((KSRESmax)*pow(TIME, KSRESn))/(pow(KSRES50, KSRESn) + pow(TIME, KSRESn)); 	    //Lungs
double KKRESUP = ((KKRESmax)*pow(TIME, KKRESn))/(pow(KKRES50, KKRESn) + pow(TIME, KKRESn)); 	    //Lungs
double KLuRESUP = ((KLuRESmax)*pow(TIME, KLuRESn))/(pow(KLuRES50, KLuRESn) + pow(TIME, KLuRESn)); 	    //Lungs
double KMRESUP  = ((KMRESmax)*pow(TIME, KMRESn))/(pow(KMRES50, KMRESn) + pow(TIME, KMRESn));	          //Muscle
double KRRESUP  = ((KRRESmax)*pow(TIME, KRRESn))/(pow(KRRES50, KRRESn) + pow(TIME, KRRESn));		        //Rest of body

double KTRESUP1 = (KTRESmax*pow(TIME,KTRESn))/(pow(KTRES50, KTRESn) + pow(TIME, KTRESn));		              
double KTRESUP2 = KTRESmax*(1-(ATRES/(ATREScap*VT))); 					                        
double KTRESUP3 = 2; 										                                                
double Kbile   = KbileC * pow(BW, 0.75);               
double Kurine  = KurineC * pow(BW, 0.75);

$CMT AA AV ALub ALut ALuRES ABrb ABrt AMb AMt AMRES ARb ARt ARRES 
AKb AKt AKRES Aurine ASb ASt ASRES ALb ALt ALRES Abile ATb ATt ATRES ATRESUP ATRESrel ADOSE 
AUCTumor 

$ODE

double APlasma  = AA + AV;
double ABlood   = AA + AV;
double ALung    = ALub+ALut+ALuRES;
double ALungt   = ALut+ALuRES;
double Arestall = ARb+ARt+ARRES;
double Aresttissue = ARt+ARRES;
double AKidney  = AKb+AKt+AKRES;
double AKidneyt = AKt+AKRES;
double ABrain   = ABrb + ABrt;
double ASpleen  = ASb+ASt+ASRES;
double ASpleent = ASt+ASRES;
double AMuscle  = AMb+AMt+AMRES;
double AMusclet = AMt+AMRES;
double ALiver   = ALb+ALt+ALRES;
double ALivert  = ALt+ALRES;
double ATumor   = ATb+ATt+ATRES;
double ATumort  = ATt+ATRES;


double CA       = AA/(VPlasma * 0.2);
double CV       = AV / (VPlasma * 0.8);
double CPlasma  = APlasma/VPlasma;
double CBlood   = ABlood/VBlood;
double CVLu     = ALub / VLub;
double CVL      = ALb/VLb;
double CLut     = ALut / VLut;
double CLung    = (ALub + ALut + ALuRES)/VLu;
double CLungt   = (ALut+ALuRES)/VLut;
double CVBr     = ABrb/VBrb;
double CBrt     = ABrt/VBrt;
double CBrain   = ABrain / VBr;
double CVK      = AKb/VKb;
double CVM      = AMb/VMb;
double CMt      = AMt/VMt;
double CMuscle  = (AMb+AMt+AMRES)/VM;
double CMusclet = (AMt+AMRES)/VMt;
double CVR      = ARb/VRb;
double CRt      = ARt/VRt;
double Crestall = (ARb+ARt+ARRES)/VR;
double Cresttissue = (ARt+ARRES)/VRt;
double CKt      = AKt/VKt;
double CKidney  = (AKb+AKt+AKRES)/VK;
double CKidneyt = (AKt+AKRES)/VKt;
double CVS      = ASb/VSb;
double CSt      = ASt/VSt;
double CSpleen  = (ASb+ASt+ASRES)/VS;
double CSpleent = (ASt+ASRES)/VSt;
double CLt      = ALt/VLt;
double CLiver = (ALb+ALt+ALRES)/VL;
double CLivert = (ALt+ALRES)/VLt;
double CVT     = ATb/VTb;
double CTt      = ATt/VTt;
double CTumort  = (ATt+ATRES)/VTt;
double CTumor   = (ATb+ATt+ATRES)/VT;


double RA      = QC * CVLu - QC * CA; 
double RV      = QL * CVL + QBr  *CVBr + QK * CVK + QM * CVM + QR * CVR + QT * CVT - QC * CV;
double RLub    =  QC * (CV - CVLu) - PALu * CVLu + (PALu * CLut)/ PLu; 
double RLut    = PALu * CVLu - (PALu * CLut)/ PLu - KLuRESUP * ALut + KLuRESrelease * ALuRES;
double RLuRES  = KLuRESUP * ALut - KLuRESrelease * ALuRES;
double RBrb    = QBr *(CA - CVBr) - PABr * CVBr + (PABr * CBrt)/ PBr;
double RBrt    = PABr * CVBr - (PABr * CBrt)/ PBr;
double RMb      = QM*(CA-CVM) - PAM*CVM + (PAM*CMt)/PM;
double RMt      = PAM*CVM - (PAM*CMt)/PM - KMRESUP*AMt + KMRESrelease*AMRES;
double RMRES    = KMRESUP*AMt-KMRESrelease*AMRES;
double RRb      = QR*(CA-CVR) - PAR*CVR + (PAR*CRt)/PR;
double RRt      = PAR*CVR - (PAR*CRt)/PR - KRRESUP*ARt + KRRESrelease*ARRES;
double RRRES    = KRRESUP*ARt-KRRESrelease*ARRES;
double Rurine   = Kurine*CVK;
double RKb      = QK*(CA-CVK) - PAK*CVK + (PAK*CKt)/PK - Rurine;
double RKt      = PAK*CVK - (PAK*CKt)/PK - KKRESUP*AKt + KKRESrelease*AKRES; 
double RKRES    = KKRESUP*AKt-KKRESrelease*AKRES;
double RSb      = QS*(CA-CVS) - PAS*CVS + (PAS*CSt)/PS; 
double RSt      = PAS*CVS - (PAS*CSt)/PS - KSRESUP*ASt + KSRESrelease*ASRES;
double RSRES    = KSRESUP*ASt-KSRESrelease*ASRES;
double Rbile    = Kbile*CLt ; 
double RLb      = QL*(CA-CVL) + QS*CVS - PAL*CVL + (PAL*CLt)/PL - KLRESUP*ALb + KLRESrelease*ALRES;
double RLt      = PAL*CVL - (PAL*CLt)/PL - Rbile;
double RLRES    = KLRESUP*ALb-KLRESrelease*ALRES;
double RTb      = QT*(CA-CVT) - PAT*CVT + (PAT*CTt)/PT; 
double RTt      = PAT*CVT - (PAT*CTt)/PT - KTRESUP1*ATt + KTRESrelease*ATRES;
double RTRES    = KTRESUP1*ATt-KTRESrelease*ATRES;
double RTRESUP  = KTRESUP1*ATt;
double RTRESrel = KTRESrelease*ATRES;

//Concentration of the chemical in blood compartment
//CA: Arterial blood concentration (mg/L or ug/ml)
dxdt_AA     = RA;
dxdt_AV     = RV;
dxdt_ALub   = RLub;
dxdt_ALut   = RLut;
dxdt_ALuRES = RLuRES;
dxdt_ABrb   = RBrb;
dxdt_ABrt   = RBrt;
dxdt_AMb    = RMb;
dxdt_AMt    = RMt;
dxdt_AMRES  = RMRES;
dxdt_ARb    = RRb;
dxdt_ARt    = RRt;
dxdt_ARRES  = RRRES;
dxdt_AKb    = RKb;
dxdt_AKt    = RKt;
dxdt_AKRES  = RKRES;
dxdt_Aurine = Rurine;
dxdt_Abile  = Rbile;
dxdt_ASb    = RSb;
dxdt_ASt    = RSt;
dxdt_ASRES  = RSRES;
dxdt_ALb    = RLb;
dxdt_ALt    = RLt;
dxdt_ALRES  = RLRES;
dxdt_Abile  = Rbile;
dxdt_ATb    = RTb;
dxdt_ATt         = RTt;
dxdt_ATRES       = RTRES;
dxdt_ATRESUP     = RTRESUP;
dxdt_ATRESrel    = RTRESrel;
dxdt_ADOSE       = 0;
dxdt_AUCTumor    = CTumor;

// {Mass balance equations}
double Tmass = APlasma + ALiver + ABrain + AKidney + ALung + Arestall + AMuscle + ASpleen + Abile + Aurine + ATumor;
double Bal   = ADOSE-Tmass;

$TABLE
capture AUCT     = AUCTumor;
capture Tumor    = ATb+ATt+ATRES;
capture Lung     = CLung;
capture Liver    = CLiver;
capture Kidney   = CKidney;
capture Spleen   = CSpleen;
capture BAL      = Bal;
'


## Build mrgsolve-based PBPK Model
mod <- mcode("NanoPBPK.code", NanoPBPK.code)

# Define the prediction function
pred.nano <- function(pars, pred=FALSE){
  
  # define the dosing parameters
  BW           = pars["BW"]/1000
  PDOSEiv      = pars["Dose"]
  TDOSE        = 1
  tinterval    = 24
  DOSEiv       = PDOSEiv * BW
  
  # define the exposure regimes
  ex.iv.1 <- ev( ID = 1, amt= DOSEiv,                ## Set up the exposure regimes
               ii=tinterval, tinf = 0.005, addl=TDOSE-1, cmt="AV", replicate = FALSE)
  ex.iv.2 <- ev( ID = 1, amt= DOSEiv,                ## Set up the exposure regimes
               ii=tinterval, tinf = 0.005, addl=TDOSE-1, cmt="ADOSE", replicate = FALSE)
  
  ex <- ex.iv.1 + ex.iv.2
  
  # define the time scales
  tsamp  = tgrid(0, tinterval*(TDOSE - 1) + 180, 0.1)
  
  # input parameters, exposure scenario and time into model
  out <- mod %>% param (pars) %>%
                 update(atol = 1E-6, rtol = 1E-3, maxsteps = 50000) %>%
                 mrgsim_d (data = ex, tgrid = tsamp)
  
  outdf = cbind.data.frame (Time    = out$time, 
                            DETumor  = (out$Tumor/DOSEiv)*100)
  
  if(pred) return(outdf)
  
  outdf<-cbind.data.frame(
  DE24  = outdf %>% filter(Time == 24) %>% select(DE24 = DETumor),
  DE168 = outdf %>% filter(Time == 168) %>% select(DE168 = DETumor),
  DEmax = outdf %>% filter(DETumor == max(DETumor)) %>% select(DEmax = DETumor)
  )
  
  
  return(outdf)
}

## Figure 3 -------------------------------------------------------------------
# read the predicted parameters
pre_pars <- read.csv(file = 'Output_pars.csv',fileEncoding = 'UTF-8-BOM')
data <- read.csv(file = 'Data.csv',fileEncoding = 'UTF-8-BOM')


pred_pars <- pre_pars%>%select(
  KTRES_n   = Pred_KTRESn, 
  KTRES_max = Pred_KTRESmax, 
  KTRES_50  = Pred_KTRES50,
  KTRES_rel  = Pred_KTRESrelease)

obs_pars <- pre_pars%>%select( 
  KTRES_n   = KTRESn, 
  KTRES_max = KTRESmax, 
  KTRES_50  = KTRES50,
  KTRES_rel  = KTRESrel) 

df_pred_pars<-pred_pars %>% gather() %>% mutate(group = 'pred')
df_obs_pars<-obs_pars %>% gather() %>% mutate(group = 'obs')

a <- rbind.data.frame(df_pred_pars,df_obs_pars)
a$group <- as.factor(a$group)


p1 <- ggplot(a, aes(x=log(value),fill = group)) +
  geom_histogram(aes(y = ..density..))+
  scale_fill_viridis(discrete=TRUE)   +
  scale_color_viridis(discrete=TRUE)  +
  scale_y_continuous(limit = c(0,1))  +
  xlab("") + ylab("") +
  facet_wrap(~key, scales='free_x')



p2<- p1+theme_bw() + 
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.ticks = element_line(colour = "black", linewidth = 0.2),
    axis.text = element_text(color='black',size=12,face='bold'),
    legend.position = 'none',
    strip.text = element_text(color='black',size=12,face='bold'),
    panel.grid.major = element_blank(),#element_line(colour = "grey70", size = 0.2),
    panel.grid.minor = element_blank()
  )

# Save the figure 4 as tiff-format
# ggsave("Fig 3.tiff", 
#        plot = p2,
#        width = 35, height = 20, units = "cm", dpi=320)




df_pred_pars %>% group_by(key) %>% 
                 summarise(median = median(value),
                           Five_th = quantile(value, 0.025),
                           ninety5_th = quantile(value, 0.975))



df_obs_pars %>% group_by(key) %>% 
                summarise(median = median(value),
                          Five_th = quantile(value, 0.025),
                          ninety5_th = quantile(value, 0.975))


summary(lm(KTRESrel ~Pred_KTRESrelease, data=pre_pars))
summary(lm(KTRES50 ~Pred_KTRES50,data=pre_pars))
summary(lm(KTRESmax ~Pred_KTRESmax,data=pre_pars))
summary(lm(KTRESn ~Pred_KTRESn,data=pre_pars))



## Figure 4-------------------------------------------------------------------
Dat <- data %>% select(ID=ID, 
                       BW = BW, 
                       Dose = Dose, 
                       QTC=QTC,
                       VTC=VTC,
                       BVT=BVT,
                       PT=PT,
                       PATC=PATC)



df_pars <- pre_pars%>%select(ID = ID, 
                             KTRESn   = Pred_KTRESn, 
                             KTRESmax = Pred_KTRESmax, 
                             KTRES50  = Pred_KTRES50,
                             KTRESrelease = Pred_KTRESrelease)

# merge all required parameters
m_df_pars <- merge (df_pars, Dat, by=c("ID"))


# crate a empty list
results <- list()

for (i in 1:dim(m_df_pars)[1]) {
  
  pars <- c(
  BW           = m_df_pars[i,6]/1000, # g to kg
  Dose         = m_df_pars[i,7],
  KTRESn       = m_df_pars[i,2],
  KTRESmax     = m_df_pars[i,3],
  KTRES50      = m_df_pars[i,4],
  KTRESrelease = ifelse(m_df_pars[i,5]<0, 0.005, m_df_pars[i,5]), # replace the negative value
  QTC          = m_df_pars[i,8],
  VTC          = m_df_pars[i,9],
  BVT          = m_df_pars[i,10],
  PT           = m_df_pars[i,11],
  PATC         = m_df_pars[i,12])
  
  results[[i]] = pred.nano(pars)

}

pred <- do.call(rbind.data.frame,results)
pred <- cbind.data.frame(ID = pre_pars[,1], pred)

ObsDat <- data %>% select(ID = ID, 
                          DE24obs = DE24, 
                          DE168obs = DE168, 
                          DEmaxobs = Demax)

# create a function to calculate the r-square
rsq <- function (x, y) cor(x, y) ^ 2


a<-merge(merge (pred, ObsDat, by=c("ID")), data %>% select(ID, Type, Shape), by=c("ID"))

summary(a)
summary(data)

## Plot for model calibration
# Install
## Plot theme for ggplot
ptheme<-theme (
  plot.background         = element_rect (fill="White"),
  text                    = element_text (family = "Times"),   # text front (Time new roman)
  panel.border            = element_rect (colour = "black", fill=NA, linewidth=2),
  panel.background        = element_rect (fill="White"),
  panel.grid.major        = element_blank(),
  panel.grid.minor        = element_blank(),
  axis.text               = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes
  axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
  legend.position='none') 

p1 <- 
  ggplot(a, aes(log10(DE24), log10(DE24obs))) + 
  geom_point  (aes(colour = factor(Shape),  shape = factor(Type)), size = 4)+
  scale_colour_brewer(palette = 1) +
  #scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  #scale_shape_manual(values=seq(16,22))+
  geom_abline (intercept = 0, 
               slope     = 1,
               color     ="black", linewidth = 1, alpha = 0.8, linetype = 2) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-2,2), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-2,2),labels = scales::math_format(10^.x)) +
  ptheme + labs (x = "", y = "")

#p1

model <- lm(log10(DE24obs) ~ log10(DE24), data=a)
summary(model)

p2 <- 
  ggplot(a, aes(log10(DE168), log10(DE168obs))) + 
  geom_point  (aes(colour = factor(Shape),  shape = factor(Type)), size = 4)+
  scale_colour_brewer(palette = 2) +
  #scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  #scale_shape_manual(values=seq(16,22))+
  geom_abline (intercept = 0, 
               slope     = 1,
               color     ="black", size = 1, alpha = 0.8, linetype = 2) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-2,2), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-2,2),labels = scales::math_format(10^.x)) +
  ptheme + labs (x = "", y = "")


#p2

model <- lm(log10(DE168obs) ~ log10(DE168), data=a)
summary(model)



p3 <- 
  ggplot(a, aes(log10(DEmax), log10(DEmaxobs))) + 
  geom_point  (aes(colour = factor(Shape),  shape = factor(Type)), size = 4)+
  scale_colour_brewer(palette = 3) +
  #scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  #scale_shape_manual(values=seq(16,22))+
  geom_abline (intercept = 0, 
               slope     = 1,
               color     ="black", size = 1, alpha = 0.8, linetype = 2) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-2,2), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-2,2),labels = scales::math_format(10^.x)) +
  ptheme + labs (x = "", y = "")


#p3

model <- lm(log10(DEmax) ~ log10(DEmaxobs), data=a)
summary(model)


# save the figure as tiff-format
# ggsave("Fig 4.tiff", 
#        plot = p1+p2+p3,
#        width = 35, height = 10, units = "cm", dpi=320)
# 

## 2- or 3- fold error
# DE24
sum((a$DE24)/(a$DE24obs)<=2 & (a$DE24)/(a$DE24obs)>=0.5)/288
sum((a$DE24)/(a$DE24obs)<=3 & (a$DE24)/(a$DE24obs)>=0.3)/288
rmse(a$DE24obs,a$DE24) 

# DE168
sum((a$DE168)/(a$DE168obs)<=2 & (a$DE168)/(a$DE168obs)>=0.5)/288
sum((a$DE168)/(a$DE168obs)<=3 & (a$DE168)/(a$DE168obs)>=0.3)/288
rmse(a$DE168obs,a$DE168) 

# DEMAX
sum((a$DEmax)/(a$DEmaxobs)<=2 & (a$DEmax)/(a$DEmaxobs)>=0.5)/288
sum((a$DEmax)/(a$DEmaxobs)<=3 & (a$DEmax)/(a$DEmaxobs)>=0.3)/288
rmse(a$DEmaxobs,a$DEmax) 



##--Figure 5---------------------------------------------------------
## -----Time vs. conc. comparison
TvsC_Dat <- read.csv(file='EvaDat_1.csv',fileEncoding = 'UTF-8-BOM')
TvsC_Dat$ID<-na.locf(TvsC_Dat$ID) # replace na values with last no-na value
colnames(TvsC_Dat) <- c("ID","Time", "DETumor","DETumor_g")

# crate a empty list
TvsC_res <- list()
pred_all <- list()

for (i in 1:dim(m_df_pars)[1]) {
  
  
  ID = m_df_pars[i,1]
  obs_Dat <- TvsC_Dat[TvsC_Dat$ID==ID,] %>% select(Time, DETumor_g)
  TW <- data[data$ID==ID, "TW"]
    
  pars <- c(
    BW           = m_df_pars[i,6]/1000, # g to kg
    Dose         = m_df_pars[i,7],
    KTRESn       = m_df_pars[i,2],
    KTRESmax     = m_df_pars[i,3],
    KTRES50      = m_df_pars[i,4],
    KTRESrelease = ifelse(m_df_pars[i,5]<0, 0.005, m_df_pars[i,5]), # replace the negative value
    QTC          = m_df_pars[i,8],
    VTC          = m_df_pars[i,9],
    BVT          = m_df_pars[i,10],
    PT           = m_df_pars[i,11],
    PATC         = m_df_pars[i,12])
  
  pred_Dat = pred.nano(pars, pred=TRUE)%>%mutate(DETumor_g=DETumor/TW)%>% select(Time, DETumor_g)
  cost <- modCost (model = pred_Dat, obs = obs_Dat, x='Time')
  TvsC_res[[i]]<-cost$residuals %>% mutate(ID=ID)%>%select(ID=ID,Time=x, obs=obs,mod=mod)
  pred_all[[i]]<-pred_Dat%>% mutate(ID=ID) 
  
}


TvsC_pred <- do.call(rbind.data.frame,TvsC_res)
pred_all_2 <- do.call(rbind.data.frame,pred_all)

TvsC_pred2 <- TvsC_pred%>%filter(obs>0 & obs<200)%>%mutate(ratios = obs/mod)

b<-merge(TvsC_pred2, data %>% select(ID, Type, Shape), by=c("ID"))


# Plot
p4 <- 
  ggplot(b, aes(log10(obs), log10(mod))) + 
  geom_point  (aes(colour = factor(Shape),  shape = factor(Type)), size = 4)+
  scale_colour_brewer(palette = 6) +
  #scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  scale_shape_manual(values=seq(15,21))+
  geom_abline (intercept = 0, 
               slope     = 1,
               color     ="black", size = 1, alpha = 0.8, linetype = 2) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-2,3), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-2,3),labels = scales::math_format(10^.x)) +
  ptheme + labs (x = "", y = "")


#p4

model <- lm(log10(obs) ~ log10(mod), data=b)
summary(model)

sum(b$ratios<=2 & b$ratios>=0.5)/length(b$ratios)
sum(b$ratios<=3 & b$ratios>=0.3)/length(b$ratios)

RMSE <- (sum((b$mod - b$obs)^2)/length(b$ratios))^0.5

p4_2 <-
  ggplot(b, aes(log10(mod), log10(ratios))) +
  geom_hline(yintercept = log10(2),linetype = 3,color   = "red", size =1) +
  geom_hline(yintercept = log10(0.5),linetype = 3,color   = "red", size =1) +
  geom_hline(yintercept = log10(3),linetype = 3,color   = "black", size =1) +
  geom_hline(yintercept = log10(0.3),linetype = 3,color   = "black", size =1) +
  geom_point  (aes(colour = factor(Shape),  shape = factor(Type)), size = 4)+
  scale_colour_brewer(palette = 6) +
  scale_shape_manual(values=seq(15,21))+
  #geom_smooth(se = FALSE) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-2,3), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-2,3),labels = scales::math_format(10^.x)) +
  ptheme + labs (x = "", y = "")


ggMarginal(p4_2, type = "histogram", margins = "y",
           yparams = list(binwidth = 0.1, fill = "#FFA488"))


p4 + ggMarginal(p4_2, type = "histogram", margins = "y",
                yparams = list(binwidth = 0.2, fill = "#FFA488")) +
  plot_layout(width=c(1,1.3))


# save figure as tiff-format
# ggsave("Fig 5.tiff", 
#        plot = p4 + ggMarginal(p4_2, type = "histogram", margins = "y",
#                               yparams = list(binwidth = 0.1, fill = "#FFA488"))+
#       plot_layout(width=c(1,1.3)),
#        width = 30, height = 10, units = "cm", dpi=320)


##-Figure 6--------------------------------------------------------------------

c<-b %>%filter(ID <150 & ID >133) #& ratios <=1 & ratios >=0.6) #%>% 
  #filter(!ID %in% c(5,8,10,16,23,24,32,46,52,62,
    #                79,87,100,101,104,109,111,116,124,125,128:131,138, 148,149))
  
d <- pred_all_2 %>% filter(ID %in% unique(c$ID))

p5<-ggplot(data=c) +
  geom_point(aes(Time, obs,color=factor(Shape),size=2)) +
  geom_line(data=d, aes(Time, (DETumor_g)),linetype="twodash", size=1)+
  scale_x_continuous(limits=c(0,24))+
  facet_wrap(~ID, scales='free_x', ncol=3)+
  labs (x = "", y = "") +
  theme (
    plot.background         = element_rect (fill="White"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    panel.border            = element_rect (colour = "black", fill=NA, linewidth=1),
    panel.background        = element_rect (fill="White"),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(),
    axis.text               = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes
    axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
    legend.position='none', 
    strip.background = element_blank(),
    strip.text.x = element_blank()
    ) 

# save the figure as tiff-format
# ggsave("Fig 6.tiff", 
#        plot = p5, 
#        width = 30, height = 20, units = "cm", dpi=320)
# 

# estimate the r-sqaured value
for (i in unique(c$ID)) {
  c1 <- c%>%filter(ID==i)
  r2 <- (cor(c1$obs, c1$mod))**2
  #lm(obs ~ mod, data=c1)
  #print(summary(model)$adj.r.squared)
  print(r2)
  }

  
##------------------------------------------------------------------------------
## Figures in Supplementary materials
for (i in 1:25) {
j = seq(1, 378, 15)[i]

c<- b %>% filter(ID %in% unique(b$ID)[seq(j,j+14)]) #& ratios <=1 & ratios >=0.6) #%>% 

#filter(!ID %in% c(5,8,10,16,23,24,32,46,52,62,
#                79,87,100,101,104,109,111,116,124,125,128:131,138, 148,149))

d <- pred_all_2 %>% filter(ID %in% unique(c$ID))

p5<-ggplot(data=c) +
  geom_point(aes(Time, obs,color=factor(Shape),size=2)) +
  geom_line(data=d, aes(Time, (DETumor_g)),linetype="twodash", size=1)+
  scale_x_continuous(limits=c(0,24))+
  facet_wrap(~ID, scales='free_x', ncol=3)+
  labs (x = "", y = "") +
  theme (
    plot.background         = element_rect (fill="White"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    panel.border            = element_rect (colour = "black", fill=NA, linewidth=1),
    panel.background        = element_rect (fill="White"),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(),
    axis.text               = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes
    axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
    legend.position='none', 
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) 

ggsave(paste("FigS", i,".tiff"), 
       plot = p5, 
       #plot_layout(width=c(1,1.3)),
       width = 30, height = 20, units = "cm", dpi=320)

}



## Supplementary materials
for (i in 1:25) {
  j = seq(1, 378, 15)[i]
  
  c<- b %>% filter(ID %in% unique(b$ID)[seq(j,j+14)]) #& ratios <=1 & ratios >=0.6) #%>% 

  print(unique(c$ID))
  }


## estimate the r-squared value 

adj_r <- list()

for (i in 1:length(unique(b$ID))) {
  dat_r <- b %>% filter(ID %in% unique(b$ID)[i])
  adj_r[[i]]<-summary(lm(obs ~mod, data=dat_r))$adj.r.squared
}
  
adj_R<-do.call(rbind,adj_r)

length(adj_R[adj_R>=0.7])





