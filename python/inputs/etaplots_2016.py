import ROOT as rt
import ijazzClass 

dirEtaScale2016 = 'IJazZ_output/FitResultsClassic/ecalEtaScale/0-999999/'
inputs = {
    'EtaScale_Cal_March2017_MC' : {
        'files' : [
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_MC_Cal_March2017_RunB.root'  , color = rt.kRed+1   , legend = 'MC (RunB)' , isMC = True, lumi = 5.77 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_MC_Cal_March2017_RunCD.root' , color = rt.kBlack   , legend = 'MC (RunCD)', isMC = True, lumi = 7.00 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_MC_Cal_March2017_RunEF.root' , color = rt.kAzure+1 , legend = 'MC (RunEF)', isMC = True, lumi = 7.15 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_MC_Cal_March2017_RunG.root'  , color = rt.kGreen+1 , legend = 'MC (RunG)' , isMC = True, lumi = 7.54 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_MC_Cal_March2017_RunH.root'  , color = rt.kYellow+1, legend = 'MC (RunH)' , isMC = True, lumi = 8.76 ),
            ],        
        },
    'EtaScale_Cal_March2017_ref' : {
        'files' : [
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_MC_Cal_March2017_RunB.root'   , color = rt.kRed+1   , legend = 'MC (RunB)', isRef = True ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_March2017_RunB.root' , color = rt.kBlack   , legend = 'RunB  ref', lumi = 5.77 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_March2017_RunCD.root', color = rt.kAzure+1 , legend = 'RunCD ref', lumi = 7.00 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_March2017_RunEF.root', color = rt.kGreen+2 , legend = 'RunEF ref', lumi = 7.15 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_March2017_RunG.root' , color = rt.kBlue-2  , legend = 'RunG  ref', lumi = 7.54 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_March2017_RunH.root' , color = rt.kYellow+1, legend = 'RunH  ref', lumi = 8.76 ),
            ],        
        },

    'EtaScale_Cal_March2017_ref_v2' : {
        'files' : [
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_MC_Cal_March2017_RunB.root'           , color = rt.kRed+1   , legend = 'MC (RunB)', isRef = True ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ref_v2_RunB.root' , color = rt.kBlack   , legend = 'RunB  ref v2', lumi = 5.77 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ref_v2_RunCD.root', color = rt.kAzure+1 , legend = 'RunCD ref v2', lumi = 7.00 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ref_v2_RunEF.root', color = rt.kGreen+2 , legend = 'RunEF ref v2', lumi = 7.15 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ref_v2_RunG.root' , color = rt.kBlue-2  , legend = 'RunG  ref v2', lumi = 7.54 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ref_v2_RunH.root' , color = rt.kYellow+1, legend = 'RunH  ref v2', lumi = 8.76 ),
            ],        
        },
    'EtaScale_Cal_March2017_ICcombv2' : {
        'files' : [
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_MC_Cal_March2017_RunB.root'           , color = rt.kRed+1   , legend = 'MC (RunB)', isRef = True ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_v2_RunB.root' , color = rt.kBlack   , legend = 'RunB  IComb v2', lumi = 5.77 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_v2_RunCD.root', color = rt.kAzure+1 , legend = 'RunCD IComb v2', lumi = 7.00 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_v2_RunEF.root', color = rt.kGreen+2 , legend = 'RunEF IComb v2', lumi = 7.15 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_v2_RunG.root' , color = rt.kBlue-2  , legend = 'RunG  IComb v2', lumi = 7.54 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_v2_RunH.root' , color = rt.kYellow+1, legend = 'RunH  IComb v2', lumi = 8.76 ),
            ],        
        },
        
    'EtaScale_Cal_March2017_ICcombv2_etascale' : {
        'files' : [
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_MC_Cal_March2017_RunB.root'           , color = rt.kRed+1   , legend = 'MC (RunB)', isRef = True ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_etaScale_v2_RunB.root' , color = rt.kBlack   , legend = 'RunB  IComb v2 (scale)', lumi = 5.77 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_etaScale_v2_RunCD.root', color = rt.kAzure+1 , legend = 'RunCD IComb v2 (scale)', lumi = 7.00 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_etaScale_v2_RunEF.root', color = rt.kGreen+2 , legend = 'RunEF IComb v2 (scale)', lumi = 7.15 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_etaScale_v2_RunG.root' , color = rt.kBlue-2  , legend = 'RunG  IComb v2 (scale)', lumi = 7.54 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_etaScale_v2_RunH.root' , color = rt.kYellow+1, legend = 'RunH  IComb v2 (scale)', lumi = 8.76 ),
            ],        
        },
        'EtaScale_Cal_March2017_ICcombv3_etascale' : {
        'files' : [
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_MC_Cal_March2017_RunB.root'           , color = rt.kRed+1   , legend = 'MC (RunB)', isRef = True ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_etaScale_v3_RunB.root' , color = rt.kBlack   , legend = 'RunB  IComb v3 (scale)', lumi = 5.77 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_etaScale_v3_RunCD.root', color = rt.kAzure+1 , legend = 'RunCD IComb v3 (scale)', lumi = 7.00 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_etaScale_v3_RunEF.root', color = rt.kGreen+2 , legend = 'RunEF IComb v3 (scale)', lumi = 7.15 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_etaScale_v3_RunG.root' , color = rt.kBlue-2  , legend = 'RunG  IComb v3 (scale)', lumi = 7.54 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_etaScale_v3_RunH.root' , color = rt.kYellow+1, legend = 'RunH  IComb v3 (scale)', lumi = 8.76 ),
            ],        
        },
    'EtaScale_Cal_March2017_ICcombv5' : {
        'files' : [
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_MC_Cal_March2017_RunB.root'           , color = rt.kRed+1   , legend = 'MC (RunB)', isRef = True ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_v5_RunB.root' , color = rt.kBlack   , legend = 'RunB  IComb v5', lumi = 5.77 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_v5_RunCD.root', color = rt.kAzure+1 , legend = 'RunCD IComb v5', lumi = 7.00 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_v5_RunEF.root', color = rt.kGreen+2 , legend = 'RunEF IComb v5', lumi = 7.15 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_v5_RunG.root' , color = rt.kBlue-2  , legend = 'RunG  IComb v5', lumi = 7.54 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_v5_RunH.root' , color = rt.kYellow+1, legend = 'RunH  IComb v5', lumi = 8.76 ),
            ],        
        },
        
    'EtaScale_Cal_March2017_ICcombv6' : {
        'files' : [
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_MC_Cal_March2017_RunB.root'           , color = rt.kRed+1   , legend = 'MC (RunB)', isRef = True ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_v6_RunB.root' , color = rt.kBlack   , legend = 'RunB  IComb v6', lumi = 5.77 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_v6_RunCD.root', color = rt.kAzure+1 , legend = 'RunCD IComb v6', lumi = 7.00 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_v6_RunEF.root', color = rt.kGreen+2 , legend = 'RunEF IComb v6', lumi = 7.15 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_v6_RunG.root' , color = rt.kBlue-2  , legend = 'RunG  IComb v6', lumi = 7.54 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_v6_RunH.root' , color = rt.kYellow+1, legend = 'RunH  IComb v6', lumi = 8.76 ),
            ],        
        },

    ### buggy versions 
    'EtaScale_Cal_March2017_ICcombv1' : {
        'files' : [
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_MC_Cal_March2017_RunB.root'             , color = rt.kRed+1   , legend = 'MC (RunB)', isRef = True ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_March2017_RunB_ICcomb_v1.root' , color = rt.kBlack   , legend = 'RunB  IComb v1', lumi = 5.77 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_March2017_RunCD_ICcomb_v1.root', color = rt.kAzure+1 , legend = 'RunCD IComb v1', lumi = 7.00 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_March2017_RunEF_ICcomb_v1.root', color = rt.kGreen+2 , legend = 'RunEF IComb v1', lumi = 7.15 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_March2017_RunG_ICcomb_v1.root' , color = rt.kBlue-2  , legend = 'RunG  IComb v1', lumi = 7.54 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_March2017_RunH_ICcomb_v1.root' , color = rt.kYellow+1, legend = 'RunH  IComb v1', lumi = 8.76 ),
            ],        
        },
    'EtaScale_Cal_March2017_ICcombv3' : {
        'files' : [
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_MC_Cal_March2017_RunB.root'             , color = rt.kRed+1   , legend = 'MC (RunB)', isRef = True ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_v3_RunB.root' , color = rt.kBlack   , legend = 'RunB  IComb v3', lumi = 5.77 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_v3_RunCD.root', color = rt.kAzure+1 , legend = 'RunCD IComb v3', lumi = 7.00 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_v3_RunEF.root', color = rt.kGreen+2 , legend = 'RunEF IComb v3', lumi = 7.15 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_v3_RunG.root' , color = rt.kBlue-2  , legend = 'RunG  IComb v3', lumi = 7.54 ),
            ijazzClass.ijazzResOP( dirEtaScale2016 + 'IJazZ_EB_Data_Cal_Mar2017_ICcomb_v3_RunH.root' , color = rt.kYellow+1, legend = 'RunH  IComb v3', lumi = 8.76 ),
            ],        
        },

    }
