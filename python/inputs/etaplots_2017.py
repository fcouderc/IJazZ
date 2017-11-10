import ROOT as rt
import ijazzClass 

dirEtaSc = 'IJazZ_output/FitResultsClassic/ecalEtaScale/0-999999/'
inputs = {
    'EtaScale_Cal_Oct2017_cand_v1' : {
        'files' : [
            ijazzClass.ijazzResOP( dirEtaSc + 'IJazZ_EB_MC_Cal_Oct2017_cand_v1-Z_RunBCDEF.root' , color = rt.kRed+1   , legend = 'MC (RunBCDEF)', 
                                   isRef = True ),
            ijazzClass.ijazzResOP( dirEtaSc + 'IJazZ_EB_Data_Cal_Oct2017_cand_v1-Z_RunB.root', color = rt.kBlack   , legend = 'RunB', lumi = 4.890 ),
            ijazzClass.ijazzResOP( dirEtaSc + 'IJazZ_EB_Data_Cal_Oct2017_cand_v1-Z_RunC.root', color = rt.kAzure+1 , legend = 'RunC', lumi = 9.932 ),
            ijazzClass.ijazzResOP( dirEtaSc + 'IJazZ_EB_Data_Cal_Oct2017_cand_v1-Z_RunD.root', color = rt.kGreen+2 , legend = 'RunD', lumi = 4.364 ),
            ijazzClass.ijazzResOP( dirEtaSc + 'IJazZ_EB_Data_Cal_Oct2017_cand_v1-Z_RunE.root', color = rt.kBlue-2  , legend = 'RunE', lumi = 9.498 ),
            ijazzClass.ijazzResOP( dirEtaSc + 'IJazZ_EB_Data_Cal_Oct2017_cand_v1-Z_RunF.root', color = rt.kYellow+1, legend = 'RunF', lumi = 4.690 ),
            ],        
        },
    'EtaScale_Cal_Oct2017_MC' : {
        'files' : [
            ijazzClass.ijazzResOP( dirEtaSc + 'IJazZ_EB_MC_Cal_Oct2017_cand_v1-Z_RunB.root', color = rt.kBlack   , legend = 'RunB', lumi = 4.890 ),
            ijazzClass.ijazzResOP( dirEtaSc + 'IJazZ_EB_MC_Cal_Oct2017_cand_v1-Z_RunC.root', color = rt.kAzure+1 , legend = 'RunC', lumi = 9.932 ),
            ijazzClass.ijazzResOP( dirEtaSc + 'IJazZ_EB_MC_Cal_Oct2017_cand_v1-Z_RunD.root', color = rt.kGreen+2 , legend = 'RunD', lumi = 4.364 ),
            ijazzClass.ijazzResOP( dirEtaSc + 'IJazZ_EB_MC_Cal_Oct2017_cand_v1-Z_RunE.root', color = rt.kBlue-2  , legend = 'RunE', lumi = 9.498 ),
            ijazzClass.ijazzResOP( dirEtaSc + 'IJazZ_EB_MC_Cal_Oct2017_cand_v1-Z_RunF.root', color = rt.kYellow+1, legend = 'RunF', lumi = 4.690 ),
            ],        
        },


    'EtaScale_Cal_Oct2017_RunB' : {
        'files' : [
            ijazzClass.ijazzResOP( dirEtaSc + 'IJazZ_EB_Data_Cal_Oct2017_cand_v1-Z_RunB.root', color = rt.kBlack   , legend = 'RunB', lumi = 4.890 ),
            ijazzClass.ijazzResOP( dirEtaSc + 'IJazZ_EB_Data_Cal_Oct2017_cand_v1-Z_RunB.root', color = rt.kAzure+1 , legend = 'RunC', lumi = 4.890 ),
            ijazzClass.ijazzResOP( dirEtaSc + 'IJazZ_EB_Data_Cal_Oct2017_cand_v1-Z_RunB.root', color = rt.kGreen+2 , legend = 'RunD', lumi = 4.890 ),
            ijazzClass.ijazzResOP( dirEtaSc + 'IJazZ_EB_Data_Cal_Oct2017_cand_v1-Z_RunB.root', color = rt.kBlue-2  , legend = 'RunE', lumi = 4.890 ),
            ijazzClass.ijazzResOP( dirEtaSc + 'IJazZ_EB_Data_Cal_Oct2017_cand_v1-Z_RunB.root', color = rt.kYellow+1, legend = 'RunF', lumi = 4.690 ),
            ],        
        },
    }
