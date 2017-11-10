# Source files to compile
FILES := ZFitterMinuit2_ND_MThread
FILES += IJazZ_fitIO
#FILES += IJazZ_tupleVar_sEoP
FILES += IJazZ_tupleVar
FILES += IJazZ_eventSelectionND
FILES += IJazZ_utils
FILES += RootUtils
FILES += EcalUtils
FILES += EcalUtilsEtaScale
FILES += EnergyCorrectionAndSmearingHgg

# Header files to use for dictionary generation
# In this module: is the same as FILES
DICTFILES :=

# Executable files
PROGRAMS := IJazZexe
PROGRAMS += IJazZetaScale
#PROGRAMS += FabAxisTest

NEEDS_ROOT  := yes
NEEDS_BOOST := yes
