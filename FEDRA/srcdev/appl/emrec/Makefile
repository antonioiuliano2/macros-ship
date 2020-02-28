#------------------------------------------------------------------------------
include $(FEDRA_ROOT)/src/config/ProjectDef.mk
include $(FEDRA_ROOT)/src/config/RootDef.mk
#------------------------------------------------------------------------------

TARGET	=	$(BIN_DIR)/emlink$(ExeSuf)      \
		$(BIN_DIR)/emalign$(ExeSuf)      \
		$(BIN_DIR)/emaligns2s$(ExeSuf)      \
		$(BIN_DIR)/emalignraw$(ExeSuf)      \
		$(BIN_DIR)/makescanset$(ExeSuf)      \
		$(BIN_DIR)/emtra$(ExeSuf)      \
		$(BIN_DIR)/emshow$(ExeSuf)      \
		$(BIN_DIR)/emrec$(ExeSuf)      \
		$(BIN_DIR)/emvtx$(ExeSuf)      \
		$(BIN_DIR)/emcheck$(ExeSuf)    \
		$(BIN_DIR)/empred$(ExeSuf)    \
		$(BIN_DIR)/emrawcorr$(ExeSuf)    \
		$(BIN_DIR)/emtrackan$(ExeSuf)
#		$(BIN_DIR)/emtraceback$(ExeSuf) \
#		$(BIN_DIR)/emtrfind$(ExeSuf)
#		$(BIN_DIR)/emcpan$(ExeSuf)



ROOTLIBS      += -lMLP
##PROJECT_LIBS += -L$(FEDRA_ROOT)/lib -lScan -lEdb -lEbase -lEmath -lEphys -lvt -lEdr -lScan -lEdr -lEIO -lEMC -lAlignment -lDataConversion -lShower -lTreePlayer
PROJECT_LIBS += $(GLIBS) -L$(FEDRA_ROOT)/lib -lScan -lEdb -lEbase -lEmath -lEphys -lvt -lEdr -lScan -lEdr -lEIO -lEMC -lAlignment -lDataConversion -lEdd -lTreePlayer -lScan -lEdb -lEbase -lEmath -lEphys -lvt -lEdr -lScan -lShower -lSpectrum -lEdr -lEIO -lEMC -lAlignment -lDataConversion -lEdd -lTreePlayer $(ROOTLIBS) $(GLIBS)

#------------------------------------------------------------------------------
include ../../config/TargetsDef.mk
#--------------------------------------------------------------------------
