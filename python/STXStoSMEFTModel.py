# Author: Jonathon Langford (ICL)
# Date: 10/2021
# Description: Model to describe how bins in STXS stage 1.2 scale using full set of dimension-6 EFT parameters
#              Equations calculated using nanoAOD reweighting using Madgraph reweighting modules
#              SMEFTsim

from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder
from math import exp
import ROOT, os, re, sys
import json
import yaml
from collections import OrderedDict as od

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Global function to extract reco category, STXS bin, decay mode and energy from process name
def getProcessInfo(bin,process):
  foundRecoCategory = bin
  foundSTXSBin = process
  foundDecay = None
  foundEnergy = "13TeV"
  #Iterate over Higgs decays
  matchedDecayString = False 
  for D in ALL_HIGGS_DECAYS:
    if matchedDecayString: continue
    if "_%s"%D in foundSTXSBin:
      foundSTXSBin = re.sub('_%s'%D,'',foundSTXSBin)
      foundDecay = D
      matchedDecayString = True
  # Also drop year tag in STXS bin name if present
  for Y in ['2016','2017','2018']:
    if "_%s"%Y in foundSTXSBin:
      foundSTXSBin = re.sub('_%s'%Y,'',foundSTXSBin)

  if not matchedDecayString: raise RuntimeError, "Validation Error: no supported decay found in process"

  return (foundRecoCategory, foundSTXSBin, foundDecay, foundEnergy)

#################################################################################################################
# STXS to EFT abstract base class: inherited classes for different stages
class STXStoSMEFTBaseModel(SMLikeHiggsModel):

  def __init__(self,fixProcesses=[]):
    SMLikeHiggsModel.__init__(self)
    self.PROCESSES = None
    self.DECAYS = None
    # Dicts to store pois + scaling functions
    self.pois = None
    self.STXSScalingTerms = None
    self.DecayScalingTerms = None
    # Options
    self.floatMass = False
    self.fixProcesses = fixProcesses #Option to fix certain STXS bins: comma separated list of STXS bins
    self.useExtendedQQHScheme=False
    self.ignoreAccCorr=False # Decides which scaling terms to load
    self.linearOnly=False

  def setPhysicsOptionsBase(self,physOptions):
    for po in physOptions:
      if po.startswith("higgsMassRange="):
        self.floatMass = True
        self.mHRange = po.replace("higgsMassRange=","").split(",")
        if len(self.mHRange) != 2:
          raise RuntimeError, "Higgs mass range definition requires two extrema"
        elif float(self.mHRange[0]) >= float(self.mHRange[1]):
          raise RuntimeError, "Extrema for Higgs mass range defined with inverterd order. Second must be larger the first"
      if po.startswith("fixProcesses="): 
        self.fixProcesses = (po.replace("fixProcesses=","")).split(",")
      if po.startswith("useExtendedQQHScheme="): 
        self.useExtendedQQHScheme = (po.replace("useExtendedQQHScheme=","") in ["yes","1","Yes","True","true"])
      if po.startswith("ignoreAccCorr="): 
        self.ignoreAccCorr = (po.replace("ignoreAccCorr=","") in ["yes","1","Yes","True","true"])
      if po.startswith("linearOnly="): 
        self.linearOnly = (po.replace("linearOnly=","") in ["yes","1","Yes","True","true"])

    #Output options to screen
    if( len( self.fixProcesses ) > 0 ): print " --> [STXStoSMEFT] Fixing following processes to SM: %s"%self.fixProcesses
    if self.useExtendedQQHScheme: print " --> [STXStoSMEFT] Use extended qqH scheme (different scaling for VBF + V(qq)H)"
    if self.ignoreAccCorr: print " --> [STXStoSMEFT] Will use inclusive (per-STXS bin) scaling functions"
    if self.linearOnly: print " --> [STXStoSMEFT] Only linear terms (Aj)"

  def doMH(self):
    if self.floatMass:
      if self.modelBuilder.out.var("MH"):
        self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
        self.modelBuilder.out.var("MH").setConstant(False)
      else:
        self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1])) 
    else:
      if self.modelBuilder.out.var("MH"):
        self.modelBuilder.out.var("MH").setVal(self.options.mass)
        self.modelBuilder.out.var("MH").setConstant(True)
      else:
        self.modelBuilder.doVar("MH[%g]" % self.options.mass)

  # Overwrite getYieldScale to extract (RECO-category,STXS bin,decay,energy)
  def getYieldScale(self,bin,process):
    if not self.DC.isSignal[process]: 
      return 1.

    # Extract process line info
    (recocat, stxsbin, decay, energy) = getProcessInfo(bin,process)

    # Return 1 (no scaling) for fixed processes and scaling for non-fixed
    if stxsbin in self.fixProcesses: 
      return 1. 
    else: 
      if not self.ignoreAccCorr:
        procStr = "%s__%s"%(stxsbin,recocat)
      else:
        procStr = stxsbin
      return self.getHiggsSignalYieldScale(procStr, decay, energy)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Extract pois from yaml file
  def extractPOIs( self, filename ):
    with open( filename, 'r' ) as fpois:
      try:
        self.pois = yaml.safe_load(fpois)
      except yaml.YAMLERROR as exc:
        print exc

  #Function to extract STXS scaling terms from json file
  def extractSTXSScalingTerms( self, filename="" ):
    if filename != "":
      with open(filename,"r") as jf: self.STXSScalingTerms = json.load(jf) 
    else: 
      self.STXSScalingTerms = {}
  
  #Function to extract decay scaling functions from file
  def extractDecayScalingTerms( self, filename="" ):
    if filename != "":
      with open(filename,"r") as jf: self.DecayScalingTerms = json.load(jf)
    else:
      self.DecayScalingTerms = {}
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Function to make scaling function in workspace from terms
  # TODO: make a RooEFTScalingFunction compiled class that takes A and B parameters as input
  def makeScalingFunction( self, what, isDecay = False ):

    # Extract terms for dict
    if what in self.STXSScalingTerms: terms = self.STXSScalingTerms[what] 
    elif what in self.DecayScalingTerms: terms = self.DecayScalingTerms[ what ]
    else:
      print " --> [WARNING] Scaling terms for %s do not exist. Setting to 1"%what
      terms = {}
      #raise ValueError("[ERROR] Scaling terms for %s do not exist"%what)

    # Loop over pois and extract the terms from scaling function
    A, B = od(), od()
    for jpoi in self.pois:
      # Interference terms: Aj
      if "A_%s"%jpoi in terms: A[jpoi] = terms["A_%s"%jpoi]
      # BSM-only terms: Bjk
      if not self.linearOnly:
	if "B_%s_2"%jpoi in terms: B[(jpoi,jpoi)] = terms["B_%s_2"%jpoi]
	# Cross terms
	for kpoi in self.pois:
	  if "B_%s_%s"%(jpoi,kpoi) in terms: B[(jpoi,kpoi)] = terms["B_%s_%s"%(jpoi,kpoi)]
        
    # TODO: Make RooEFTScalingFunction
    #self.modelBuilder.factory_("RooEFTScalingFunction::scaling_%s(A,B)"%what)
    # Replace below code
    termIdx = 0
    termList = []
    for poi,coeff in A.iteritems():
      termStr = "prod::term_%s_%g(%g,%g,%s)"%(what,termIdx,coeff,self.pois[poi]['multiplier'],poi)
      self.modelBuilder.factory_( termStr )
      termList.append("term_%s_%g"%(what,termIdx))
      termIdx += 1

    for pois,coeff in B.iteritems():
      jpoi, kpoi = pois[0], pois[1]
      termStr = "prod::term_%s_%g(%g,%g,%s,%g,%s)"%(what,termIdx,coeff,self.pois[jpoi]['multiplier'],jpoi,self.pois[kpoi]['multiplier'],kpoi)
      self.modelBuilder.factory_( termStr )
      termList.append("term_%s_%g"%(what,termIdx))
      termIdx += 1
 
    # Split up terms into sizeable chunks to sum: 15 terms max
    sumTerms = od()
    sumIdx = -1
    for termIdx in range(len(termList)): 
      if termIdx % 10 == 0: sumIdx +=1
      #Add terms to sum
      if( what in self.DecayScalingTerms )&( what != "tot" ): sumString = "scaling_partial_%s_%s"%(what,sumIdx)
      else: sumString = "scaling_%s_%s"%(what,sumIdx)

      if sumString in sumTerms: sumTerms[sumString] += "%s,"%termList[termIdx]
      else: sumTerms[ sumString ] = "%s,"%termList[termIdx]

    #Add sizeable sums as RooAdditions
    for key, value in sumTerms.iteritems(): self.modelBuilder.factory_( "sum::%s(%s)"%(key,value[:-1]) )

    #Define string for total: 1 + sizeable sums
    if( isDecay )&( what != "tot" ): totalStr = "sum::scaling_partial_%s(1,"%what
    else: totalStr = "sum::scaling_%s(1,"%what
    for key in sumTerms: totalStr += "%s,"%key
    totalStr = totalStr[:-1]+")"
          
    #Add scaling function as RooAddition into model
    self.modelBuilder.factory_( totalStr )

  #Function to make BR scaling functions: partial width/total width
  def makeBRScalingFunction( self, what ): self.modelBuilder.factory_( 'expr::scaling_BR_%s("@0/@1", scaling_partial_%s, scaling_tot)'%(what,what) )
    
#################################################################################################################
# Define inherited classes: AllStageToEFT and StageXToEFT

# Combination of different stages
class STXSToSMEFTModel(STXStoSMEFTBaseModel):
  def __init__(self):
    STXStoSMEFTBaseModel.__init__(self)

  def setPhysicsOptions(self,physOptions):
    self.setPhysicsOptionsBase(physOptions)
  
  def doParametersOfInterest(self):
    if self.floatMass: print " --> [WARNING] Floating Higgs mass selected. STXStoSMEFT model assumes MH=125.0 GeV"
    self.doMH()
    self.SMH = SMHiggsBuilder(self.modelBuilder)
    
    #Read in parameters of interest from yaml file
    self.extractPOIs("/nfs/dust/cms/user/nigamova/CombineUpdate/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/data/eft/smeft_pois.yaml")
    #self.extractPOIs("/vols/cms/jl2117/postdoc/hcomb/EFT/Oct21/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/data/eft/smeft_pois_thesis.yaml")
    # Create list of pois and build RooRealVars
    POIs = ','.join(self.pois.keys())
    for poi in self.pois: self.modelBuilder.doVar("%s[%g,%g,%g]"%(poi,self.pois[poi]['val'],self.pois[poi]['min'],self.pois[poi]['max']))
    self.modelBuilder.doSet("POI",POIs)

    #set up model
    self.setup()

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def setup(self):
 
    # Extract scaling terms from json files: inclusive vs reco-level
    if not self.ignoreAccCorr:
      self.extractSTXSScalingTerms(filename="/vols/cms/jl2117/postdoc/hcomb/EFT/Oct21/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/data/eft/hgg_recolevel_pass2.json")    
      self.extractDecayScalingTerms(filename="/vols/cms/jl2117/postdoc/hcomb/EFT/Oct21/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/data/eft/decay_inclusive_thesis.json")    
    else:
      self.extractSTXSScalingTerms(filename="/nfs/dust/cms/user/nigamova/CombineUpdate/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/data/eft/VH.json")    
      self.extractDecayScalingTerms(filename="/nfs/dust/cms/user/nigamova/CombineUpdate/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/data/eft/decay_inclusive_thesis.json")    

    # Make total scaling function for decay side
    self.makeScalingFunction("tot", isDecay=True)
    
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def getHiggsSignalYieldScale(self,production,decay,energy):

    # Function to convert troublesome procs into viable one for HC combination
    production = convert_to_STXS(production,decay)

    name = "stxstosmeft_scaling_%s_%s_%s"%(production,decay,energy)

    if self.modelBuilder.out.function(name) == None:

      # Build scaling functions if they do not exist
      if self.modelBuilder.out.function("scaling_%s"%production) == None:
        self.makeScalingFunction(production)
      XSscal = "scaling_%s"%production

      if self.modelBuilder.out.function("scaling_BR_%s"%decay) == None:
        self.makeScalingFunction(decay, isDecay=True)
        self.makeBRScalingFunction(decay)
      BRscal = "scaling_BR_%s"%decay
  
      #Combine XS and BR scaling: incuding theory unc if option selected
      self.modelBuilder.factory_("prod::%s(%s)"%(name,",".join([XSscal,BRscal])))
      
    return name
 

#################################################################################################################
# Function to convert troublesome procs to the names in the json files
def convert_to_STXS( _production, _decay ):
  # Add string replace functions
  return _production


#################################################################################################################
# Instantiation of STXStoSMEFT model
STXStoSMEFT = STXSToSMEFTModel()

