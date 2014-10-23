#
# Python package implementing the parametric model described in ApJ paper
#
# $Source: /nfs/slac/g/glast/ground/cvs/users/oxon/PyCRFlux/src/ext/ParamModel.py,v $
# $Author: oxon $ $Date: 2009/05/19 17:04:10 $ $Revision: 1.1 $
#

#
# Provides convencience classes that implements calculations of cross sections
# based on the parametric model published in ApJ. The classes are wrappers around
# the c functions implemented in cparamlib.
#

import sys
import math
import cparamlib

# --------------------------------------------------------------------------------
#  Helper function definitions
# --------------------------------------------------------------------------------

#
# Function: fcmp
#
# Compares two floating point numbers
#
def fcmp(x, y):
    d = math.fabs(x - y)
    if (d < 1.0e-15):
        return 1
    else:
        return 0

# --------------------------------------------------------------------------------
#  Class definitions
# --------------------------------------------------------------------------------

#
# Class ParamModel
#
# Description:
# The class ParamModel is a python representation of the parameteric model implemented
# in cparamlib. It is just a convenience class which wraps the implemented functions.
# On instantiation, the parameters are caluclated for the given Tp and particle.
#
# Inclusive cross section calculated with calls to sigma_incl_<process> which takes the
# secondary particle energy E and proton kinetic energy as arguments. If Tp is different
# than from instantiation or previous call, then parameters are recalculated and Tp stored.
#
class ParamModel:

    #
    # Constructor
    #
    def __init__(self, Tp, particle = None):
        if (particle == None):
            # default to gamma-rays
            self.particle = cparamlib.ID_GAMMA
        else:
            self.particle = particle

        self.Tp = Tp

        # create struct for parameter storage
        self.params = cparamlib.PARAMSET()

        # initialize parameters for the given Tp
        self.param_all(self.Tp)

    #
    # Method param_all
    #
    # Calculates the parameters for all four components at the given Tp
    # This is a wrapper around cparamlib.<particle>_param(Tp, params)
    #
    def param_all(self, Tp):
        if (self.particle == cparamlib.ID_GAMMA):
            cparamlib.gamma_param(Tp, self.params)
        elif (self.particle == cparamlib.ID_ELECTRON):
            cparamlib.elec_param(Tp, self.params)
        elif (self.particle == cparamlib.ID_POSITRON):
            cparamlib.posi_param(Tp, self.params)
        elif (self.particle == cparamlib.ID_NUE):
            cparamlib.nue_param(Tp, self.params)
        elif (self.particle == cparamlib.ID_NUMU):
            cparamlib.numu_param(Tp, self.params)
        elif (self.particle == cparamlib.ID_ANTINUE):
            cparamlib.antinue_param(Tp, self.params)
        elif (self.particle == cparamlib.ID_ANTINUMU):
            cparamlib.antinumu_param(Tp, self.params)

    #
    # Method param_nd
    #
    # Calculates the parameters for the non-diffraction parameter set at the given Tp
    # This is a wrapper around cparamlib.<particle>_param_nd(Tp, params)
    #
    def param_nd(self, Tp):
        if (self.particle == cparamlib.ID_GAMMA):
            cparamlib.gamma_param_nd(Tp, self.params)
        elif (self.particle == cparamlib.ID_ELECTRON):
            cparamlib.elec_param_nd(Tp, self.params)
        elif (self.particle == cparamlib.ID_POSITRON):
            cparamlib.posi_param_nd(Tp, self.params)
        elif (self.particle == cparamlib.ID_NUE):
            cparamlib.nue_param_nd(Tp, self.params)
        elif (self.particle == cparamlib.ID_NUMU):
            cparamlib.numu_param_nd(Tp, self.params)
        elif (self.particle == cparamlib.ID_ANTINUE):
            cparamlib.antinue_param_nd(Tp, self.params)
        elif (self.particle == cparamlib.ID_ANTINUMU):
            cparamlib.antinumu_param_nd(Tp, self.params)

    #
    # Method param_diff
    #
    # Calculates the parameters for the diffraction parameter set at the given Tp
    # This is a wrapper around cparamlib.<particle>_param_diff(Tp, params)
    #
    def param_diff(self, Tp):
        if (self.particle == cparamlib.ID_GAMMA):
            cparamlib.gamma_param_diff(Tp, self.params)
        elif (self.particle == cparamlib.ID_ELECTRON):
            cparamlib.elec_param_diff(Tp, self.params)
        elif (self.particle == cparamlib.ID_POSITRON):
            cparamlib.posi_param_diff(Tp, self.params)
        elif (self.particle == cparamlib.ID_NUE):
            cparamlib.nue_param_diff(Tp, self.params)
        elif (self.particle == cparamlib.ID_NUMU):
            cparamlib.numu_param_diff(Tp, self.params)
        elif (self.particle == cparamlib.ID_ANTINUE):
            cparamlib.antinue_param_diff(Tp, self.params)
        elif (self.particle == cparamlib.ID_ANTINUMU):
            cparamlib.antinumu_param_diff(Tp, self.params)

    #
    # Method param_delta
    #
    # Calculates the parameters for the delta(1232) parameter set at the given Tp
    # This is a wrapper around cparamlib.<particle>_param_delta(Tp, params)
    #
    def param_delta(self, Tp):
        if (self.particle == cparamlib.ID_GAMMA):
            cparamlib.gamma_param_delta(Tp, self.params)
        elif (self.particle == cparamlib.ID_ELECTRON):
            cparamlib.elec_param_delta(Tp, self.params)
        elif (self.particle == cparamlib.ID_POSITRON):
            cparamlib.posi_param_delta(Tp, self.params)
        elif (self.particle == cparamlib.ID_NUE):
            cparamlib.nue_param_delta(Tp, self.params)
        elif (self.particle == cparamlib.ID_NUMU):
            cparamlib.numu_param_delta(Tp, self.params)
        elif (self.particle == cparamlib.ID_ANTINUE):
            cparamlib.antinue_param_delta(Tp, self.params)
        elif (self.particle == cparamlib.ID_ANTINUMU):
            cparamlib.antinumu_param_delta(Tp, self.params)
    
    #
    # Method param_res
    #
    # Calculates the parameters for the res(1600) parameter set at the given Tp
    # This is a wrapper around cparamlib.<particle>_param_res(Tp, params)
    #
    def param_res(self, Tp):
        if (self.particle == cparamlib.ID_GAMMA):
            cparamlib.gamma_param_res(Tp, self.params)
        elif (self.particle == cparamlib.ID_ELECTRON):
            cparamlib.elec_param_res(Tp, self.params)
        elif (self.particle == cparamlib.ID_POSITRON):
            cparamlib.posi_param_res(Tp, self.params)
        elif (self.particle == cparamlib.ID_NUE):
            cparamlib.nue_param_res(Tp, self.params)
        elif (self.particle == cparamlib.ID_NUMU):
            cparamlib.numu_param_res(Tp, self.params)
        elif (self.particle == cparamlib.ID_ANTINUE):
            cparamlib.antinue_param_res(Tp, self.params)
        elif (self.particle == cparamlib.ID_ANTINUMU):
            cparamlib.antinumu_param_res(Tp, self.params)

    #
    # Method sigma_incl_nd
    #
    # Calculates the inclusive non-diffraction cross section at the given E and Tp
    # This is a wrapper around cparamlib.sigma_incl_nd(particle, E, Tp, params)
    #
    def sigma_incl_nd(self, E, Tp):
        if (Tp != self.Tp):
            self.Tp = Tp
            self.param_nd(self.Tp)

        sigma = cparamlib.sigma_incl_nd(self.particle, E, self.Tp, self.params)
            
        return sigma
            
    #
    # Method sigma_incl_diff
    #
    # Calculates the inclusive diffraction cross section at the given E and Tp
    # This is a wrapper around cparamlib.sigma_incl_diff(particle, E, Tp, params)
    #
    def sigma_incl_diff(self, E, Tp):
        if (Tp != self.Tp):
            self.Tp = Tp
            self.param_diff(self.Tp)

        sigma = cparamlib.sigma_incl_diff(self.particle, E, self.Tp, self.params)
            
        return sigma

    #
    # Method sigma_incl_delta
    #
    # Calculates the inclusive delta(1232) cross section at the given E and Tp
    # This is a wrapper around cparamlib.sigma_incl_delta(particle, E, Tp, params)
    #
    def sigma_incl_delta(self, E, Tp):
        if (Tp != self.Tp):
            self.Tp = Tp
            self.param_delta(self.Tp)

        sigma = cparamlib.sigma_incl_delta(self.particle, E, self.Tp, self.params)
            
        return sigma

    #
    # Method sigma_incl_res
    #
    # Calculates the inclusive res(1600) cross section at given E and Tp
    # This is a wrapper around cparamlib.sigma_incl_res(particle, E, Tp, params)
    #
    def sigma_incl_res(self, E, Tp):
        if (Tp != self.Tp):
            self.Tp = Tp
            self.param_res(self.Tp)

        sigma = cparamlib.sigma_incl_res(self.particle, E, self.Tp, self.params)
            
        return sigma

    #
    # Method sigma_incl_tot
    #
    # Calculates the total inclusive cross section, i.e. from all four components,
    # at the given E and Tp
    # This is a wrapper around cparamlib.sigma_incl_tot(particle, E, Tp, params)
    #
    def sigma_incl_tot(self, E, Tp):
        if (Tp != self.Tp):
            self.Tp = Tp
            self.param_all(self.Tp)
            
        sigma = cparamlib.sigma_incl_tot(self.particle, E, self.Tp, self.params)
            
        return sigma

#
# Class AngularParamModel
#
# Description:
# The class AngularParamModel is a python representation of the parameteric model for
# the angular distribution implemented in cparamlib. It is jus ta convenience class
# which wraps the implemented functions. On instantiation, the parameters are caluclated
# for the given Tp, E and particle.
#
# Calculate the differential cross section d^2sigma/dlogEdpT for gamma rays at given
# pT, E and Tp. Consecutive calls checks if E and Tp changed and recalculates parameters
# only as neccesary.
#
class AngularParamModel:

    #
    # Constructor
    #
    def __init__(self, E, Tp, particle = None):
        if (particle == None):
            # default to gamma-rays
            self.particle = cparamlib.ID_GAMMA
        else:
            self.particle = particle

        self.Tp = Tp
        self.E = E
        self.chTp_all = 1
        self.chTp_nr = 1
        self.chTp_delta = 1
        self.chTp_res = 1

        # create struct for parameter storage
        self.pt_params = cparamlib.PARAMSET_PT()
            
        # initialize parameters for the given Tp and E
        self.param_pt_all(self.E, self.Tp)
        
    #
    # Method param_pt_all
    #
    # Calculates the pT parameters for all three components at the given Tp and E
    # This is a wrapper around cparamlib.<particle>_pt_param(E, Tp, params, flag)
    #
    def param_pt_all(self, E, Tp):
        if (self.particle == cparamlib.ID_GAMMA):
            cparamlib.gamma_pt_param(E, Tp, self.pt_params, self.chTp_all)
                
    #
    # Method param_pt_nr
    #
    # Creates the list containing the non-resonance pT parameter set at the given Tp
    # This is a wrapper around cparamlib.<particle>_pt_param_nr(E, Tp, params, flag)
    #
    def param_pt_nr(self, E, Tp):
        if (self.particle == cparamlib.ID_GAMMA):
            cparamlib.gamma_pt_param_nr(E, Tp, self.pt_params, self.chTp_nr)

    #
    # Method param_pt_delta
    #
    # Creates the list containing the Delta(1232) pT parameter set at the given Tp
    # This is a wrapper around cparamlib.<particle>_pt_param_delta(E, Tp, params, flag)
    #
    def param_pt_delta(self, E, Tp):
        if (self.particle == cparamlib.ID_GAMMA):
            cparamlib.gamma_pt_param_delta(E, Tp, self.pt_params, self.chTp_delta)

    #
    # Method param_pt_res
    #
    # Creates the list containing the res(1232) pT parameter set at the given Tp
    # This is a wrapper around cparamlib.<particle>_pt_param_res(E, Tp, params, flag)
    #
    def param_pt_res(self, E, Tp):
        if (self.particle == cparamlib.ID_GAMMA):
            cparamlib.gamma_pt_param_res(E, Tp, self.pt_params, self.chTp_res)

    #
    # Method sigma_pt_nr
    #
    # Calculates the differential cross section d^2sigma/dlogEdpT for the non-resonance
    # component at given pT, E and Tp
    # This is a wrapper around cparamlib.sigma_pt_nr(particle, E, Tp, params)
    #
    def sigma_pt_nr(self, pT, E, Tp):
        # check whether Tp has changed from the previous call
        # if it has, set the chTp flag to 1 and update the stored Tp value
        if fcmp(Tp, self.Tp):
            self.chTp_nr = 0
        else:
            print "Tp changed from", self.Tp, "to", Tp
            self.Tp = Tp
            self.chTp_nr = 1

        # now check if E has changed since the last call or if chTp is set
        if ((fcmp(E, self.E) == 0) or (self.chTp_nr == 1)):
            self.E = E
            self.param_pt_nr(self.E, self.Tp)

        sigma = cparamlib.sigma_pt_nr(self.particle, pT, self.E, self.Tp, self.pt_params)
            
        return sigma

    #
    # Method sigma_pt_delta
    #
    # Calculates the differential cross section d^2sigma/dlogEdpT for Delta(1232)
    # at given pT, E and Tp
    # This is a wrapper around cparamlib.sigma_pt_delta(particle, E, Tp, params)
    #
    def sigma_pt_delta(self, pT, E, Tp):
        # check whether Tp has changed from the previous call
        # if it has, set the chTp flag to 1 and update the stored Tp value
        if fcmp(Tp, self.Tp):
            self.chTp_delta = 0
        else:
            print "Tp changed from", self.Tp, "to", Tp
            self.Tp = Tp
            self.chTp_delta = 1
        
        # now check if E has changed since the last call or if chTp is set
        if ((fcmp(E, self.E) == 0) or (self.chTp_delta == 1)):
            self.E = E
            self.param_pt_delta(self.E, self.Tp)

        sigma = cparamlib.sigma_pt_delta(self.particle, pT, self.E, self.Tp, self.pt_params)
            
        return sigma

    #
    # Method sigma_pt_res
    #
    # Calculates the differential cross section d^2sigma/dlogEdpT for res(1600)
    # at given pT, E and Tp
    # This is a wrapper around cparamlib.sigma_pt_res(particle, E, Tp, params)
    #
    def sigma_pt_res(self, pT, E, Tp):
        # first check whether Tp has changed from the previous call
        # if it has, set the chTp flag to 1 and update the stored Tp value
        if fcmp(Tp, self.Tp):
            self.chTp_res = 0
        else:
            print "Tp changed from", self.Tp, "to", Tp
            self.Tp = Tp
            self.chTp_res = 1
        
        # now check if E has changed since the last call or if chTp is set
        if ((fcmp(E, self.E) == 0) or (self.chTp_res == 1)):
            self.E = E
            self.param_pt_res(self.E, self.Tp)

        sigma = cparamlib.sigma_pt_res(self.particle, pT, self.E, self.Tp, self.pt_params)
            
        return sigma

    #
    # Method sigma_pt_tot
    #
    # Calculates the differential cross section d^2sigma/dlogEdpT from all three
    # components at the given pT, E and Tp
    # This is a wrapper around cparamlib.sigma_pt_tot(particle, E, Tp, params)
    #
    def sigma_pt_tot(self, pT, E, Tp):
        # first check whether Tp has changed from the previous call
        # if it has, set the chTp flag to 1 and update the stored Tp value
        if fcmp(Tp, self.Tp):
            self.chTp_all = 0
        else:
            print "Tp changed from", self.Tp, "to", Tp
            self.Tp = Tp
            self.chTp_all = 1
        
        # now check if E has changed since the last call or if chTp is set
        if ((fcmp(E, self.E) == 0) or (self.chTp_all == 1)):
            self.E = E
            self.param_pt_all(self.E, self.Tp)

        sigma = cparamlib.sigma_pt_tot(self.particle, pT, self.E, self.Tp, self.pt_params)
            
        return sigma
    
        
if (__name__ == "__main__"):
    # sys.exit(0)
    import ParamModel as mymodel   
    ID_GAMMA=0    
    Tp=[1.0,2.0,3.0,4.0,5.0]
    E=[0.5,1.0,1.5,2.0,2.5,3.0]
    mymodel.ParamModel(Tp,ID_GAMMA)
    print mymodel.sigma_incl_nd(Tp,E)    

