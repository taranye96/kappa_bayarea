#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 00:56:43 2021

@author: tnye
"""

###############################################################################
# Module with functions to call the following GMMs: Abrahamson et al. 2014,
# Boore 2014 and Boore 2020.
###############################################################################


def ask14(IM,M,Rrup,vs30=760,ztor=7.13,rake=0.0,dip=90.0,width=10.0,z1pt0=0.05):
    """
    Computes PGA or PGV with Abrahamson et al. (2014) GMPE using OpenQuake engine.
    
        Inputs:
            imt(string): IM (PGA or PGV)
            M(float): Magnitude
            Rrup(float): Rrup in km
            vs30(float): Vs30 in m/s
            ztor(float): Depth to top of rupture. Default: 7.13 from Annemarie. ASSUMPTION IS RRUP > ZTOR
            rake(float): Rake of fault. Default is 0 degrees.
            dip(float): Dip of fault. Default is 90 degrees. 
            width(float): Width of fault. Default is 10 km.
            z1pt0(float): Soil depth to Vs = 1km/s, in m.  Default is 50.
            
        Return:
            lmean_ask14(float): ln Mean PGA in units of %g
            sd_ask14(float): Standard deviation of PGA
            
    """

    import numpy as np
    import pandas as pd
    from openquake.hazardlib.gsim.abrahamson_2014 import AbrahamsonEtAl2014
    from openquake.hazardlib import const, imt
    from openquake.hazardlib.gsim.base import RuptureContext
    from openquake.hazardlib.gsim.base import DistancesContext
    from openquake.hazardlib.gsim.base import SitesContext
    
    # Initiate model
    ASK14  = AbrahamsonEtAl2014()
    
    # Define intensity measure
    if IM == 'PGA':
        IMT = imt.PGA()
    elif IM == 'PGV':
        IMT = imt.PGV()
    
    # Initiate the rupture, distances, and sites objects:
    rctx = RuptureContext()
    dctx = DistancesContext()
    sctx = SitesContext()
    sctx_rock = SitesContext()
    
    # Fill the rupture context...assuming rake is 0, dip is 90
    rctx.rake = rake
    rctx.dip = dip
    rctx.ztor = ztor
    rctx.width = width 
    
    # Scenario I: If M and Rrup are both single values:
    #   Then set the magnitude as a float, and the rrup/distance as an array
    #   of one value
    if isinstance(M,float) & isinstance(Rrup,float):
        rctx.mag = M
        dctx.rrup = np.logspace(np.log10(Rrup),np.log10(Rrup),1)

        # Then compute everything else...
        #   Assuming average ztor, get rjb:
        dctx.rjb = np.sqrt(dctx.rrup**2 - rctx.ztor**2)
        dctx.rhypo = dctx.rrup
        dctx.rx = dctx.rjb
        dctx.ry0 = dctx.rx
        
        #   Set site parameters
        sctx.vs30 = np.ones_like(dctx.rrup) * vs30
        sctx.vs30measured = np.full_like(dctx.rrup, False, dtype='bool')
        sctx.z1pt0 = np.ones_like(dctx.rrup) * z1pt0
        
        #  Compute prediction
        lnmean_ask14, sd_ask14 = ASK14.get_mean_and_stddevs(
            sctx, rctx, dctx, IMT, [const.StdDev.TOTAL])
            
        return(lnmean_ask14, sd_ask14)


def boore2014(M,Rrup,vs30=760,ztor=7.13,rake=0.0):
    """
    Computes PGA with Boore et al. (2014) GMPE using OpenQuake engine.
    
        Inputs:
            M(float): Magnitude
            Rrup(float): Rrup in km
            vs30(float): Vs30 in m/s
            ztor(float): Depth to top of rupture. Default: 7.13 from Annemarie. ASSUMPTION IS RRUP > ZTOR
            rake(float): Rake of fault. Default is 0 degrees.
            dip(float): Dip of fault. Default is 90 degrees. 
            width(float): Width of fault. Default is 10 km.
            
        Return:
            lmean_boore14(float): Mean PGA
            sd_boore14(float): Standard deviation of PGA
            
    """

    import numpy as np
    import pandas as pd
    from openquake.hazardlib.gsim.boore_2014 import BooreEtAl2014
    from openquake.hazardlib import const, imt
    from openquake.hazardlib.gsim.base import RuptureContext
    from openquake.hazardlib.gsim.base import DistancesContext
    from openquake.hazardlib.gsim.base import SitesContext
    
    # Initiate model
    boore  = BooreEtAl2014()
    
    # Define intensity measure
    IMT = imt.PGA()
    
    # Initiate the rupture, distances, and sites objects:
    rctx = RuptureContext()
    dctx = DistancesContext()
    sctx = SitesContext()
    sctx_rock = SitesContext()
    
    # Fill the rupture context...assuming rake is 0, dip is 90
    rctx.rake = rake 
    rctx.ztor = ztor
    
    # Scenario I: If M and Rrup are both single values:
    #   Then set the magnitude as a float, and the rrup/distance as an array
    #   of one value
    if isinstance(M,float) & isinstance(Rrup,float):
        rctx.mag = M
        dctx.rrup = np.logspace(np.log10(Rrup),np.log10(Rrup),1)

        # Then compute everything else...
        #   Assuming average ztor, get rjb:
        dctx.rjb = np.sqrt(dctx.rrup**2 - rctx.ztor**2)
        
        #   Set site parameters
        sctx.vs30 = np.ones_like(dctx.rrup) * vs30
        sctx.vs30measured = np.full_like(dctx.rrup, False, dtype='bool')
        
        #  Compute prediction
        lnmean_boore14, sd_boore14 = boore.get_mean_and_stddevs(
            sctx, rctx, dctx, IMT, [const.StdDev.TOTAL])
            
        return(lnmean_boore14, sd_boore14)
    

def oq_boore2020(M,Rrup,predictive_parameter='pga',vs30=760,ztor=7.13,rake=0.0,dip=90.0,width=10.0,z1pt0 = 0.05):
    '''
    Compute the predicted ground motions with Abrahamson, Silva, and Kamai 2014 model
        from OpenQuake engine.  Assuming all events are a point source. (Sent to me by Avigyan Chaterjee)
    Input:
        M:                      Float or array with magnitudes to compute
        Rrup:                   Float or array with rrups - if it's an array, it should be np.logspace(log10(start),log10(stop),num)
        predictive_parameter:   Predictive parameter to compute: 'pga','pgv', or float with SA period (i.e., 1.0).  Default: 'pga'
        vs30:                   Value or array with Vs30 to use.  Default: 760. 
        ztor:                   Depth to top of rupture. Default: 7.13 from Annemarie. ASSUMPTION IS RRUP > ZTOR!!!!
        rake:                   Rake.  Default: 0.0 degrees.
        dip:                    Dip.  Default: 90.0 degrees.
        width:                  Fault width.  Default: 10.0
        z1pt0:                  Soil depth to Vs = 1.0km/s, in km.  Default: 0.05.
    Output:
        lmean_b14:            Mean ground motion. If M and Rrup floats returns float, if M float and Rrup array returns array like Rrup,
                                    if M array and Rrup float returns array like M, if M and rrup arrays returns array like len(M) x len(Rrup)
        sd_boore14:               Standard deviation. If M and Rrup floats returns float, if M float and Rrup array returns array like Rrup,
                                    if M array and Rrup float returns array like M, if M and rrup arrays returns array like len(M) x len(Rrup)
        
    '''
    from boore_2020 import BooreEtAl2020
    from openquake.hazardlib import imt, const
    from openquake.hazardlib.gsim.base import RuptureContext
    from openquake.hazardlib.gsim.base import DistancesContext
    from openquake.hazardlib.gsim.base import SitesContext
    import numpy as np
    
    # Initiate which model:
    B20  = BooreEtAl2020()

    # Predictive parameter:
    if predictive_parameter=='pga':
        IMT = imt.PGA()
    elif predictive_parameter=='pgv':
        IMT = imt.PGV()
    else:
        IMT = imt.SA(predictive_parameter)
        
    # Initiate the rupture, distances, and sites objects:
    rctx = RuptureContext()
    dctx = DistancesContext()
    sctx = SitesContext()
    sctx_rock = SitesContext()

    # Fill the rupture context...assuming rake is 0, dip is 90,
    rctx.rake = rake
    rctx.dip = dip
    rctx.ztor = ztor
    rctx.width = width   
    
    # Scenario I: If M and Rrup are both single values:
    #   Then set the magnitude as a float, and the rrup/distance as an array
    #   of one value
    if isinstance(M,float) & isinstance(Rrup,float):
        rctx.mag = M
        dctx.rrup = np.logspace(np.log10(Rrup),np.log10(Rrup),1)

        # Then compute everything else...
	#    Assuming average ztor, get rjb:
        dctx.rjb = np.sqrt(dctx.rrup**2 - rctx.ztor**2)
        dctx.rhypo = dctx.rrup
        dctx.rx = dctx.rjb
        dctx.ry0 = dctx.rx
        
        #   Set site parameters
        sctx.vs30 = np.ones_like(dctx.rrup) * vs30
        sctx.vs30measured = np.full_like(dctx.rrup, False, dtype='bool')
        sctx.z1pt0 = np.ones_like(dctx.rrup) * z1pt0
        
        # Compute prediction
        lmean_b20, sd_boore20 = B20.get_mean_and_stddevs(
            sctx, rctx, dctx, IMT, [const.StdDev.TOTAL])
            
        return lmean_b20, sd_boore20

    
    # Scenario II: If M is a single value and Rrup is an array:
    if isinstance(M,float) & isinstance(Rrup,np.ndarray):
        # Set them as intended...Rrup should be in logspace
        rctx.mag = M
        dctx.rrup = Rrup
        
	# Assuming average ztor, get rjb:
        dctx.rjb = np.sqrt(dctx.rrup**2 - rctx.ztor**2)
        dctx.rhypo = dctx.rrup
        dctx.rx = dctx.rjb
        dctx.ry0 = dctx.rx
        
        sctx.vs30 = np.ones_like(dctx.rrup) * vs30
        sctx.vs30measured = np.full_like(dctx.rrup, False, dtype='bool')
        sctx.z1pt0 = np.ones_like(dctx.rrup) * z1pt0
        
        lmean_b20, sd_boore20 = B20.get_mean_and_stddevs(
            sctx, rctx, dctx, IMT, [const.StdDev.TOTAL])
            
        return lmean_b20, sd_boore20
    
    
    # Scenario III: If M is an array and Rrup is a single value:
    if isinstance(M,np.ndarray) & isinstance(Rrup,float):
        # Set dctx to be a single value array, like in scenario I:
        dctx.rrup = np.logspace(np.log10(Rrup),np.log10(Rrup),1)
        
        # The rest of dctx depends only on rrup, as wella s site, so populate those:
	# Assuming average ztor, get rjb:
        dctx.rjb = np.sqrt(dctx.rrup**2 - rctx.ztor**2)
        dctx.rhypo = dctx.rrup
        dctx.rx = dctx.rjb
        dctx.ry0 = dctx.rx
        
        # Site: 
        sctx.vs30 = np.ones_like(dctx.rrup) * vs30
        sctx.vs30measured = np.full_like(dctx.rrup, False, dtype='bool')
        sctx.z1pt0 = np.ones_like(dctx.rrup) * z1pt0
        
        # But rctx depends on M and can only take a float, so will have to run many times.
        # Initiate mean and std lists:
        lmean_b20 = np.zeros_like(M)
        sd_boore20 = np.zeros_like(M)
        
        # Then loop over M's for rctx:
        for iMag in range(len(M)):
            # Set mag:
            rctx.mag = M[iMag]
            
            # Set 
            i_lmean_b20, i_sd_boore20 = B20.get_mean_and_stddevs(
                sctx, rctx, dctx, IMT, [const.StdDev.TOTAL])
                
            lmean_b20[iMag] = i_lmean_b20
            sd_boore20[iMag] = i_sd_boore20[0]
            
        return lmean_b20,sd_boore20
    
     
    # If both M and Rrup are arrays: 
    if isinstance(M,np.ndarray) & isinstance(Rrup,np.ndarray):
        # Set dctx to be its array as intended:
        dctx.rrup = Rrup
        
        # The rest of dctx depends only on rrup, as wella s site, so populate those:
	# Assuming average ztor, get rjb:	
        #dctx.rjb = np.log10(np.sqrt((10**dctx.rrup)**2 - rctx.ztor**2))
        dctx.rjb = dctx.rrup
        dctx.rhypo = dctx.rrup
        dctx.rx = dctx.rjb
        dctx.ry0 = dctx.rx
        
        # Site: 
        sctx.vs30 = np.ones_like(dctx.rrup) * vs30
        sctx.vs30measured = np.full_like(dctx.rrup, False, dtype='bool')
        sctx.z1pt0 = np.ones_like(dctx.rrup) * z1pt0
        
        # But rctx depends on M and can only take a float, so will have to run many times.
        # Initiate mean and std lists:
        lmean_b20 = np.zeros((len(M),len(Rrup)))
        sd_boore20 = np.zeros((len(M),len(Rrup)))
        
        # Then loop over M's for rctx:
        for iMag in range(len(M)):
            # Set mag:
            rctx.mag = M[iMag]
            
            # Set 
            i_lmean_b20, i_sd_boore20 = B20.get_mean_and_stddevs(
                sctx, rctx, dctx, IMT, [const.StdDev.TOTAL])
                
            lmean_b20[iMag,:] = i_lmean_b20
            sd_boore20[iMag,:] = i_sd_boore20[0]
        
        return lmean_b20,sd_boore20
    
