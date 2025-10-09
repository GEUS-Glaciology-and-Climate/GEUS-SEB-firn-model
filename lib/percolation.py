import numpy as np
from numba import jit
import matplotlib.pyplot as plt


@jit(nopython=True)
def calc_darcy_fluxes(pslwc, psnowc, psnic, pdgrain, prhofirn, zdtime):
    # calc_darcy_fluxes: Calculates the amount of water (mm weq) that each layer
    # can transmit to the next one according to a Darcy flow.
    #
    #   syntax:
    #   [darcy_fluxes] = calc_darcy_fluxes (pslwc, psnowc, psnic, pdgrain, prhofirn, c)
    #
    #   input:
    #         psnowc, psnic, pslwc - vectors of respectively snow, ice and
    #         liquid water part for each subsurface layer (mm weq). Their sum
    #         for each layer should always be equal to the layer fixed water
    #         eq. thickness.
    #
    #         pdgrain - Vector of layer grain size. see graingrowth def.
    #
    #         prhofirn - vector of firn (snow) density (kg/m**3)
    #
    #         c - Structure containing all the physical, site-depant or user
    #         defined variables.
    #
    #   output:
    #         darcy_fluxes - vector containing the amount of water (mm weq) that each layer
    #         can transmit to the next one according to a Darcy flow.
    #
    #   This script was originally developped by Peter Langen (pla@dmi.dk) and
    #   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to python by
    #   Baptiste Vandecrux (bav@geus.dk).
    #
    #   Aug. 2016
    # =========================================================================

    darcy_fluxes = np.zeros_like(pslwc)

    for jk in range(len(pslwc) - 1):  # 1:c.jpgrnd-1
        # Special treatment of top layer in case of no snow
        if (jk == 0) & (psnowc[jk] < 1e-12):
            # In case of no snow, the calculation below of K crashes. Instead we say
            # * if layer below may receive water, it gets all there is in layer 1
            # * if layer may not receive water, it all runs off.
            #
            # As long as we set darcy_fluxes[1] = pslwc[1] then this will
            # be taken care of in the perc_runoff code.
            darcy_fluxes[jk] = pslwc[jk]
        else:
            # Other layers do (should!) not come into a situation where there is only ice and water
            qlim = qlimF(
                pslwc[jk], pslwc[jk + 1],
                psnowc[jk], psnowc[jk + 1],
                psnic[jk], psnic[jk + 1],
                prhofirn[jk], prhofirn[jk + 1],
                pdgrain[jk], pdgrain[jk + 1],
            )
            Theta1 = ThetaF(pslwc[jk], psnowc[jk], prhofirn[jk])
            Theta2 = ThetaF(pslwc[jk + 1], psnowc[jk + 1], prhofirn[jk + 1])
            h1 = hHirF(Theta1, pdgrain[jk])
            h2 = hHirF(Theta2, pdgrain[jk + 1])

            delz1 = 999.8395 * (psnowc[jk] / prhofirn[jk] + psnic[jk] / 900)  # layer 1
            delz2 = 999.8395 * (
                psnowc[jk + 1] / prhofirn[jk + 1] + psnic[jk + 1] / 900
            )  # layer 2
            delz = (delz1 + delz2) / 2  # Midpoint-to-midpoint distance

            dhdz = (h2 - h1) / delz

            k1 = kF(Theta1, pdgrain[jk], prhofirn[jk], psnic[jk], psnowc[jk])
            k2 = kF(
                Theta2, pdgrain[jk + 1], prhofirn[jk + 1], psnic[jk + 1], psnowc[jk + 1]
            )
            k = (k1 + k2) / 2  # Average k between the two layers

            # Initial flux according to Hirashima eqn [1]
            q0 = max(k * (dhdz + 1), 0)

            # Total time-step flux according to Hirashima eqn (23)
            # and make sure it doesn't exceed available water in
            # the upper of the two layers:
            if qlim > 0:
                darcy_fluxes[jk] = min(
                    pslwc[jk], qlim * (1 - np.exp(-q0 / qlim * zdtime))
                )
            else:
                darcy_fluxes[jk] = 0
    return darcy_fluxes


@jit(nopython=True)
def CLliqF(rhosin):
    """Compute irreducible liquid water content per pore-space volume.

    Calculates the irreducible liquid water content following the
    parameterization of Coléou and Lesaffre (1998). The input snow density
    is capped at 830 kg/m³ to avoid nonphysical values near pore close-off.

    Args:
        rhosin (float): Snow density in kg/m³.

    Returns:
        float: Irreducible liquid water content per pore-space volume (fraction).
    """
    rhos = min(rhosin, 830)
    P = 1 - rhos / 900
    Wm = 0.057 * P / max(1 - P, 1e-12) + 0.017
    return Wm / (P / max(1 - P, 1e-12)) * 900 / 999.8395 / max(1 - Wm, 1e-12)


def hetero_percol(prhofirn, psnowc, psnic, pslwc, pdgrain, c):
    # hetero_percol
    #   syntax:
    # [prhofirn, psnowc, psnic, pslwc, pdgrain] = hetero_percol (prhofirn, psnowc, psnic, pslwc, pdgrain, c)
    #   input:
    #         prhofirn - vector of firn (snow) density (kg/m**3)
    #         psnowc, psnic, pslwc - vectors of respectively snow, ice and
    #         liquid water part for each subsurface layer (mm weq). Their sum
    #         for each layer should always be equal to the layer fixed water
    #         eq. thickness.
    #         ptsoil - vector of subsurface temperature (K)
    #         pdgrain - Vector of layer grain size. see graingrowth def.
    #         pTdeep - lower boundary temperature (K). Defined by the constant
    #         T_ice_AWS so far.
    #         zrogl - total amount of run off (mm weq)
    #         c - Structure containing all the physical, site-depant or user
    #         defined variables.
    #   output:
    #          updated [prhofirn, psnowc, psnic, pslwc, ptsoil, pdgrain, zrogl]
    #   This script was originally developped by Baptiste Vandecrux (bav@geus.dk)
    # =========================================================================
    avail_water = zeros(size(prhofirn))

    # in the formulation of Marchenko et al. 2017, the liquid water is
    # distributed according to a probability funcion accross the layers before
    # the standard percolation scheme takes over. In theory the heterogeneous
    # flow could start from any layer at depth, routing water even deeper
    # That is why we commented:
    # for jk = 1:c.jpgrnd-1
    # and work with all the water located in the 1st layer
    jk = 1

    # calculating the amount of water that can be held by capillary forces
    if c.calc_CLliq:
        liqmaxloc = CLliqF(prhofirn[jk])
    else:
        liqmaxloc = c.liqmax

    liqmaxM = liqmaxloc * 999.8395 / 900 * (900 / prhofirn[jk] - 1)
    potret = max(liqmaxM * psnowc[jk], 0)
    # and what is in excess
    liqexcess = pslwc[jk] - potret
    avail_water[jk] = max(liqexcess, 0)

    # if there is available water
    if avail_water[jk] > 1e-12:
        # we draw from a binomial distribution of mode hetero_percol_p to
        # know whether this layer leads to heterogeneous percolation
        if np.random.binomial(1, c.hetero_percol_p):
            percol_water = c.hetero_percol_frac * avail_water[jk]
            # note: dflux is positive when it leaves the layer

            # determine random depth between current layer and maximum depth range
            thickness_act = psnowc * (999.8395 / prhofirn) + psnic * (999.8395 / 900)
            depth_act = np.zeros_like(thickness_act)
            depth_act[1] = thickness_act[1] / 2
            for i in range(1, len(pslwc)):  # 2:size(thickness_act,1)
                depth_act[i] = np.sum(thickness_act[:i]) + thickness_act[i] / 2

            depth_current = depth_act[jk]
            depth_dest = c.hetero_percol_dist  # * rand(1,1) + depth_current
            # find index of layer at closest depth to the destination of
            # percolation
            ind_dest = np.argmin(np.abs(depth_act - depth_dest))

            # sav.depth_current = depth_current
            # sav.depth_dest = depth_dest

            # In the percolation scheme of Marchenko et al. (2017) the meltwater is
            # redistributed from the surface down to a specified destination depth
            # according to a probability def

            # in an uniform distribution def, the subsurface layers receive water
            # proportional to their thickness:

            # pslwc[jk] = pslwc[jk] - avail_water[jk]
            # frac_received = thickness_act(jk:ind_dest)./sum(thickness_act(jk:ind_dest))
            # pslwc(jk:ind_dest) = pslwc(jk:ind_dest) + avail_water[jk]*frac_received

            # an alternative is to go through the stratigraphy and stop when there is a
            # gradient in grain size or an ice layer
            for ii in range(jk, ind_dest):  # = jk:ind_dest-1
                # here we test all the layer through which the pipe travels
                if prhofirn(ii + 1) > 800:
                    # fprintf('ice - ')
                    break
                    # elseif psnic(ii+1)> c.ice_threshold
                    # plus relative to frozen mass
                    # if there is too much ie in the next layer
                    #                     fprintf('ice - ')
                    #                     break
                    #                 elseif ThetaF(pslwc(ii+1), psnowc(ii+1), prhofirn(ii+1), c) >= 1
                    # if there is a saturated layer piping can't go through
                    # it
                    #                     fprintf('sat - ')
                    #                     break
                    # PERMEABILITY? h?
                    #                 elseif pdgrain(ii)-pdgrain(ii+1) < c.crit_diff_grain
                    #                     fprintf('grain - ')
                    #                     break

            # if ind_dest ~= ii:
            #     print('Piping stopped at layer#i instead of#i\n',ii, ind_dest)
            ind_dest = ii
            jj = ind_dest
            while (jj > jk) & (percol_water > 1e-12):
                # Calculate water in destination layer, when this is at saturation (Theta = 1):
                plsat = psnowc[jj] * 999.8395 / 900 * (900 / prhofirn[jj] - 1)
                # Do not allow flux to be greater than plsat-pl in next layer.
                dflux = min(max(0, plsat - pslwc[jj]), percol_water)

                pslwc[jk] = pslwc[jk] - dflux
                pslwc[jj] = pslwc[jj] + dflux
                percol_water = percol_water - dflux
                jj = jj - 1
    return pslwc


@jit(nopython=True)
def hHirF(Theta, d):
    """Compute hydraulic suction head using Hirashima et al. (2010).

    Calculates the hydraulic suction head h (in meters) based on the
    volumetric water content and grain diameter, following equations from
    Hirashima et al. (2010).

    Args:
        Theta (float): Volumetric water content (dimensionless, 0–1).
        d (float): Mean grain diameter in meters.

    Returns:
        float: Hydraulic suction head in meters.
    """
    alpha = 7.3 * np.exp(1.9)
    n = nHirF(d)
    m = 1 - 1 / n
    Theta_nozero = max(Theta, 1e-12)
    return 1 / alpha * (Theta_nozero ** (-1 / m) - 1) ** (1 / n)



@jit(nopython=True)
def kF(Theta, d, rhos, pi, ps):
    """Compute vertical hydraulic conductivity accounting for ice lenses.

    Calculates the unsaturated vertical hydraulic conductivity of snow using:
    - Hirashima et al. (2010) for the unsaturated correction (eq. 10–11)
    - Calonne et al. (2012) for saturated permeability parameterization (Shimizu-form, Calonne-units)
    - Colbeck (1975, eq. 32) for vertical flow reduction due to horizontal ice lenses

    Args:
        Theta (float): Volumetric water content (dimensionless, 0–1).
        d (float): Mean grain diameter in millimeters.
        rhos (float): Snow density in kg/m³.
        pi (float): Thickness of ice in the layer (m).
        ps (float): Thickness of snow in the layer (m).

    Returns:
        float: Effective vertical hydraulic conductivity (m/s), reduced by ice layering.
    """
    # Saturated permeability (Calonne et al. 2012, based on Shimizu 1970)
    ks = 3 * (d / 2000) ** 2 * 9.82 / 1.79e-06 * np.exp(-0.013 * rhos)

    # Unsaturated correction factor (Hirashima et al. 2010, eq. 10)
    n = nHirF(d)
    m = 1 - 1 / n  # Hirashima eq. 9
    kr = Theta ** 0.5 * (1.0 - (1.0 - Theta ** (1.0 / m)) ** m) ** 2  # Hirashima eq. 10

    # Unsaturated hydraulic conductivity (Hirashima eq. 11)
    k11 = kr * ks

    # Ice layering correction (Colbeck 1975, eq. 32)
    Hs = ps / rhos  # total snow height [m]
    Hi = pi / 900   # total ice height [m]
    fsnow = Hs / (Hs + Hi)  # snow fraction of the layer
    whwice = 0.1  # scaling constant from Colbeck

    if k11 > 0:
        k22factor = fsnow + (1 - fsnow) * (1.0 + whwice) / (0.0 / k11 + whwice)  # Colbeck eq. 32
    else:
        k22factor = 1  # ensures kF = 0 when k11 = 0

    # Effective hydraulic conductivity (snow–ice layered system)
    return k11 / k22factor


@jit(nopython=True)
def nHirF(d):
    """Compute van Genuchten parameter n from grain size.

    Calculates the shape parameter n used in the van Genuchten model
    based on grain diameter, following Hirashima et al. (2010), eq. 17.

    Args:
        d (float): Mean grain diameter in millimeters.

    Returns:
        float: van Genuchten n parameter (dimensionless).
    """
    return 15.68 * np.exp(-0.46 * d) + 1



@jit(nopython=True)
def qlimF(pl1, pl2, ps1, ps2, pi1, pi2, rhos1, rhos2, d1, d2):
    """Compute water flux limit (qlim) between layers following Hirashima et al. (2010).

    Calculates the flux limit qlim at the interface between two snow layers using
    the approach described in Hirashima et al. (2010), eq. 20. The method tests
    for hydraulic equilibrium between the layers by comparing suction heads and
    iteratively solves for qlim using a bisection method if needed.

    Args:
        pl1 (float): Liquid water in top layer (mm w.e.).
        pl2 (float): Liquid water in bottom layer (mm w.e.).
        ps1 (float): Snow thickness in top layer (m).
        ps2 (float): Snow thickness in bottom layer (m).
        pi1 (float): Ice thickness in top layer (m).
        pi2 (float): Ice thickness in bottom layer (m).
        rhos1 (float): Snow density in top layer (kg/m³).
        rhos2 (float): Snow density in bottom layer (kg/m³).
        d1 (float): Grain size in top layer (mm).
        d2 (float): Grain size in bottom layer (mm).

    Returns:
        float: Water flux limit qlim (mm w.e.).
    """
    # Physical layer thicknesses (m), using density to convert ps to physical height
    delz1 = 999.8395 * (ps1 / rhos1 + pi1 / 900)  # Layer 1 height
    delz2 = 999.8395 * (ps2 / rhos2 + pi2 / 900)  # Layer 2 height
    delz = (delz1 + delz2) / 2  # Midpoint-to-midpoint distance

    # First, test if total transfer of pl1 achieves equilibrium
    qlim = pl1
    pl1test = pl1 - qlim
    pl2test = pl2 + qlim
    Theta1 = ThetaF(pl1test, ps1, rhos1)
    Theta2 = ThetaF(pl2test, ps2, rhos2)
    h1 = hHirF(Theta1, d1)
    h2 = hHirF(Theta2, d2)
    diff = h1 - (h2 + delz)  # Hirashima eq. 20: check if equilibrium is possible

    if diff < 0:
        Conv = 1  # No solution: all water can be transferred, return full pl1
        qlimOut = pl1
    else:
        Conv = 0  # Solution exists: iterate to find correct qlim

    if Conv == 0:
        # Start bisection with initial guess = half of water
        qlim = pl1 / 2
        qlimL = 0
        qlimR = pl1

        Nmax = 11  # Number of iterations (sufficient for ~1e-3 precision)
        for i in range(Nmax):
            plt.plot([i], [qlim], marker='o', color='k')  # Optional: track convergence

            pl1test = pl1 - qlim
            pl2test = pl2 + qlim
            Theta1 = ThetaF(pl1test, ps1, rhos1)
            Theta2 = ThetaF(pl2test, ps2, rhos2)
            h1 = hHirF(Theta1, d1)
            h2 = hHirF(Theta2, d2)
            diff = h1 - (h2 + delz)  # Hirashima eq. 20

            # Adjust interval depending on sign of difference
            if diff > 0:
                qlimR = qlim  # Moved too much
            else:
                qlimL = qlim  # Moved too little

            qlim = (qlimR + qlimL) / 2  # New midpoint

        qlimOut = qlim  # Final result after iteration

    return qlimOut


@jit(nopython=True)
def ThetaF(pl, ps, rhos):
    """Compute effective water saturation (Θ) following Coléou & Lesaffre (1998).

    Calculates the effective saturation (Θ) based on the layer’s liquid water content,
    snow mass, and density. The function caps Θ between 0 and 1, and uses the
    Coléou & Lesaffre (1998) parameterization for the local irreducible water content.

    Args:
        pl (float): Liquid water amount in the layer (mm w.e.).
        ps (float): Snow mass in the layer (kg/m², or equivalent thickness × density).
        rhos (float): Snow density (kg/m³).

    Returns:
        float: Effective water saturation Θ (dimensionless, 0–1).
    """
    # Irreducible liquid water content (Coléou & Lesaffre 1998)
    liqmaxloc = CLliqF(rhos)

    if ps > 1e-12:
        # Compute effective saturation Θ (bounded between 0 and 1)
        ThetaF = min(
            1,
            max(
                0,
                (
                    pl / ps * rhos * 900 / 999.8395 / max(1e-12, 900 - rhos)
                    - liqmaxloc
                )
                / (1 - liqmaxloc),
            ),
        )
    else:
        # No snow: set Θ = 1 to avoid division by zero
        ThetaF = 1

    return ThetaF



@jit(nopython=True)
def perc_runoff_new(prhofirn, psnowc, psnic, pslwc, pdgrain, zdtime):
    """Compute liquid water percolation and runoff through firn layers.

    Calculates percolation and runoff of meltwater through a vertical snow/firn column
    using either a simple bucket scheme or a Darcy-like percolation scheme. Each layer's
    snow, ice, and water masses are updated based on vertical water movement. The function
    allows mass change in individual layers and merges/splits layers as needed to prevent
    numerical instability.

    Args:
        prhofirn (array-like): Snow/firn density for each subsurface layer (kg/m³).
        psnowc (array-like): Snow content per layer (m w.e.).
        psnic (array-like): Ice content per layer (m w.e.).
        pslwc (array-like): Liquid water content per layer (m w.e.).
        pdgrain (array-like): Mean grain size per layer (mm).
        zdtime (float): Time step duration (s).

    Returns:
        tuple: Updated values for:
            - prhofirn: snow/firn density (kg/m³)
            - psnowc: snow content (m w.e.)
            - psnic: ice content (m w.e.)
            - pslwc: liquid water content (m w.e.)
            - pdgrain: grain size (mm)
            - zdtime: time step duration (s)
    
    Notes:
        - This version does not conserve total mass per layer; water redistributes
          and changes total layer mass.
        - Layer merging and splitting ensures numerical stability (e.g., no negative mass).
        - Originally developed in FORTRAN by Peter Langen and Robert S. Fausto,
          translated to Python by Baptiste Vandecrux.

    Reference:
        Based on water percolation and retention schemes described in:
        Hirashima et al. (2010), Colbeck (1975), and Coléou & Lesaffre (1998).
    """
    # ======================================================
    # Perform liquid water percolation and runoff processes
    # ======================================================

    # *** Potential Darcy fluxes calculated across interfaces ***
    # (if non-Darcy mode is used, fluxes will be calculated inline below)
    # The flux is limited such that no more water than available in the upper layer
    # is redistributed downward.

    if np.sum(pslwc) < 1e-12:
        return prhofirn, psnowc, psnic, pslwc, pdgrain, 0

    # if do_no_darcy==0:
    darcy_fluxes = calc_darcy_fluxes(pslwc, psnowc, psnic, pdgrain, prhofirn, zdtime)

    # Update BV2017: t_runoff is now calculated outside of the subsurface scheme
    # Bottom layer: Receive from above. Give excess liquid to runoff.
    jk = len(pslwc) - 1
    zrogl = 0
    # Remove runoff from layer (and add to runoff-box zrogl)

    liqmaxloc = CLliqF(prhofirn[jk])
    liqmaxM = liqmaxloc * 999.8395 / 900 * (900 / prhofirn[jk] - 1)
    potret = max(liqmaxM * psnowc[jk], 0)
    # Calculate liqexcess from layer water. Make sure it is not negative:
    liqexcess = max(pslwc[jk] - potret, 0)

    t_runoff = 28513.796102033262
    do_no_darcy = 0
    avoid_runoff = 0

    liqro = liqexcess / t_runoff * zdtime
    # Take runoff from water content and give to runoff box (Update PLA)
    zrogl = zrogl + liqro
    pslwc[jk] = pslwc[jk] - liqro

    for jk in range(len(pslwc) - 2, -1, -1):  # c.jpgrnd-1:-1:1
        # BV2017 removing percolation blocking criteria
        if ThetaF(pslwc[jk + 1], psnowc[jk + 1], prhofirn[jk + 1]) >= 1:
            # if next layer already saturated, then we do not allow flux
            # into that layer, but instead activate runoff
            dflux = 0
            do_runoff = 1
        else:
            if do_no_darcy:
                # Potential Darcy-like flux in case of do-no-Darcy,
                # basically just all liqexcess:
                # if (c.calc_CLliq):
                liqmaxloc = CLliqF(prhofirn[jk])
                # else:
                #    liqmaxloc = c.liqmax
                liqmaxM = liqmaxloc * 999.8395 / 900 * (900 / prhofirn[jk] - 1)
                potret = max(liqmaxM * psnowc[jk], 0)
                liqexcess = pslwc[jk] - potret
                darcy_fluxes[jk] = max(liqexcess, 0)

            # Calculate water in next layer, when this is at saturation (Theta = 1):
            plsat = psnowc[jk + 1] * 999.8395 / 900 * (900 / prhofirn[jk + 1] - 1)

            # Do not allow flux to be greater than plsat-pl in next layer.
            # Also, do not allow it to be negative
            dflux = max(min(darcy_fluxes[jk], plsat - pslwc[jk + 1]), 0)

            # Update BV2017: if no snow to hold liquid water then runoff is
            # turned on
            if (darcy_fluxes[jk] >= plsat - pslwc[jk + 1]) | (psnowc[jk] < 1e-12):
                # There is enough darcy flow to fill up layer below (and perhaps more)
                # so runoff should be activated (on what's left here after
                # we have moved down enough to fill layer below):
                do_runoff = 1
            else:
                # No runoff from this layer, since darcy flow is not large enough
                # to fill up layer below
                do_runoff = 0

        if dflux > 0:
            # Yes: Darcy flow to next layer
            # Update BV2017: Now temperature is not updated. water
            # percolating is at 0degC cannot change temperature until phase
            # change occurs in refreezing scheme
            pslwc[jk + 1] = pslwc[jk + 1] + dflux
            pslwc[jk] = pslwc[jk] - dflux

        if do_runoff & (~avoid_runoff):
            # Yes: Remove runoff from layer (and add to runoff-box zrogl)
            # PLA Darcy 2016
            # if c.calc_CLliq:
            liqmaxloc = CLliqF(prhofirn[jk])
            # else:
            #     liqmaxloc = c.liqmax

            liqmaxM = liqmaxloc * 999.8395 / 900 * (900 / prhofirn[jk] - 1)
            potret = max(liqmaxM * psnowc[jk], 0)

            # Calculate liqexcess from layer water. Make sure it is not negative:
            # Update BV2017: since dflux has already left slwc[jk] so no
            # need to remove it anymore (compared to old version)
            liqexcess = max(pslwc[jk] - potret, 0)

            # Update BV2017: for all layer, if there is no snow to hold
            # liquid water then runoff is instantaneous (not allowing
            # subglacial lakes to form
            if psnowc[jk] < 1e-12:
                # Here in layer: If there is no snow, run off immediately
                # in other word surface runoff is instantaneous
                liqro = liqexcess
            else:
                # Theta = ThetaF(pslwc[jk], psnowc[jk], prhofirn[jk])
                # liqro_darcy = kF(Theta,pdgrain[jk],prhofirn[jk],psnic[jk],psnowc[jk], c) * c.ElevGrad
                # old version based on Zuo and Oerlemans (1996)
                liqro_darcy = liqexcess / t_runoff * zdtime

                liqro = min(liqro_darcy, liqexcess)
                pore_space = psnowc[jk] * 999.8395 * (1 / prhofirn[jk] - 1 / 900)
                if pore_space < 1e-12: pore_space = 0

                # Update PLA
                if pslwc[jk] - liqro > pore_space:
                    liqro = pslwc[jk] - pore_space

            # Take runoff from water content and give to runoff box
            zrogl = zrogl + liqro

            pslwc[jk] = pslwc[jk] - liqro
    return prhofirn, psnowc, psnic, pslwc, pdgrain, zrogl
