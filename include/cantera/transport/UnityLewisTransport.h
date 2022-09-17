/**
 *  @file UnityLewisTransport.h
 *    Headers for the UnityLewisTransport object, which models transport
 *    properties in ideal gas solutions using the unity Lewis number
 *    approximation
 *    (see \ref tranprops and \link Cantera::UnityLewisTransport UnityLewisTransport \endlink) .
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_UNITYLEWISTRAN_H
#define CT_UNITYLEWISTRAN_H

#include "MixTransport.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{
//! Class UnityLewisTransport implements the unity Lewis number approximation
//! for the mixture-averaged species diffusion coefficients. Mixture-averaged
//! transport properties for viscosity and thermal conductivity are inherited
//! from the MixTransport class.
//! @ingroup tranprops
class UnityLewisTransport : public MixTransport
{
public:
//    UnityLewisTransport() {}

    virtual std::string transportType() const {
        return "UnityLewis";
    }

    //! Returns the unity Lewis number approximation based diffusion
    //! coefficients [m^2/s].
    /*!
     * Returns the unity Lewis number approximation based diffusion coefficients
     * for a gas, appropriate for calculating the mass averaged diffusive flux
     * with respect to the mass averaged velocity using gradients of the mole
     * fraction.
     *
     * \f[
     *     D^\prime_{km} = \frac{\lambda}{\rho c_p}
     * \f]
     *
     * In order to obtain the expected behavior from a unity Lewis number model,
     * this formulation requires that the correction velocity be computed as
     *
     * \f[
     *     V_c = \sum \frac{W_k}{\overline{W}} D^\prime_{km} \nabla X_k
     * \f]
     *
     * @param[out] d  Vector of diffusion coefficients for each species (m^2/s).
     * length m_nsp.
     */
    virtual void getMixDiffCoeffs(double* const d) {
        double Dm = thermalConductivity() / (m_thermo->density() * m_thermo->cp_mass());
        for (size_t k = 0; k < m_nsp; k++) {
            d[k] = Dm;
        }
    }

    //! Not implemented for unity Lewis number approximation
    virtual void getMixDiffCoeffsMole(double* const d){
        throw NotImplementedError("UnityLewisTransport::getMixDiffCoeffsMole");
    }

    //! Returns the unity Lewis number approximation based diffusion
    //! coefficients [m^2/s].
    /*!
     * These are the coefficients for calculating the diffusive mass fluxes
     * from the species mass fraction gradients, computed as
     *
     * \f[
     *     D_{km} = \frac{\lambda}{\rho c_p}
     * \f]
     *
     * @param[out] d  Vector of diffusion coefficients for each species (m^2/s).
     * length m_nsp.
     */
    virtual void getMixDiffCoeffsMass(double* const d){
        double Dm = thermalConductivity() / (m_thermo->density() * m_thermo->cp_mass());
        double Le_i[53] ={
            0.27443569898605347,
            0.15385763347148895,
            0.5543738603591919,
            0.8452105522155762,
            0.56460040807724,
            0.713638424873352,
            0.8536510467529297,
            0.8592507839202881,
            0.5604416728019714,
            0.5199534893035889,
            0.7710865139961243,
            0.7710865139961243,
            0.7878470420837402,
            0.790300726890564,
            0.8603102564811707,
            1.0825186967849731,
            1.0451807975769043,
            1.053497314453125,
            1.068619966506958,
            1.068619966506958,
            1.074398398399353,
            1.027429461479187,
            1.0375499725341797,
            1.0472034215927124,
            1.059842586517334,
            1.1490812301635742,
            1.1581830978393555,
            0.6890178918838501,
            1.2294784784317017,
            1.2294784784317017,
            0.5888339281082153,
            0.5265711545944214,
            0.5083094239234924,
            0.74753338098526,
            0.8286235332489014,
            0.8666804432868958,
            0.9892193675041199,
            1.094902515411377,
            0.8575202822685242,
            0.8238334655761719,
            1.0616075992584229,
            1.070987582206726,
            0.6890308856964111,
            1.0902538299560547,
            1.0902538299560547,
            1.0902538299560547,
            1.085353970527649,
            0.7689386606216431,
            0.8779736757278442,
            1.4854788780212402,
            1.4919153451919556,
            1.2350292205810547,
            1.2404004335403442,
        };

        for (size_t k = 0; k < m_nsp; k++) {
            d[k] = Dm;
        }
    }
};
}
#endif
