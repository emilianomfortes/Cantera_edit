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
        double Le_i[10] ={
            0.32400885224342346,
            0.18087716400623322,
            0.6437716484069824,
            0.9815124273300171,
            0.6555334329605103,
            0.9272491931915283,
            0.9917637705802917,
            0.9982198476791382,
            1.021479845046997,
            0.9413276314735413,
        };

        for (size_t k = 0; k < m_nsp; k++) {
            d[k] = Dm/Le_i[k];
        }
    }
};
}
#endif
