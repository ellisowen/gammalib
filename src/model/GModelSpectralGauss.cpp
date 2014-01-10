/***************************************************************************
 *         GModelSpectralGauss.cpp - Spectral Gaussian model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2013 by Juergen Knoedlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpectralGauss.cpp
 * @brief Gaussian spectral model class implementation
 * @author Christoph Deil & Ellis Owen
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GModelSpectralGauss.hpp"
#include "GModelSpectralRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpectralGauss    g_spectral_gauss_seed;
const GModelSpectralRegistry g_spectral_gauss_registry(&g_spectral_gauss_seed);

/* __ Method name definitions ____________________________________________ */
#define G_FLUX                "GModelSpectralGauss::flux(GEnergy&, GEnergy&)"
#define G_EFLUX              "GModelSpectralGauss::eflux(GEnergy&, GEnergy&)"
#define G_MC     "GModelSpectralGauss::mc(GEnergy&, GEnergy&, GTime&, GRan&)"
#define G_READ                      "GModelSpectralGauss::read(GXmlElement&)"
#define G_WRITE                    "GModelSpectralGauss::write(GXmlElement&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GModelSpectralGauss::GModelSpectralGauss(void) : GModelSpectral()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs constant spectral model by extracting information from an XML
 * element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpectralGauss::GModelSpectralGauss(const GXmlElement& xml) :
                     GModelSpectral()
{
    // Initialise members
    init_members();

    // Read information from XML element
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] prefactor Power law pre factor (in ph/cm2/s/MeV).
 * @param[in] mean Mean energy.
 * @param[in] sigma Energy width.
 ***************************************************************************/
GModelSpectralGauss::GModelSpectralGauss(const double&  prefactor,
                                         const GEnergy& mean,
                                         const GEnergy& sigma) :
                           GModelSpectral()
{
    // Initialise members
    init_members();

    // Set parameters
    m_norm.value(prefactor);
    m_mean.value(mean.MeV());
    m_sigma.value(sigma.MeV());

    // Autoscale parameters
    autoscale();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Spectral Gaussian model.
 ***************************************************************************/
GModelSpectralGauss::GModelSpectralGauss(const GModelSpectralGauss& model) :
                     GModelSpectral(model)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(model);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GModelSpectralGauss::~GModelSpectralGauss(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] model Gauss spectral model.
 * @return Gauss spectral model.
 ***************************************************************************/
GModelSpectralGauss& GModelSpectralGauss::operator=(const GModelSpectralGauss& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpectral::operator=(model);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(model);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear Gaussian spectral model
 ***************************************************************************/
void GModelSpectralGauss::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelSpectral::free_members();

    // Initialise members
    this->GModelSpectral::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone Gaussian spectral model
 *
 * @return Pointer to deep copy of Gaussian spectral model.
 ***************************************************************************/
GModelSpectralGauss* GModelSpectralGauss::clone(void) const
{
    // Clone Gaussian spectral model
    return new GModelSpectralGauss(*this);
}


/***********************************************************************//**
 * @brief Evaluate model value
 *
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @return Model value (ph/cm2/s/MeV).
 *
 * Evaluates
 *
 * TODO: write formula
 ***************************************************************************/
double GModelSpectralGauss::eval(const GEnergy& srcEng,
                                 const GTime&   srcTime) const
{
    double energy = srcEng.MeV();
    double norm = m_norm.value();
    double mean = m_mean.value();
    double sigma = m_sigma.value();

    // Compute function value
    double term1 = (norm / sigma) * gammalib::inv_sqrt2pi;
    double term2 = (energy - mean) * (energy - mean) / (2 * sigma * sigma);
    double value = term1 * std::exp(- term2);

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate model value and gradient
 *
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @return Model value (ph/cm2/s/MeV).
 *
 * This method simply calls the eval() method as no analytical gradients will
 * be computed. See the eval() method for details.
 ***************************************************************************/
double GModelSpectralGauss::eval_gradients(const GEnergy& srcEng,
                                           const GTime&   srcTime)
{
    // Return
    return eval(srcEng, srcTime);
}


/***********************************************************************//**
 * @brief Returns model photon flux between [emin, emax] (ph/cm2/s)
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @return Photon flux (ph/cm2/s).
 *
 * Computes
 *
 * \f[
 *    \int_{\tt emin}^{\tt emax} S_{\rm E}(E | t) dE
 * \f]
 *
 * where
 * - [@p emin, @p emax] is an energy interval, and
 * - \f$S_{\rm E}(E | t)\f$ is the spectral model (ph/cm2/s/MeV).
 * The integration is done analytically.
 ***************************************************************************/
double GModelSpectralGauss::flux(const GEnergy& emin,
                                 const GEnergy& emax) const
{
    // TODO: implement as in Gaussian::integral from Gaussian.cxx in the Fermi ScienceTools
    // Initialise flux
    double flux = 0.0;
    
    // Compute only if integration range is valid
    if (emin < emax) {

        // Compute flux for a constant model
        flux = m_norm.value() * (emax.MeV() - emin.MeV());
    
    } // endif: integration range was valid

    // Return
    return flux;
}


/***********************************************************************//**
 * @brief Returns model energy flux between [emin, emax] (erg/cm2/s)
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @return Energy flux (erg/cm2/s).
 *
 * Computes
 *
 * \f[
 *    \int_{\tt emin}^{\tt emax} S_{\rm E}(E | t) E \, dE
 * \f]
 *
 * where
 * - [@p emin, @p emax] is an energy interval, and
 * - \f$S_{\rm E}(E | t)\f$ is the spectral model (ph/cm2/s/MeV).
 * The integration is done analytically.
 ***************************************************************************/
double GModelSpectralGauss::eflux(const GEnergy& emin,
                                  const GEnergy& emax) const
{
    // TODO Implement using numerical integration as in GModelSpectralExpPlaw::eflux
    // Initialise flux
    double eflux = 0.0;
    
    // Compute only if integration range is valid
    if (emin < emax) {

        // Compute flux for a constant model
        eflux = m_norm.value() * 0.5 * (emax.MeV()*emax.MeV() - 
                                        emin.MeV()*emin.MeV());

        // Convert from MeV/cm2/s to erg/cm2/s
        eflux *= gammalib::MeV2erg;
    
    } // endif: integration range was valid

    // Return
    return eflux;
}


/***********************************************************************//**
 * @brief Returns MC energy between [emin, emax]
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @param[in] time True photon arrival time.
 * @param[in,out] ran Random number generator.
 * @return Energy.
 *
 * @exception GException::erange_invalid
 *            Energy range is invalid (emin < emax required).
 *
 * Returns Monte Carlo energy by randomly drawing from a constant between
 * the minimum and maximum photon energy.
 ***************************************************************************/
GEnergy GModelSpectralGauss::mc(const GEnergy& emin,
                                const GEnergy& emax,
                                const GTime&   time,
                                GRan&          ran) const
{
    // Throw an exception if energy range is invalid
    if (emin >= emax) {
        throw GException::erange_invalid(G_MC, emin.MeV(), emax.MeV(),
              "Minimum energy < maximum energy required.");
    }

    // TODO: implement
    // Allocate energy
    GEnergy energy;

    // Get uniform value between 0 and 1
    double u = ran.uniform();

    // Map into [emin,emax] range
    energy.MeV((emax.MeV() - emin.MeV()) * u + emin.MeV());

    // Return energy
    return energy;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element containing Gaussian model information.
 *
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter name found in XML element.
 *
 * Read the spectral Gaussian information from an XML element.
 * The format of the XML elements is:
 *
 *     <spectrum type="Gaussian">
 *       <parameter name="Prefactor" scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Mean"      scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Sigma"     scale=".." value=".." min=".." max=".." free=".."/>
 *     </spectrum>
 *
 ***************************************************************************/
void GModelSpectralGauss::read(const GXmlElement& xml)
{
  const int n_pars = 3;

  // Verify that XML element has exactly 3 parameters
  if (xml.elements() != n_pars || xml.elements("parameter") != n_pars) {
      throw GException::model_invalid_parnum(G_READ, xml,
            "Gaussian model requires exactly 3 parameters.");
  }

  // Extract model parameters
  int npar[] = {0, 0, 0};
  for (int i = 0; i < n_pars; ++i) {

      // Get parameter element
      const GXmlElement* par = xml.element("parameter", i);

      // Handle prefactor
      if (par->attribute("name") == "Prefactor") {
          m_norm.read(*par);
          npar[0]++;
      }

      // Handle mean
      else if (par->attribute("name") == "Mean") {
          m_mean.read(*par);
          npar[1]++;
      }

      // Handle sigma
      else if (par->attribute("name") == "Sigma") {
          m_sigma.read(*par);
          npar[2]++;
      }

  } // endfor: looped over all parameters

  // Verify that all parameters were found
  if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1) {
      throw GException::model_invalid_parnames(G_READ, xml,
            "Require \"Prefactor\", \"Mean\", and \"Sigma\""
            " parameters.");
  }

  // Return
  return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * @exception GException::model_invalid_spectral
 *            Existing XML element is not of type "ConstantValue"
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters or nodes found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Writes the spectral information into an XML element with the format
 *
 *     <spectrum type="ConstantValue">
 *       <parameter name="Value" scale="1" min="0" max="1000" value="1" free="1"/>
 *     </spectrum>
 ***************************************************************************/
void GModelSpectralGauss::write(GXmlElement& xml) const
{
    // TODO: implement as in GModelSpectralExpPlaw::write

    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", "ConstantValue");
    }

    // Verify model type
    if (xml.attribute("type") != "ConstantValue") {
        throw GException::model_invalid_spectral(G_WRITE, xml.attribute("type"),
              "Spectral model is not of type \"ConstantValue\".");
    }

    // If XML element has 0 nodes then append parameter node. The name
    // of the node is "Value" as this is the Fermi-LAT standard.
    // We thus assure that the XML files will be compatible with
    // Fermi-LAT.
    if (xml.elements() == 0) {
        xml.append(GXmlElement("parameter name=\"Value\""));
    }

    // Verify that XML element has exactly 1 parameter
    if (xml.elements() != 1 || xml.elements("parameter") != 1) {
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Spectral constant requires exactly 1 parameter.");
    }

    // Get parameter element
    GXmlElement* par = xml.element("parameter", 0);

    // Set parameyter
    if (par->attribute("name") == "Normalization" ||
        par->attribute("name") == "Value") {
        m_norm.write(*par);
    }
    else {
        throw GException::model_invalid_parnames(G_WRITE, xml,
                          "Spectral constant requires either"
                          " \"Normalization\" or \"Value\" parameter.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print spectral model information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing spectral model information.
 ***************************************************************************/
std::string GModelSpectralGauss::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpectralGauss ===");

        // Append model content
        result.append("\n"+gammalib::parformat("Number of parameters"));
        result.append(gammalib::str(size()));
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+m_pars[i]->print(chatter));
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelSpectralGauss::init_members(void)
{
    // Initialise constant normalisation
    m_norm.clear();
    m_norm.name("Value");
    m_norm.scale(1.0);
    m_norm.value(1.0);         // default: 1.0
    m_norm.range(0.0, 1000.0); // range:   [0, 1000]
    m_norm.free();
    m_norm.gradient(0.0);
    m_norm.has_grad(true);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Spectral constant model.
 ***************************************************************************/
void GModelSpectralGauss::copy_members(const GModelSpectralGauss& model)
{
    // Copy members
    m_norm = model.m_norm;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralGauss::free_members(void)
{
    // Return
    return;
}
