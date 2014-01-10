/***************************************************************************
 *         GModelSpectralGauss.hpp - Spectral Gauss model class            *
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
 * @file GModelSpectralGauss.hpp
 * @brief Gauss spectral model class interface definition
 * @author Christoph Deil
 */

#ifndef GMODELSPECTRALGAUSS_HPP
#define GMODELSPECTRALGAUSS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelPar.hpp"
#include "GModelSpectral.hpp"
#include "GEnergy.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelSpectralGauss
 *
 * @brief Gauss spectral model class
 *
 * This class implements a Gauss spectrum. The model is defined by
 *
 * \f[
 *    S_{\rm E}(E | t) = {\tt m\_norm}
 * \f]
 *
 * where
 * \f${\tt m\_norm}\f$ is the normalization constant in units of
 * ph/cm2/s/MeV.
 ***************************************************************************/
class GModelSpectralGauss : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralGauss(void);
    explicit GModelSpectralGauss(const GXmlElement& xml);
    explicit GModelSpectralGauss(const double& value);
    GModelSpectralGauss(const GModelSpectralGauss& model);
    virtual ~GModelSpectralGauss(void);

    // Operators
    virtual GModelSpectralGauss& operator=(const GModelSpectralGauss& model);

    // Implemented pure virtual base class methods
    virtual void                 clear(void);
    virtual GModelSpectralGauss* clone(void) const;
    virtual std::string          type(void) const;
    virtual double               eval(const GEnergy& srcEng,
                                      const GTime&   srcTime) const;
    virtual double               eval_gradients(const GEnergy& srcEng,
                                                const GTime&   srcTime);
    virtual double               flux(const GEnergy& emin,
                                      const GEnergy& emax) const;
    virtual double               eflux(const GEnergy& emin,
                                       const GEnergy& emax) const;
    virtual GEnergy              mc(const GEnergy& emin,
                                    const GEnergy& emax,
                                    const GTime&   time,
                                    GRan&          ran) const;
    virtual void                 read(const GXmlElement& xml);
    virtual void                 write(GXmlElement& xml) const;
    virtual std::string          print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double value(void) const;
    void   value(const double& value);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralGauss& model);
    void free_members(void);

    // Protected members
    GModelPar m_norm;  //!< Normalization factor
};


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "ConstantValue".
 *
 * Returns the type of the Gauss spectral model.
 ***************************************************************************/
inline
std::string GModelSpectralGauss::type(void) const
{
    return "Gauss";
}


/***********************************************************************//**
 * @brief Return model value
 *
 * @return Model value (ph/cm2/s/MeV).
 *
 * Returns the model value.
 ***************************************************************************/
inline
double GModelSpectralGauss::value(void) const
{
    return (m_norm.value());
}


/***********************************************************************//**
 * @brief Set model value
 *
 * @param[in] value Model value (ph/cm2/s/MeV).
 *
 * Sets the model value.
 ***************************************************************************/
inline
void GModelSpectralGauss::value(const double& value)
{
    m_norm.value(value);
    return;
}

#endif /* GMODELSPECTRALGAUSS_HPP */
