/***************************************************************************
 *     GModelSpectralExpPlaw.i  -  Exponential cut off power law model     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpectralExpPlaw.i
 * @brief GModelSpectralExpPlaw class SWIG interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralExpPlaw.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralExpPlaw
 *
 * @brief Implements an exponentially cut off power law.
 ***************************************************************************/
class GModelSpectralExpPlaw  : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralExpPlaw(void);
    explicit GModelSpectralExpPlaw(const double& norm, const double& index,
                                   const double& ecut);
    explicit GModelSpectralExpPlaw(const GXmlElement& xml);
    GModelSpectralExpPlaw(const GModelSpectralExpPlaw& model);
    virtual ~GModelSpectralExpPlaw(void);

    // Implemented pure virtual methods
    void                   clear(void);
    GModelSpectralExpPlaw* clone(void) const;
    int                    size(void) const { return m_npars; }
    std::string            type(void) const { return "ExpCutoff"; }
    double                 eval(const GEnergy& srcEng);
    double                 eval_gradients(const GEnergy& srcEng);
    void                   read(const GXmlElement& xml);
    void                   write(GXmlElement& xml) const;

    // Other methods
    void   autoscale(void);
    double norm(void) const { return m_norm.real_value(); }
    double index(void) const { return m_index.real_value(); }
    double ecut(void) const { return m_ecut.real_value(); }
    double pivot(void) const { return m_pivot.real_value(); }
};


/***********************************************************************//**
 * @brief GModelSpectralExpPlaw class extension
 ***************************************************************************/
%extend GModelSpectralExpPlaw {
    char *__str__() {
        static std::string result = self->print();
        return ((char*)result.c_str());
    }
    GModelSpectralExpPlaw copy() {
        return (*self);
    }
};
