/***************************************************************************
 *               GLATResponse.hpp  -  GLAST LAT Response class             *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GLATResponse.hpp
 * @brief GLATResponse class interface definition.
 * @author J. Knodlseder
 */

#ifndef GLATRESPONSE_HPP
#define GLATRESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include "GNodeArray.hpp"
#include "GVector.hpp"
#include "GSkyDir.hpp"
#include "GResponse.hpp"
#include "GLATResponseTable.hpp"
#include "GFits.hpp"
#include "GFitsHDU.hpp"

/* __ Namespaces _________________________________________________________ */


/***********************************************************************//**
 * @class GLATResponse
 *
 * @brief Interface for the GLAST LAT instrument response function classes.
 ***************************************************************************/
class GLATResponse : public GResponse {

public:
    // Constructors and destructors
    GLATResponse();
    GLATResponse(const GLATResponse& rsp);
    ~GLATResponse();

    // Operators
    GLATResponse& operator= (const GLATResponse & rsp);

    // IRF Methods
    double irf(const GSkyDir& obsDir, const double& obsEng,
               const GSkyDir& srcDir, const double& srcEng,
               const GSkyDir& instPntDir, const double& instPosAng,
               const double& time);

    // Aeff Methods
    double aeff(const GSkyDir& obsDir, const double& obsEng,
                const GSkyDir& srcDir, const double& srcEng,
                const GSkyDir& instPntDir, const double& instPosAng,
                const double& time);

    // PSF Methods
    double  psf(const GSkyDir& obsDir, const double& obsEng,
                const GSkyDir& srcDir, const double& srcEng,
                const GSkyDir& instPntDir, const double& instPosAng,
                const double& time);
    double  psf(const double& delta, const double& logE, const double& ctheta);
    GVector psf(const GVector& delta, const double& logE, const double& ctheta);

    // Edisp Methods
    double edisp(const GSkyDir& obsDir, const double& obsEng,
                 const GSkyDir& srcDir, const double& srcEng,
                 const GSkyDir& instPntDir, const double& instPosAng,
                 const double& time);

    // Other Methods
    void          set_caldb(const std::string& caldb);
    void          load(const std::string& rspname, const std::string& rsptype);
    void          save(const std::string& rspname) const;
    GLATResponse* clone(void) const;

private:
    // Private Effective Area methods
    void    aeff_init(void);
    void    aeff_append(GFits& file) const;
    void    aeff_init_members(void);
    void    aeff_copy_members(const GLATResponse& rsp);
    void    aeff_free_members(void);

    // Private PSF methods
    void    psf_init(void);
    void    psf_append(GFits& file) const;
    double  psf_scale(const double& energy) const;
    double  psf_scale_log(const double& logE) const;
    double  psf_base_value(const double& u, const double& gamma) const;
    GVector psf_base_vector(const GVector& u, const double& gamma) const;
    void    psf_init_members(void);
    void    psf_copy_members(const GLATResponse& rsp);
    void    psf_free_members(void);

    // Private energy dissipation methods
    void    edisp_init(void);
    void    edisp_append(GFits& file) const;
    void    edisp_init_members(void);
    void    edisp_copy_members(const GLATResponse& rsp);
    void    edisp_free_members(void);

    // Other private methods
    void    init_members(void);
    void    copy_members(const GLATResponse& rsp);
    void    free_members(void);
    GVector get_fits_vector(const GFitsHDU* hdu, const std::string& colname, int row = 0);

    // Private Aeff data
    GLATResponseTable m_aeff_bins;       //!< Aeff energy and cos theta binning

    // Private PSF data
    GLATResponseTable m_psf_bins;        //!< PSF energy and cos theta binning
    GNodeArray        m_psf_angle;       //!< PSF vector binning
    double            m_psf_angle_max;   //!< PSF maximum angular distance covered
    double*           m_psf;             //!< PSF vector array
    double*           m_norm;            //!< PSF normalization array
    double*           m_sigma;           //!< PSF sigma array

    // Private Edisp data
    GLATResponseTable m_edisp_bins;      //!< Edisp energy and cos theta binning

    // Other private data
    std::string m_rsptype;   //!< Response type ('front','back')
    int         m_section;   //!< PSF section (0=front, 1=back)

};

#endif /* GLATRESPONSE_HPP */
