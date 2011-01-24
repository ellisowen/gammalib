/***************************************************************************
 *              GGti.i  -  Good time interval class python I/F             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GGti.i
 * @brief GGti class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GGti.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GGti
 *
 * @brief Interface for the GTI class.
 *
 * This class holds a list of Good Time Intervals, i.e. time intervals that
 * are valid for science analysis. Times are stored using the GTime class.
 ***************************************************************************/
class GGti {
public:
    // Constructors and destructors
    GGti(void);
    GGti(const GGti& gti);
    virtual ~GGti(void);

    // Methods
    void   clear(void);
    void   add(const GTime& tstart, const GTime& tstop);
    void   append(const GTime& tstart, const GTime& tstop);
    void   insert(const GTime& tstart, const GTime& tstop);
	void   load(const std::string& filename,
                const std::string& extname = "GTI");
	void   save(const std::string& filename, bool clobber,
                const std::string& extname = "GTI") const;
    void   read(GFitsTable* hdu);
    void   write(GFits* file, const std::string& extname = "GTI") const;
    int    size(void) const { return m_num; }
	GTime  tstart(void) const { return m_tstart; }
	GTime  tstop(void) const { return m_tstop; }
	GTime  tstart(int inx) const;
	GTime  tstop(int inx) const;
	double telapse(void) const { return m_telapse; }
	double ontime(void) const { return m_ontime; }
    bool   isin(const GTime& t) const;
};


/***********************************************************************//**
 * @brief GGti class extension
 ***************************************************************************/
%extend GGti {
    /*
    char *__str__() {
        return tochar(self->print());
    }
    */
    GGti copy() {
        return (*self);
    }
};
