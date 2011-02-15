/***************************************************************************
 *                        GEnergy.hpp - Energy class                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GEnergy.hpp
 * @brief Energy value class definition.
 * @author J. Knodlseder
 */

#ifndef GENERGY_HPP
#define GENERGY_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"


/***********************************************************************//**
 * @class GEnergy
 *
 * @brief Class that handles energies in a unit independent way.
 *
 * The GEnergy class stores an energy value in units of MeV and implements
 * methods that provide automatic conversion of the energy values in other
 * units. This makes instrument specific implementations more robust and
 * reduces the risk of unit errors. The method also holds on request the
 * log10 of the energy in MeV, which is only recomputed when the energy
 * value is changed. This considerably reduces the number of log10 operations
 * in many applications and speeds up the performance.
 ***************************************************************************/
class GEnergy {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GEnergy& eng);
    friend GLog&         operator<< (GLog& log, const GEnergy& eng);

    // Operator friends
    friend GEnergy operator+ (const GEnergy &a, const GEnergy &b);
    friend GEnergy operator- (const GEnergy &a, const GEnergy &b);
    friend GEnergy operator* (const double &a, const GEnergy &b);
    friend GEnergy operator* (const GEnergy &a, const double &b);
    friend GEnergy operator/ (const GEnergy &a, const double &b);
    friend bool    operator== (const GEnergy &a, const GEnergy &b);
    friend bool    operator!= (const GEnergy &a, const GEnergy &b);
    friend bool    operator< (const GEnergy &a, const GEnergy &b);
    friend bool    operator<= (const GEnergy &a, const GEnergy &b);
    friend bool    operator> (const GEnergy &a, const GEnergy &b);
    friend bool    operator>= (const GEnergy &a, const GEnergy &b);

public:
    // Constructors and destructors
    GEnergy(void);
    GEnergy(const GEnergy& eng);
    virtual ~GEnergy(void);
 
    // Operators
    GEnergy& operator= (const GEnergy& eng);
    GEnergy& operator+= (const GEnergy& eng);
    GEnergy& operator-= (const GEnergy& eng);

    // Methods
    void        clear(void);
    double      erg(void) const;
    double      keV(void) const;
    double      MeV(void) const;
    double      GeV(void) const;
    double      TeV(void) const;
    double      log10keV(void) const;
    double      log10MeV(void) const;
    double      log10GeV(void) const;
    double      log10TeV(void) const;
    void        erg(const double& eng);
    void        keV(const double& eng);
    void        MeV(const double& eng);
    void        GeV(const double& eng);
    void        TeV(const double& eng);
    void        log10keV(const double& eng);
    void        log10MeV(const double& eng);
    void        log10GeV(const double& eng);
    void        log10TeV(const double& eng);
    std::string print(void) const;
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GEnergy& eng);
    void free_members(void);

    // Protected data members
    double         m_energy;      //!< Energy in MeV
    mutable double m_elog10;      //!< log10 of energy in MeV
    mutable bool   m_has_log10;   //!< log10 of energy is valid
};

/***************************************************************************
 *                               Inline friends                            *
 ***************************************************************************/
inline
GEnergy& GEnergy::operator+= (const GEnergy& eng)
{
    m_energy    += eng.m_energy;
    m_has_log10  = false;
    return *this;
}
inline
GEnergy& GEnergy::operator-= (const GEnergy& eng)
{
    m_energy     += eng.m_energy;
    m_has_log10  = false;
    return *this;
}
inline
GEnergy operator+ (const GEnergy& a, const GEnergy& b)
{
    GEnergy result;
    result.m_energy = a.m_energy + b.m_energy;
    return result;
}
inline
GEnergy operator- (const GEnergy& a, const GEnergy& b)
{
    GEnergy result;
    result.m_energy = a.m_energy - b.m_energy;
    return result;
}
inline
GEnergy operator* (const double& a, const GEnergy& b)
{
    GEnergy result;
    result.m_energy = a * b.m_energy;
    return result;
}
inline
GEnergy operator* (const GEnergy& a, const double& b)
{
    GEnergy result;
    result.m_energy = b * a.m_energy;
    return result;
}
inline
GEnergy operator/ (const GEnergy& a, const double& b)
{
    GEnergy result;
    result.m_energy = a.m_energy / b;
    return result;
}
inline
bool operator== (const GEnergy &a, const GEnergy &b)
{
    return (a.m_energy == b.m_energy);
}
inline
bool operator!= (const GEnergy &a, const GEnergy &b)
{
    return (a.m_energy != b.m_energy);
}
inline
bool operator< (const GEnergy &a, const GEnergy &b)
{
    return (a.m_energy < b.m_energy);
}
inline
bool operator<= (const GEnergy &a, const GEnergy &b)
{
    return (a.m_energy <= b.m_energy);
}
inline
bool operator> (const GEnergy &a, const GEnergy &b)
{
    return (a.m_energy > b.m_energy);
}
inline
bool operator>= (const GEnergy &a, const GEnergy &b)
{
    return (a.m_energy >= b.m_energy);
}

#endif /* GENERGY_HPP */
