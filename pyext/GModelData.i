/***************************************************************************
 *      GModelData.i  -  Abstract virtual data model class python I/F      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelData.i
 * @brief GModelData class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelData.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelData
 *
 * @brief Abstract virtual sky model class python interface defintion.
 ***************************************************************************/
class GModelData : public GModel {
public:
    // Constructors and destructors
    GModelData(void);
    explicit GModelData(const GXmlElement& xml);
    GModelData(const GModelData& model);
    virtual ~GModelData(void);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GModelData* clone(void) const = 0;
    virtual std::string type(void) const = 0;
    virtual std::string print(void) const = 0;
    virtual double      eval(const GEvent& event, const GObservation& obs) = 0;
    virtual double      eval_gradients(const GEvent& event, const GObservation& obs) = 0;
    virtual void        read(const GXmlElement& xml) = 0;
    virtual void        write(GXmlElement& xml) const = 0;

    // Implemented pure virtual methods
    int size(void) const;

    // Other methods

protected:
    virtual void set_pointers(void) = 0;
};


/***********************************************************************//**
 * @brief GModelData class extension
 ***************************************************************************/
%extend GModelData {
/*
    GModelPar __getitem__(int index) {
    if (index >= 0 && index < self->size())
        return (*self)(index);
    else
        throw GException::out_of_range("__getitem__(int)", index, self->size());
    }
    void __setitem__(int index, const GModelPar& val) {
        if (index>=0 && index < self->size())
            (*self)(index) = val;
        else
            throw GException::out_of_range("__setitem__(int)", index, self->size());
    }
*/
};