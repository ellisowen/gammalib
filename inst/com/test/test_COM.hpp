/***************************************************************************
 *                  test_COM.hpp  -  Test COMPTEL classes                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file test_COM.hpp
 * @brief Definition of COMTEL classes unit tests
 * @author Juergen Knoedlseder
 */

#ifndef TEST_COM_HPP
#define TEST_COM_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCOMLib.hpp"


/***********************************************************************//**
 * @class TestGCOMResponse
 *
 * @brief Test suite for COMPTEL response class testing
 *
 * This class defines a unit test suite for the COM response class.
 ***************************************************************************/
class TestGCOMResponse : public GTestSuite {
public:
    // Constructors and destructors
    TestGCOMResponse(void) : GTestSuite() {}
    virtual ~TestGCOMResponse(void) {}

    // Methods
    virtual void set(void);
    void         test_iaq_response(void);
};


/***********************************************************************//**
 * @class TestGCOMObservation
 *
 * @brief Test suite for COMPTEL observation testing
 *
 * This class defines a unit test suite for testing of GCOMObservation.
 ***************************************************************************/
class TestGCOMObservation : public GTestSuite {
public:
    // Constructors and destructors
    TestGCOMObservation(void) : GTestSuite() {}
    virtual ~TestGCOMObservation(void) {}

    // Methods
    virtual void set(void);
    void         test_binned_obs(void);
    void         test_event_bin(void);
    void         test_event_cube(void);
};

#endif /* TEST_COM_HPP */
