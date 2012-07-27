/***************************************************************************
 *           GTestEventList.hpp  -  Test event atom container class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Jean-Baptiste Cayrou                             *
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

#ifndef GTESTEVENTLIST_HPP
#define GTESTEVENTLIST_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GEventList.hpp"
#include "GTestEventAtom.hpp"
#include "GTestRoi.hpp"

class GTestEventList : public GEventList {

public:
    // Constructors and destructors
    GTestEventList(void): GEventList() {
        init_members();
        return;
    }
    GTestEventList(const GTestEventList& list) : GEventList(list){
        init_members();
        copy_members(list);
        return;
    }
    virtual ~GTestEventList(void){
        free_members();
        return;
    }

    // Operators
    virtual GTestEventList&       operator=(const GTestEventList& list){
        // Execute only if object is not identical
        if (this != &list) {

             // Copy base class members
            this->GEventList::operator=(list);

            // Free members
            free_members();

             // Initialise members
            init_members();

            // Copy members
            copy_members(list);

        } // endif: object was not identical

        // Return this object
        return *this;
    }
    
    virtual GTestEventAtom* operator[](const int& index){
            // Optionally check if the index is valid
        #if defined(G_RANGE_CHECK)
            if (index < 0 || index >= size())
                throw GException::out_of_range(G_OPERATOR, index, 0, size()-1);
        #endif
    
        // Return pointer
        return (&(m_events[index]));
    }
    virtual const GTestEventAtom* operator[](const int& index) const {
        // Optionally check if the index is valid
        #if defined(G_RANGE_CHECK)
            if (index < 0 || index >= size())
                throw GException::out_of_range(G_OPERATOR, index, 0, size()-1);
        #endif
        
            // Return pointer
            return (&(m_events[index]));
    }

    // Implemented pure virtual base class methods
    virtual void           clear(void){
        free_members();
        this->GEventList::free_members();
        this->GEvents::free_members();

        // Initialise members
        this->GEvents::init_members();
        this->GEventList::init_members();
        init_members();
        
        return;
    }
    virtual GTestEventList* clone(void) const{
        return new GTestEventList(*this);
    }
    virtual int            size(void) const { return m_events.size(); }
    virtual void           load(const std::string& filename){};
    virtual void           save(const std::string& filename, bool clobber = false) const{};
    virtual void           read(const GFits& file){}
    virtual void           write(GFits& file) const{}
    virtual int            number(void) const { return m_events.size(); }
    virtual void           roi(const GRoi& roi){
            // Cast ROI dynamically
            const GTestRoi* ptr = dynamic_cast<const GTestRoi*>(&roi);
    
             // Throw exception if ROI is not of correct type
            if (ptr == NULL)
                throw;
    
             // Set ROI
            m_roi = *ptr;
            
            return;
    }
    
    virtual const GTestRoi& roi(void) const { return m_roi; }
    
    std::string print(void) const{
        
          // Initialise result string
        std::string result;
        
        // Append header
        result.append("=== GTestEventList ===");
        result.append("\n"+parformat("Number of events")+str(number()));
        
        return result;
    
    }

        // Implement other methods
    void append(const GTestEventAtom& event){
        m_events.push_back(event);
        return;
    }
    

protected:
    // Protected methods
    void init_members(void){
        m_roi.clear();
        m_events.clear();
        return;
    }
    
    void copy_members(const GTestEventList& list){
            m_roi = list.m_roi;
            m_events = list.m_events;
            return;
    }
    void free_members(void){ return; }
    
    virtual void set_energies(void) { return; }
    virtual void set_times(void) { return; }
    
    GTestRoi                    m_roi;     //!< Region of interest
    std::vector<GTestEventAtom> m_events;  //!< Events
};

#endif /* GTestEVENTLIST_HPP */
