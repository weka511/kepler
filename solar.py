# Copyright (C) 2015-2017 Greenweaves Software Pty Ltd

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this software.  If not, see <http://www.gnu.org/licenses/>

'''
Model for solar irradiation. References:
    Solar Radiation on Mars, Joseph Appelbaum & Dennis Flood,
    Lewis Research Center, NASA
    GA Goosse H., P.Y. Barriat, W. Lefebvre, M.F. Loutre and V. Zunz, (2008-2010).
    Introduction to climate dynamics and climate modeling.
    Online textbook available at http://www.climate.be/textbook.
'''
import math as m

def sin_declination(obliquity,true_longitude):
    '''
    Sine of declination
    Appelbaum & Flood equation (7)
    Goosse et al equation (2.22)
    Parameters:
         obliquity
         true_longitude        
    '''
    return m.sin(obliquity) * m.sin(true_longitude)    

def hour_angle(T):
    '''
    Hour angle
    Appelbaum & Flood equation (8)
    Parameters:
         T     Time in Planetary hours
    '''
    return m.radians(15*T-180)

def cos_zenith_angle(obliquity,true_longitude,latitude,T):
    '''
    Cosine of zenith angle
    Appelbaum & Flood equation (6)
    Goosse et al equation (2.21)
    See also Derivation of the solar geometric 
    relationships using vector analysis by Alistair Sproul

    Renewable Energy 32 (2007) 1187-1205
    '''
    sin_decl=sin_declination(obliquity,true_longitude)
    #print ('D',m.degrees(m.asin(sin_decl)))
    cos_decl=m.sqrt(1-sin_decl*sin_decl)
    return m.sin(latitude)*sin_decl + m.cos(latitude)*cos_decl *  m.cos(hour_angle(T))


    
class Solar:
    '''
    Model solar irradiation.
    
    Attributes:
       planet  The planet that we are processing
       S       Solar constant at the mean Sun-Earth distance of l AU, in N/m2
    '''
    def __init__(self,planet,S = 1371):
        self.planet=planet
        self.S = S   # Solar constant at the mean Sun-Earth distance of l AU, in N/m2
                     # Appelbaum & Flood        
        

    def beam_irradience(self,r):
        '''
           Beam Irradience in W/m2
           Appelbaum & Flood equation (1)
           Goosse et al equation (2.18)
        '''        
        return self.S/(r*r)
 
    def surface_irradience(self,true_longitude,latitude,T):
        '''
        Beam irradience on a horizonal surface
           Appelbaum & Flood equations (5) & (6)
           Goosse et al equation (2.21)
        '''
        return max(0,
                   cos_zenith_angle(
                       self.planet.obliquity,
                       true_longitude,
                       latitude,T) *       \
                   self.beam_irradience(
                       self.planet.instantaneous_distance(true_longitude)))
    
    def surface_irradience_daily(self,true_longitude,latitude):
        '''
        Daily Insolation on a horizontal surface
        Goosse et al equation (2.26)
        '''
        ha = self.hour_angle_sunrise_sunset(true_longitude,latitude)
        sin_decl = sin_declination(self.planet.obliquity,true_longitude)
        cos_decl = m.sqrt(1-sin_decl*sin_decl)
        beam_irradience = self.beam_irradience(self.planet.instantaneous_distance(true_longitude))
        return beam_irradience * (86400/m.pi) * (ha*m.sin(latitude)*sin_decl + m.cos(latitude)*cos_decl*m.sin(ha))
    
    def hour_angle_sunrise_sunset(self,true_longitude,latitude,sunset=True):
        '''
        Hour angle for Sunrise & Sunset
        Goosse et al equation (2.24)
        '''
        sin_decl = sin_declination(self.planet.obliquity,true_longitude)
        tan_declination = sin_decl/m.sqrt(1-sin_decl*sin_decl)
        prod=tan_declination*m.tan(latitude)
        if abs(prod)<1:
            return  m.acos(-prod)*( 1 if sunset else -1)
        else: # See Goosse et al, pargraph preceding figure 2.11
            if latitude>0:
                if prod>1:
                    if 0<true_longitude and true_longitude<m.pi:
                        return m.pi
            if latitude<0:
                if prod>1:
                    if m.pi<true_longitude and true_longitude<2*m.pi:
                        return m.pi
            return 0
    
    def length_of_day(self,true_longitude,latitude):
        return max(0,
                   (24/m.pi)*self.hour_angle_sunrise_sunset(true_longitude,latitude))

if __name__=='__main__':
    import unittest
    
    obliquity       = m.radians(23.4)
    
    class TestDeclination(unittest.TestCase):
        def test_winter_solstice(self):            
            true_longitude = 3*m.pi/2
            self.assertAlmostEqual(-obliquity,
                                   m.asin(sin_declination(obliquity,true_longitude)),
                                   places=1)
        def test_summer_solstice(self):
            true_longitude = m.pi/2
            self.assertAlmostEqual(obliquity,
                                   m.asin(sin_declination(obliquity,true_longitude)),
                                   places=1) 
        def test_vernal_equinox(self):
            true_longitude = m.pi
            self.assertAlmostEqual(0,
                                   m.asin(sin_declination(obliquity,true_longitude)),
                                   places=1)
        def test_autumn_equinox(self):
            true_longitude = 0
            self.assertAlmostEqual(0,
                                   m.asin(sin_declination(obliquity,true_longitude)),
                                   places=1)        
        
    class TestZenithAngle(unittest.TestCase):
        def test_winter_solstice(self):
            true_longitude = 3*m.pi/2
            latitude        = 0
            self.assertAlmostEqual(obliquity+latitude,
                                   m.acos(cos_zenith_angle(obliquity,true_longitude,latitude,12)),places=1)
            
        def test_summer_solstice(self):
            true_longitude = m.pi/2
            latitude        = 0
            self.assertAlmostEqual(obliquity+latitude,
                                   m.acos(cos_zenith_angle(obliquity,true_longitude,latitude,12)),places=1)            
                                   
        def test_autumn_equinox(self):
            true_longitude = 0
            latitude        = 0
            self.assertAlmostEqual(latitude,
                                   m.acos(cos_zenith_angle(obliquity,true_longitude,latitude,12)),
                                   places=1)
            
        def test_spring_equinox(self):
            true_longitude = m.pi
            latitude        = m.pi/2
            self.assertAlmostEqual(latitude,
                                   m.acos(cos_zenith_angle(obliquity,true_longitude,latitude,12)),
                                   places=1) 
            
        #def test_winter_solstice_sp(self):
            #true_longitude = 3*m.pi/2
            #latitude        = -m.pi/2
            #self.assertAlmostEqual(obliquity+latitude,
                                   #m.acos(cos_zenith_angle(obliquity,true_longitude,latitude,12)),places=1)         
        #def test_poles(self):
            #latitude = -m.pi/2+obliquity
            #true_longitude = -m.pi/2
            #self.assertAlmostEqual(42,
                                   #m.acos(cos_zenith_angle(obliquity,true_longitude,latitude,12)),
                                   #places=1)            
    try:
        unittest.main()
    except SystemExit as inst:
        if inst.args[0]: # raised by sys.exit(True) when tests failed
            raise     