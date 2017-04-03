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
Model for solar irradiation, based on Solar Radiation on Mars, 
 Joseph Appelbaum & Dennis Flood, Lewis Research Center, NASA 
'''


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
        
#   Beam Irradience in W/m2
#   Appelbaum & Flood equation (1)
    def beam_irradience(self,r):
        return self.S/(r*r)
 
#   Beam irradience on a horizonal surface
#   Appelbaum & Flood equations (5) & (6)
    def surface_irradience(self,true_longitude,latitude,T):
        cos_zenith_angle = self.planet.cos_zenith_angle(true_longitude,latitude,T)
        r=self.planet.instantaneous_distance(true_longitude)
        beam_irradience = self.beam_irradience(r)
        return max(0,cos_zenith_angle*beam_irradience)
    
