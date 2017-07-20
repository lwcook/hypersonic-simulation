import numpy as np
import math

import geometry

class WCV(object):

    def __init__(self):

        self.geometry = geometry.Geometry(fidelity=13)

        self.exposed_planform_area = 0
        self.internal_volume = 0
        self.inlet_capture_width = 0
        self.fuel_volume = 0
        self.mass = calculateMassFromScaling()



    def calculateMassFromScaling(self):

    # Scale the vehicle mass based on total surface area of vehicle

        mass_model_type = 1;

    #         wing_area_base = 58.3872;
    #         wing_area_ratio = self.wing_area/wing_area_base;
    #         surface_area_base = 93.887;
    #         surface_area_ratio = self.surface_area/surface_area_base;

        exposed_area_base = 20.1272;
        exposed_area_ratio = self.exposed_planform_area/exposed_area_base;

        internal_volume_base = 33.5431;
        internal_volume_ratio = self.internal_volume/internal_volume_base;

        scramjet_width_base = 0.64;
        scramjet_ratio = self.inlet_capture_width / scramjet_width_base;

        fuel_volume_base = 12.5952;
        fuel_volume_ratio = self.fuel_volume / fuel_volume_base;

        # Empty mass in kg
        m_fuselage_base = 2019;
        m_wings_base = 209;
        m_scramjets_base = 473;

        m_tank = 135;
        m_systems = 707;
        m_landing_gear = 142;

        m_fuselage = m_fuselage_base * internal_volume_ratio;
        m_wings = m_wings_base * exposed_area_ratio;
        m_scramjets = m_scramjets_base * scramjet_ratio;
        m_other = (m_tank*fuel_volume_ratio + m_systems + m_landing_gear);

    #         m_body_base = m_fuselage + m_tank + m_systems + m_landing_gear + m_scramjets; # reusable mass
    #         m_fuselage = m_body_base*surface_area_ratio;
    #         m_wings = m_wings_base*surface_area_ratio;

            m_vehicle = m_scramjets + m_fuselage + m_wings + m_other ; # Reusable vehicle mass

            # Payload and fuel mass
            m_payload = 1000;

            fuel_density = 71; # 71 kg/m^3 for hydrogen
            m_fuel = self.fuel_volume*fuel_density;

    #         fuselage_specific_mass = 61.5; # kg/m^3
    #         wing_specific_mass = 148.3; # kg/m^3 
    #         fuselage_mass = fuselage_specific_mass*self.internal_volume;
    #         wing_mass = wing_specific_mass*self.planform_area;
    #         scramjet_mass = 0.25*structural_mass;

        return m_vehicle, m_fuel, m_payload
