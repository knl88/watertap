#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
This module contains a zero-order representation of a membrane aerated biofilm reactor unit.
"""

import pyomo.environ as pyo
from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class
from watertap.core import build_sido_reactive, ZeroOrderBaseData, constant_intensity

# Some more information about this module
__author__ = "Kim Leirvik"


@declare_process_block_class("MBBRZO")
class MBBRZOData(ZeroOrderBaseData):
    """
    This module contains a zero-order representation of a moving bed biofilm reactor unit.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "mbbr"

        build_sido_reactive(self)
        constant_intensity(self)
        
        self.media_specific_area = Var(
            initialize=500.0,
            units=pyunits.m**2 / pyunits.m**3,
            doc="Specific surface area of media",
        )
        # Different default fill fractions retained
        self._perf_var_dict["Media specific area"] = self.media_specific_area

        self.fill_fraction = Var(
            initialize=0.50,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Volumetric media fill fraction",
        )
        self._perf_var_dict["Fill fraction"] = self.fill_fraction
        # unit variables
        self.volume = Var(
            initialize=1, bounds=(0, None), units=pyunits.m**3, doc="Reactor volume"
        )
        self.hydraulic_retention_time = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.hr,
            doc="Hydraulic residence time",
        )

        self._perf_var_dict["Hydraulic retention time"] = self.hydraulic_retention_time

        @self.Constraint(
            doc="Constraint for reactor volume based on hydraulic residence time"
        )
        def eq_reactor_volume(b):
            return b.volume == (
                pyunits.convert(
                    b.get_inlet_flow(0), to_units=pyunits.m**3 / pyunits.hour
                )
                * b.hydraulic_retention_time
            )

        self._perf_var_dict["Reactor volume"] = self.volume

        self.media_area = Var(
            initialize=250000,
            units=pyunits.m**2,
            bounds=(0, None),
            doc="Sizing variable for effective reactor area",
        )
        if self.config.process_subtype.startswith("bod-stage"):
            self.bod_removal_rate = Var(
                units=pyunits.g / pyunits.m**2 / pyunits.day,
                bounds=(0, None),
                doc="BOD removal rate per day",
            )

            @self.Constraint(
                self.flowsheet().time,
                doc="Constraint for effective reactor area",
            )
            def media_area_constraint(b, t):
                return b.media_area == pyunits.convert(
                    b.properties_in[t].flow_mass_comp["bod"]
                    / (b.bod_removal_rate / b.reaction_conversion[0, "bod_removal"]),
                    to_units=pyunits.m**2,
                )

            self._perf_var_dict["BOD Removal Rate"] = self.bod_removal_rate

            self._fixed_perf_vars.append(self.bod_removal_rate)
        elif self.config.process_subtype.startswith("nitrification-stage"):
            self.nitrogen_removal_rate = Var(
                units=pyunits.g / pyunits.m**2 / pyunits.day,
                bounds=(0, None),
                doc="Nitrogen removal rate per day",
            )

            @self.Constraint(
                self.flowsheet().time,
                doc="Constraint for effective reactor area",
            )
            def media_area_constraint(b, t):
                return b.media_area == pyunits.convert(
                    b.properties_in[t].flow_mass_comp["ammonium_as_nitrogen"]
                    / (
                        b.nitrogen_removal_rate
                        / b.reaction_conversion[0, "nitrification"]
                    ),
                    to_units=pyunits.m**2,
                )

            self._perf_var_dict["Nitrogen Removal Rate"] = self.nitrogen_removal_rate

            self._fixed_perf_vars.append(self.nitrogen_removal_rate)

        self._perf_var_dict["Media Area"] = self.media_area

        @self.Constraint(
            doc="Constraint for reactor volume based on media area and fill fraction"
        )
        def volume_constraint(b):
            return b.volume == pyunits.convert(
                b.media_area / (b.media_specific_area * b.fill_fraction),
                to_units=pyunits.m**3,
            )
