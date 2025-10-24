#%%
import pyomo.environ as pyo
from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
    value,
    Var,
    assert_optimal_termination,
)
from pyomo.network import Arc
import watertap.unit_models.zero_order as unit_ZO
from pyomo.util.check_units import assert_units_consistent
from watertap.unit_models.zero_order import MBBRZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.costing.zero_order_costing import ZeroOrderCosting
from pyomo.util.check_units import assert_units_consistent
from idaes.models.unit_models import Feed, Product
from idaes.core import FlowsheetBlock
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core import UnitModelCostingBlock

#%%
m = ConcreteModel()
m.db = Database()

m.fs = FlowsheetBlock(dynamic=False)
m.fs.prop = WaterParameterBlock(
        water_source="wastewater",
        solute_list=[
            #"tss",
            #"tds",
            #"nitrogen",
            #"phosphates",
            #"phosphorus",
            "bod",
            "nitrate",
            "ammonium_as_nitrogen",
        ],
    )

    # unit models
m.fs.feed = unit_ZO.FeedZO(property_package=m.fs.prop)
m.fs.unit = MBBRZO(property_package=m.fs.prop, database=m.db, process_subtype="bod-stage")
m.fs.unit.load_parameters_from_database(use_default_removal=True)

m.fs.feed.flow_vol[0].fix(30000 * pyo.units.m**3 / pyo.units.day)
m.fs.feed.conc_mass_comp[0, "bod"].fix(140 * pyo.units.g / pyo.units.m**3)
m.fs.feed.conc_mass_comp[0, "nitrate"].fix(0)
m.fs.feed.conc_mass_comp[0, "ammonium_as_nitrogen"].fix(100.7 * pyo.units.g / pyo.units.m**3)
#m.fs.unit.volume.fix(1050 * pyo.units.m**3)
m.fs.unit.media_specific_area.fix(500)
m.fs.unit.fill_fraction.fix(0.50)
#m.fs.unit.bod_removal_rate.fix()
#m.fs.unit.reaction_conversion[0, "bod_removal"].fix()
m.fs.feed_to_mbbr = Arc(source=m.fs.feed.outlet, destination = m.fs.unit.inlet)
pyo.TransformationFactory("network.expand_arcs").apply_to(m)
#%%
m.fs.unit.extent_of_reaction.pprint()

# %%
m.fs.unit.report()
#%%
degrees_of_freedom(m)
#%%
solver = get_solver()
results = solver.solve(m)
# %%
m.fs.unit.report()

# %%
