# Copyright (C) 2019 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Functions for validating Maud objects."""

import warnings

from maud import data_model


def validate_maud_input(mi: data_model.MaudInput):
    """Check that priors, experiments and kinetic model are consistent."""

    def missing_coords_message(level, missing_coords):
        return f"{level} labels in coords but not priors: {missing_coords}"

    def missing_priors_message(level, param, missing):
        return f"{level} labels in {param} priors but not coords: {missing}"

    ps = mi.priors
    cs = mi.stan_coords
    should_match_1d = [
        (ps.priors_km, "mic_id", cs.km_mics),
        (ps.priors_km, "enzyme_id", cs.km_enzs),
        (ps.priors_kcat, "enzyme_id", cs.enzymes),
    ]
    should_match_2d = [
        (
            ps.priors_drain,
            ["experiment_id", "drain_id"],
            [cs.experiments, cs.drains],
        ),
        (
            ps.priors_conc_unbalanced,
            ["experiment_id", "mic_id"],
            [cs.experiments, cs.unbalanced_mics],
        ),
    ]
    for p, level, c in should_match_1d:
        param = p.parameter_name
        if len(p.location) == 0:
            assert len(c) == 0, missing_coords_message(level, c)
            continue
        if len(c) == 0:
            ix_ps = p.location.index
            assert len(ix_ps) == 0, missing_priors_message(level, param, ix_ps)
            continue
        ix_ps = set(p.location.index.get_level_values(level))
        ix_cs = set(c)
        assert ix_ps.issubset(ix_cs), missing_coords_message(level, ix_ps - ix_cs)
        assert ix_cs.issubset(ix_ps), missing_priors_message(
            level, param, ix_cs - ix_ps
        )
    for p, levels, cs in should_match_2d:
        param = p.parameter_name
        ix_ps_row = set(p.location.index.values)
        ix_ps_col = set(p.location.columns.values)
        ix_cs_row = set(cs[0])
        ix_cs_col = set(cs[1])
        if not ix_ps_row.issubset(ix_cs_row):
            msg = missing_coords_message(levels[0], ix_ps_row - ix_cs_row)
            warnings.warn(msg)
        if not ix_ps_col.issubset(ix_cs_col):
            msg = missing_coords_message(levels[1], ix_ps_col - ix_cs_col)
            warnings.warn(msg)
        if not ix_cs_col.issubset(ix_ps_col):
            msg = missing_priors_message(levels[1], param, ix_cs_col - ix_ps_col)
            warnings.warn(msg)
