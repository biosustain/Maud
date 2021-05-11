"""Integration test from code generation to correct sampling."""

import os
import shutil
import tempfile

from maud.io import load_maud_input_from_toml
from maud.sampling import sample


here = os.path.dirname(__file__)
data_path = os.path.join(here, "../data")


def test_linear():
    """Tests linear model.

    tests from code generation to sampling of the linear model by computing 50
    samples after 50 warmups and checking if the sampled median is within the
    94% CI of a precomputed set.
    """
    expected = {
        # marginal 3% and 97% intervals from a control run. Use the following
        # code to generate a dictionary like this from a CmdStanMCMC object
        # called 'fit':
        #
        # draws = fit.draws_pd()
        # expected = dict(zip(
        #     draws.columns, map(list, draws.quantile([0.03, 0.97]).T.values)
        # ))
        "lp__": [-29.316072, -16.076300999999997],
        "accept_stat__": [0.96072929, 0.99988951],
        "stepsize__": [0.133839, 0.133839],
        "treedepth__": [5.0, 5.0],
        "n_leapfrog__": [31.0, 31.0],
        "divergent__": [0.0, 0.0],
        "energy__": [26.507623, 44.448424],
        "dgf_z[1]": [-1.6088815, 1.7438237],
        "dgf_z[2]": [-1.9428797, 1.8600326],
        "kcat[1]": [0.7354878399999999, 3.3920119],
        "kcat[2]": [0.7169074799999999, 3.5383422],
        "kcat[3]": [0.5103703599999999, 2.3251939000000004],
        "km[1]": [0.29829006999999996, 2.1484104000000004],
        "km[2]": [0.39007882, 3.3307122000000002],
        "km[3]": [0.30081803999999995, 2.1872641],
        "km[4]": [0.30032217, 2.956041],
        "km[5]": [0.35870995999999994, 2.4696105],
        "km[6]": [0.35748458, 3.1736108000000005],
        "conc_enzyme[1,1]": [0.90327701, 1.1091029000000001],
        "conc_enzyme[2,1]": [1.3650513999999998, 1.6433149],
        "conc_enzyme[1,2]": [0.9169324400000001, 1.0900779],
        "conc_enzyme[2,2]": [0.45675294, 0.55617683],
        "conc_enzyme[1,3]": [0.9258319500000001, 1.087983],
        "conc_enzyme[2,3]": [1.3852606, 1.6246428],
        "conc_unbalanced[1,1]": [1.8437724, 2.2289229],
        "conc_unbalanced[2,1]": [1.83488, 2.1967661],
        "conc_unbalanced[1,2]": [0.92532526, 1.0731987],
        "conc_unbalanced[2,2]": [0.9100531199999999, 1.0846023999999999],
        "ki[1]": [0.35249536, 2.0479711000000003],
        "diss_t[1]": [0.39240253999999997, 2.9900479000000004],
        "diss_r[1]": [0.3673301, 2.4169647000000003],
        "transfer_constant[1]": [0.26424406999999994, 2.6506271],
        "transfer_constant[2]": [0.27446419000000005, 2.8339571],
        "dgf[1]": [-1.0804408, -0.9128092999999999],
        "dgf[2]": [-2.0971405, -1.9069949],
        "conc[1,1]": [1.8437724, 2.2289229],
        "conc[2,1]": [1.83488, 2.1967661],
        "conc[1,2]": [0.92532526, 1.0731987],
        "conc[2,2]": [0.9100531199999999, 1.0846023999999999],
        "conc[1,3]": [1.4225884000000002, 1.8135276],
        "conc[2,3]": [1.5866924, 1.9812213],
        "conc[1,4]": [1.2548785, 1.6608011999999999],
        "conc[2,4]": [1.0856881999999999, 1.3590284],
        "flux[1,1]": [0.058571784999999994, 0.26232753],
        "flux[2,1]": [0.044527502, 0.2225318],
        "flux[1,2]": [0.058571784999999994, 0.26232753],
        "flux[2,2]": [0.044527502, 0.2225318],
        "flux[1,3]": [0.058571781999999996, 0.26232753],
        "flux[2,3]": [0.044527502, 0.2225318],
        "keq[1]": [1.0, 1.0],
        "keq[2]": [1.4224362, 1.5845171],
        "keq[3]": [1.0, 1.0],
        "log_lik_flux[1]": [-1.2944431, 1.3392327],
        "log_lik_flux[2]": [0.014465301999999961, 1.376049],
        "log_lik_conc[1]": [-0.73380201, 0.85230665],
        "log_lik_conc[2]": [-0.55109339, 1.0460312],
        "log_lik_conc[3]": [-0.9655723700000002, 1.3825698],
        "log_lik_conc[4]": [0.5349423899999999, 2.07629],
        "log_lik_conc[5]": [-0.014904273000000003, 0.7953745],
        "log_lik_conc[6]": [-0.50145479, 1.1199048],
        "log_lik_conc[7]": [-0.87171379, 1.3808676000000002],
        "log_lik_conc[8]": [-0.07020213300000022, 2.0754574999999997],
    }

    linear_input = os.path.join(data_path, "linear")
    temp_directory = tempfile.mkdtemp(dir=data_path)
    mi = load_maud_input_from_toml(linear_input)
    mi.config.cmdstanpy_config.update(
        {
            "chains": 1,
            "iter_sampling": 50,
            "iter_warmup": 50,
            "save_warmup": False,
        }
    )
    fit = sample(mi, output_dir=temp_directory)
    samples_test = fit.draws_pd()
    # Check that each output column (other than the diagnostic ones) is
    # statistically similar to its matching control column.
    test_mean = samples_test.mean()
    cols = [c for c in expected.keys() if not c.endswith("__")]
    for col in cols:
        assert col in samples_test.columns, col + " is not present in test"
        assert test_mean[col] >= expected[col][0], col + " is too low."
        assert test_mean[col] <= expected[col][1], col + " is too high."
    # Delete temporary directory
    shutil.rmtree(temp_directory)
