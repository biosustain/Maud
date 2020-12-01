"""Integration test from code generation to correct sampling."""

import os
import shutil
import tempfile

import maud.sampling as sampling


here = os.path.dirname(__file__)
data_path = os.path.join(here, '..','data')


def test_linear():
    """Tests linear model.

    tests from code generation to sampling of the linear model by computing 200
    samples after 200 warmups and checking if the sampled median is within the
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
        "formation_energy_z[1]": [-1.6088815, 1.7438237],
        "formation_energy_z[2]": [-1.9428797, 1.8600326],
        "kcat[1]": [0.7354878399999999, 3.3920119],
        "kcat[2]": [0.7169074799999999, 3.5383422],
        "kcat[3]": [0.5103703599999999, 2.3251939000000004],
        "km[1]": [0.29829006999999996, 2.1484104000000004],
        "km[2]": [0.39007882, 3.3307122000000002],
        "km[3]": [0.30081803999999995, 2.1872641],
        "km[4]": [0.30032217, 2.956041],
        "km[5]": [0.35870995999999994, 2.4696105],
        "km[6]": [0.35748458, 3.1736108000000005],
        "enzyme[1,1]": [0.90327701, 1.1091029000000001],
        "enzyme[2,1]": [0.45675294, 0.55617683],
        "enzyme[1,2]": [0.9169324400000001, 1.0900779],
        "enzyme[2,2]": [1.3650513999999998, 1.6433149],
        "enzyme[1,3]": [0.9258319500000001, 1.087983],
        "enzyme[2,3]": [1.3852606, 1.6246428],
        "conc_unbalanced[1,1]": [1.8437724, 2.2289229],
        "conc_unbalanced[2,1]": [1.83488, 2.1967661],
        "conc_unbalanced[1,2]": [0.92532526, 1.0731987],
        "conc_unbalanced[2,2]": [0.9100531199999999, 1.0846023999999999],
        "ki[1]": [0.35249536, 2.0479711000000003],
        "dissociation_constant_t[1]": [0.39240253999999997, 2.9900479000000004],
        "dissociation_constant_r[1]": [0.3673301, 2.4169647000000003],
        "transfer_constant[1]": [0.26424406999999994, 2.6506271],
        "transfer_constant[2]": [0.27446419000000005, 2.8339571],
        "formation_energy[1]": [-1.0804408, -0.9128092999999999],
        "formation_energy[2]": [-2.0971405, -1.9069949],
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
        "keq[1]": [1.4224362, 1.5845171],
        "keq[2]": [1.0, 1.0],
        "keq[3]": [1.0, 1.0],
        "log_lik[1]": [-1.2944431, 1.3392327],
        "log_lik[2]": [0.014465301999999961, 1.376049],
        "log_lik[3]": [-0.73380201, 0.85230665],
        "log_lik[4]": [-0.55109339, 1.0460312],
        "log_lik[5]": [-0.9655723700000002, 1.3825698],
        "log_lik[6]": [0.5349423899999999, 2.07629],
        "log_lik[7]": [-0.014904273000000003, 0.7953745],
        "log_lik[8]": [-0.50145479, 1.1199048],
        "log_lik[9]": [-0.87171379, 1.3808676000000002],
        "log_lik[10]": [-0.07020213300000022, 2.0754574999999997],
    }

    linear_input = os.path.join(data_path, "linear.toml")
    temp_directory = tempfile.mkdtemp(dir=data_path)
    linear_input_values = {
        "abs_tol": 1e-6,
        "rel_tol": 1e-6,
        "max_num_steps": int(1e9),
        "likelihood": 1,
        "n_samples": 50,
        "n_warmup": 50,
        "n_chains": 1,
        "n_cores": 1,
        "timepoint": 500,
        "output_dir": temp_directory,
        "data_path": linear_input,
        "threads_per_chain": 1,
        "save_warmup": False,
    }
    fit = sampling.sample(**linear_input_values)
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
