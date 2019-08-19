import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import os


def plot_pairs(infd, var_names, coords=None, log_transform=True, png_location=None):
    plt.close('all')
    axes = az.plot_pair(infd, var_names, coords, figsize=(25, 20), divergences=True)
    if log_transform:
        for ax in axes.ravel():
            if ax.get_xticks()[1] >= 0:
                ax.semilogx()
            if ax.get_yticks()[1] >= 0:
                ax.semilogy()
    if png_location is not None:
        plt.savefig(png_location, format='png')
    return axes


def plot_posterior(infd, var_names, coords=None, log_transform=True, png_location=None):
    plt.close('all')
    axes = az.plot_posterior(infd, var_names, coords, figsize=(25, 20), kind='histogram')
    if log_transform:
        for ax in axes.ravel():
            if ax.get_xticks()[1] >= 0:
                ax.semilogx()
    if png_location is not None:
        plt.savefig(png_location, format='png')
    return axes


def plot_pairs_for_reaction_kinetics(infd, reaction, figure_path):
    var_names = ['kinetic_parameter']
    all_coords = infd.posterior.coords
    kp_names = list(filter(lambda i: reaction in i, all_coords['parameter_names'].to_series()))
    coords = {'parameter_names': kp_names}
    axes = plot_pairs(infd, var_names, coords)
    out_path = os.path.join(figure_path, f'{var_names[0]}_correlation_{reaction}.png')
    axes = plot_pairs(infd, var_names, coords, png_location=out_path)
    return axes


def plot_correlation_for_fluxs(infd, experiment, figure_path):
    var_names = ['flux']
    all_coords = infd.posterior.coords
    kp_names = list(filter(lambda i: experiment in i, all_coords['experiments'].to_series()))
    coords = {'experiments': kp_names}
    axes = plot_pairs(infd, var_names, coords)
    out_path = os.path.join(figure_path, f'{var_names[0]}_corellation_{experiment}.png')
    axes = plot_pairs(infd, var_names, coords, log_transform=False, png_location=out_path)
    return axes

def plot_distribution_for_fluxs(infd, experiment, figure_path):
    var_names = ['flux']
    all_coords = infd.posterior.coords
    kp_names = list(filter(lambda i: experiment in i, all_coords['experiments'].to_series()))
    coords = {'experiments': kp_names}
    axes = plot_pairs(infd, var_names, coords)
    out_path = os.path.join(figure_path, f'{var_names[0]}_distribution_{experiment}.png')
    axes = plot_posterior(infd, var_names, coords, log_transform=False, png_location=out_path)
    return axes


def plot_pairs_condition_concentrations(infd, experiment, figure_path, type):
    if type == 'inference':
        var_names = ['concentration']
    elif type == 'relative_inference':
        var_names = ['scaled_concentration']
    all_coords = infd.posterior.coords
    experiment_names = list(filter(lambda i: experiment in i, all_coords['experiments'].to_series()))
    coords = {'experiments': experiment_names}
    out_path = os.path.join(figure_path, f'{var_names[0]}_correlation_{experiment}.png')
    axes = plot_pairs(infd, var_names, coords, png_location=out_path)
    return axes


def plot_distributions_condition_concentrations(infd, experiment, figure_path, type):
    plt.close('all')
    n_rows = round(len(infd.posterior.metabolites)/3+0.5)
    fig, axes = plt.subplots(n_rows, 3,
                             tight_layout=True,
                             figsize=(25, (5*n_rows)))

    metabolite_list = [x for x in infd.posterior.metabolites.to_series()]

    min_bin_list = []
    max_bin_list = []

    for met in metabolite_list:
        if type == 'inference':
            min_bin_list.append(np.min(infd.posterior.concentration.loc[{
                                                    'metabolites': [met],
                                                    'experiments': [experiment]
                                                    }]))
            max_bin_list.append(np.max(infd.posterior.concentration.loc[{
                                                    'metabolites': [met],
                                                    'experiments': [experiment]
                                                    }]))
        elif type == 'relative_inference':
            min_bin_list.append(np.min(infd.posterior.scaled_concentration.loc[{
                                                    'metabolites': [met],
                                                    'experiments': [experiment]
                                                    }]))
            max_bin_list.append(np.max(infd.posterior.scaled_concentration.loc[{
                                                    'metabolites': [met],
                                                    'experiments': [experiment]
                                                    }]))

    bin_list = [np.logspace(np.log10(x), np.log10(y), 20)
                for x, y in zip(min_bin_list,
                                max_bin_list)]

    for ax, met, bin in zip(axes.ravel(),
                            metabolite_list,
                            bin_list):
        if type == 'inference':
            val = infd.posterior.concentration.loc[{'experiments': experiment,
                                                    'metabolites': met}]
        elif type == 'relative_inference':
            val = infd.posterior.scaled_concentration.loc[
                                                {'experiments': experiment,
                                                 'metabolites': met}]

        ax.hist(val.to_series().ravel(), bins=bin)
        if ax.get_xticks()[1] >= 0:
            ax.semilogx()
        ax.set_title(f'{experiment} {met}')

    for ax in axes.ravel()[len(infd.posterior.metabolites):]:
        ax.axis('off')

    out_path = os.path.join(figure_path, f'{experiment}_posterior_metabolite_distributions.png')
    fig.savefig(out_path, format='png')

    return axes


def plot_distributions(model_name: str, model_type: str='inference'):
    inference_type_template = {
        'relative_inference': f'relative_inference_{model_name}.nc',
        'inference': f'model_inference_{model_name}.nc'
    }

    data_folder = f'../data/out'
    here = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(here, data_folder)

    paths = {
        'input_file': os.path.join(data_path,
                                    inference_type_template[model_type]),
        'output_figure_path': os.path.join(data_path,
                                        f'{model_name}_{model_type}_plots')
    }

    if not os.path.exists(paths['output_figure_path']):
        os.mkdir(paths['output_figure_path'])

    if not os.path.exists(data_path):
        raise Exception(f'{data_path} does not exist, please check name or run simulation to plot')

    data = az.from_netcdf(paths['input_file'])

    experiment_list = [x for x in data.posterior.experiments.to_series()]
    reaction_list = [x for x in data.posterior.reactions.to_series()]

    # Plotting concentration distributions and inter-condition correlations
    for cond in experiment_list:
        plt.figure()
        plot_pairs_condition_concentrations(data,
                                            cond,
                                            paths['output_figure_path'],
                                            model_type)

        plt.figure()
        plot_distributions_condition_concentrations(data,
                                                    cond,
                                                    paths['output_figure_path'],
                                                    model_type)

        plt.figure()
        plot_correlation_for_fluxs(data,
                                    cond,
                                    paths['output_figure_path'])

        plt.figure()
        plot_distribution_for_fluxs(data,
                                    cond,
                                    paths['output_figure_path'])

    # Plotting kinetic parameter distributions and inter-reaction correlations
    for rxn in reaction_list:
        plt.figure()
        plot_pairs_for_reaction_kinetics(data,
                                         rxn,
                                         paths['output_figure_path'])
