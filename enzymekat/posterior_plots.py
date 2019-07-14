import arviz as az
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

os.chdir('/Users/nicow/Documents/enzymeKAT/enzymekat')
data = az.from_netcdf('../data/out/infd_relative_metabolomics.nc')
data_tr= az.from_netcdf('../data/out/infd_training.nc')

data_sum = az.summary(data)

def plot_pairs(infd, var_names, coords=None, log_transform=True, png_location=None):
    plt.close('all')
    axes = az.plot_pair(infd, var_names, coords, figsize=(25, 20))
    if log_transform:
        for ax in axes.ravel():
            if ax.get_xticks()[1] >= 0:
                ax.semilogx()
            if ax.get_yticks()[1] >= 0:
                ax.semilogy()
    if png_location is not None:
        plt.savefig(png_location)
    return axes

def plot_pairs_for_reaction(infd, reaction, png_location=None):
    var_names = ['kinetic_parameters']
    all_coords = infd.posterior.coords
    kp_names = list(filter(lambda i: reaction in i, all_coords['kinetic_parameter_names'].to_series()))
    coords = {'kinetic_parameter_names': kp_names}
    axes = plot_pairs(infd, var_names, coords, log_transform=True, png_location=png_location)
    return axes

plt.close('all')

axes = az.plot_pair(data, var_names=['metabolite_scaling'])
plt.savefig('metabolite_scaling.png')



    met_list = ['metabolite' in x for x in data_tr_sum.index]
flux_list = ['flux' in x for x in data_tr_sum.index]


plot_pairs_for_reaction(data, 'RXN4', 'pair_plot_RXN4');


    infd = az.from_cmdstan(
        "../data/out/model_output_relative_metabolomics.csv",
        coords={
            'kinetic_parameter_names': data_tr.posterior.kinetic_parameter_names,
            'measurement_names': data_tr.posterior.measurement_names},
        dims={
            'kinetic_parameters': ['kinetic_parameter_names'],
            'metabolite_concentration_hat': ['measurement_names']
        })
    infd.to_netcdf("../data/out/infd_relative_metabolomics.nc")


# <md>
# converting the chains into a single array

# <codecell>

metabolite_posterior_draws = np.stack(data.posterior['metabolite_concentration_hat'][0])
flux_posterior_draws = np.stack(data.posterior['flux_hat'][0])
kinetic_posterior_draws = np.stack(data.posterior['kinetic_parameters'][0])

metabolite_posterior_draws_tr = np.stack(data_tr.posterior['metabolite_concentration_hat'][0])
flux_posterior_draws_tr = np.stack(data_tr.posterior['flux_hat'][0])
kinetic_posterior_draws_tr = np.stack(data_tr.posterior['kinetic_parameters'][0])


for i in np.arange(1, 4):
    kinetic_posterior_draws = np.append(kinetic_posterior_draws, np.stack(np.array(data.posterior['kinetic_parameters'][i])), axis=0)
    metabolite_posterior_draws = np.append(metabolite_posterior_draws, np.stack(np.array(data.posterior['metabolite_concentration_hat'][i])), axis=0)
    flux_posterior_draws = np.append(flux_posterior_draws, np.stack(np.array(data.posterior['flux_hat'][i])), axis=0)
    kinetic_posterior_draws_tr = np.append(kinetic_posterior_draws_tr, np.stack(np.array(data_tr.posterior['kinetic_parameters'][i])), axis=0)
    metabolite_posterior_draws_tr = np.append(metabolite_posterior_draws_tr, np.stack(np.array(data_tr.posterior['metabolite_concentration_hat'][i])), axis=0)
    flux_posterior_draws_tr = np.append(flux_posterior_draws_tr, np.stack(np.array(data_tr.posterior['flux_hat'][i])), axis=0)

# <md>
# preparing the data for the sampled data

# <codecell>

kinetic_posterior_df = pd.DataFrame(kinetic_posterior_draws)
kinetic_posterior_df = np.log(kinetic_posterior_df)

metabolite_posterior_df = pd.DataFrame(metabolite_posterior_draws)

flux_posterior_df = pd.DataFrame(flux_posterior_draws)

# <md>
# preparing the array for the training data

# <codecell>

kinetic_posterior_tr_df = pd.DataFrame(kinetic_posterior_draws_tr)
kinetic_posterior_tr_df = np.log(kinetic_posterior_tr_df)

metabolite_posterior_tr_df = pd.DataFrame(metabolite_posterior_draws_tr)

flux_posterior_tr_df = pd.DataFrame(flux_posterior_draws_tr)


# <md>
# Plotting the kinetic posterior distribution with the true value distribution additionally plotted

# <codecell>
fig, axes = plt.subplots(nrows=4, ncols=7, figsize=(15, 20));

for ax, dat, training_dat, names in zip(axes.flatten(), kinetic_posterior_df.values.T, kinetic_posterior_tr_df.values.T, data.posterior.kinetic_parameter_names.values):
    bins = np.linspace(min(np.concatenate((dat,training_dat))), max(np.concatenate((dat,training_dat))), 30)
    ax.hist(dat, alpha=0.5, bins=bins)
    ax.hist(training_dat, alpha=0.5, bins=bins)
    ax.set_title(names)

fig.tight_layout()

fig.savefig('hist_kinetic_parameters_rel_big.png')


# <codecell>
fig1, axes1 = plt.subplots(nrows=1, ncols=5, figsize=(15, 5));

for ax, dat, training_dat, names in zip(axes1.flatten(), metabolite_posterior_df.values.T, metabolite_posterior_tr_df.values.T, data.posterior.measurement_names.values):
    bins = np.linspace(min(np.concatenate((dat,training_dat))), max(np.concatenate((dat,training_dat))), 30)
    ax.hist(dat, alpha=0.5, bins=bins)
    ax.hist(training_dat, alpha=0.5, bins=bins)
    ax.set_title(names)
    ax.set_xlabel('Concentration mol/L')

fig1.tight_layout()

fig1.savefig('hist_metabolites_rel_big.png')


# <codecell>
fig2, axes2 = plt.subplots(nrows=1, ncols=6, figsize=(15, 5));

for ax, dat, training_dat, name in zip(axes2.flatten(), flux_posterior_df.values.T, flux_posterior_tr_df.values.T, ['Reaction ' + str(x) for x in range(1,7)]):
    bins = np.linspace(min(np.concatenate((dat,training_dat))), max(np.concatenate((dat,training_dat))), 30)
    ax.hist(dat, alpha=0.5, bins=bins)
    ax.hist(training_dat, alpha=0.5, bins=bins)
    ax.set_title(name)
    ax.set_xlabel('Flux mol/s')

fig2.tight_layout()

fig2.savefig('hist_fluxes_rel_big.png')
