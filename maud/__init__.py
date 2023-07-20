import os

import cmdstanpy

from maud.running_stan import (
    CMDSTAN_VERSION,
    STAN_FILES_FOLDER,
    load_stan_model,
)

# on Windows specifically, we should point cmdstanpy to the repackaged
# CmdStan if it exists. This lets cmdstanpy handle the TBB path for us.
local_cmdstan = STAN_FILES_FOLDER / f"cmdstan-{CMDSTAN_VERSION}"
if os.path.exists(local_cmdstan):
    cmdstanpy.set_cmdstan_path(str((local_cmdstan.resolve())))

MODEL = load_stan_model("model", cpp_options={}, stanc_options={})
OUT_OF_SAMPLE_MODEL = load_stan_model(
    "out_of_sample_model", cpp_options={}, stanc_options={}
)
