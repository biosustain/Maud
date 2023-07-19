from maud.running_stan import load_stan_model

MODEL = load_stan_model("model", cpp_options={}, stanc_options={})
OUT_OF_SAMPLE_MODEL = load_stan_model(
    "out_of_sample_model", cpp_options={}, stanc_options={}
)
