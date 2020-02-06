# 10. Use only modular rate laws

Date: 2020-02-05

## Status

Accepted

## Context

The `code_generation.py` module is probably the ugliest part of Maud,
particularly the function `create_fluxes_function`. It could be cleaned up a
bit by removing non-modular enzyme mechanisms.

Non-modular mechanisms aren't used in any of the models we are currently
testing, but for historical reasons the logic in `code_generation.py` basically
assumes that they are the default. The modular rate law is kind of an
add-on. This change would let us move towards switching this around, having the
modular rate law as the default and working in other mechanisms as and when
necessary.

## Decision

We remove non-modular rate laws from the `code_generation.py` module.

## Consequences

We won't be able to use non-modular rate laws without writing them back in.
