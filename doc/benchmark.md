\page benchmark Benchmarking results for Vagabond

This page acts more like a lab notebook than a documentation page.

## Vagabond v0.0.0.589bcee

This tests the helpfulness of a default bond-phi angle of 


## Vagabond v0.0.0.eae5928

Reverted the 3-monomer to 2-monomer fit, back to 3 monomers compared to
v0.0.0.92138c5.

Average Rwork/Rfree: 28.04 / 30.24% (**update me**)

## Vagabond v0.0.0.92138c5

Reverted the refinement of libration compared to 0.0.0.06408992

Average Rwork/Rfree: 28.18 / 30.39% (**update me**)

## Vagabond v0.0.0.06408992

Added 'undo if worse' which is applied to the global fit stage. This will
force crystals to correct with intramolecular motion if intermolecular motion
makes the Rwork worse.

Directory is called: 06408992
Average Rwork/Rfree: 28.20 / 30.40% (**update me**)

## Vagabond v0.0.0.4b44bb9c

Changed 3 monomer fit to 2 monomer fit, refining libration against B factor
targets from electron density and minor modifications to refinement in FlexLocal.

Average Rwork/Rfree: 28.41 / 30.63% (**update me**)


## Vagabond v0.0.0.95d54531

Benchmarking Vagabond after B factor fit of 3 monomers at a time to real electron
density results. This is not as good as the previous benchmark, but we are now
fitting directly to the electron density, which is harder.

Directory is called: reflex3

Average Rwork/Rfree: 28.30 / 30.45%

## Vagabond v0.0.0.25ac0445

Fitting against the original positions and B factor from Refmac. Just after
the removal of the non-cube FFT grid bug.

Directory is called: antgone

Average Rwork/Rfree: 28.04 / 30.24%
