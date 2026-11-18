 PK Trial Simulator
wh
## The model
The body is three connected compartments: plasma (blood), tumor, and heart.
The drug is dosed into plasma, distributes to the tumor and heart, flows back,
and is eliminated from plasma.

    dCp/dt = -(k_pt + k_ph + k_el)*Cp + k_tp*Ct + k_hp*Ch
    dCt/dt =  k_pt*Cp - k_tp*Ct
    dCh/dt =  k_ph*Cp - k_hp*Ch

## The virtual trial
Each virtual patient gets slightly different rate constants (log-normal
variation). For 1,000 patients, the script checks whether the tumor reaches
an effective concentration while the heart stays below a toxicity limit.
The share of patients meeting both is the predicted probability of success.

## Run it
    pip install numpy scipy matplotlib
    python simulate.py

## Output
- Probability of success printed to the console
- concentration_curves.png: drug levels over 48 hours for a typical patient

## Notes
Parameters are illustrative, not fitted to a real drug.
