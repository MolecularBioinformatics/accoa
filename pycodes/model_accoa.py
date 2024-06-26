# ode models for the acetylation of histones

def model_acetylcoa(t, y, pars):
    '''
    A simple model for the labelling of Acetyl-CoA and histone acetylation

    '''
    ac, l_ac = y

    d_ac_dt = pars['k0'].value - pars['k_de'].value * ac
    d_l_ac_dt = pars['k1'].value - pars['k_de'].value * l_ac
    return [d_ac_dt, d_l_ac_dt]


def model_acetylation_1site_uncorr(t, y, pars):
    '''
    A simple model for the labelling of Acetyl-CoA and histone acetylation

    '''
    ac, l_ac = 0.0, 1.0
    non_ac, na_ac12, na_ac13 = y

    k_d1 = pars['k_d1'].value
    k_a1 = pars['k_a1'].value

    d_nonac_dt = k_d1 * (na_ac12 + na_ac13) - (k_a1) * (non_ac * (ac + l_ac))

    d_na_ac12_dt = k_a1 * non_ac * ac - k_d1 * na_ac12

    d_na_ac13_dt = k_a1 * non_ac * l_ac - k_d1 * na_ac13

    return [d_nonac_dt, d_na_ac12_dt, d_na_ac13_dt]


def model_acetylation_1site(t, y, pars):
    '''
    A simple model for the labelling of Acetyl-CoA and histone acetylation

    '''
    ac, l_ac, non_ac, na_ac12, na_ac13 = y

    k_d1 = pars['k_d1'].value
    k_a1 = pars['k_a1'].value

    d_ac_dt = pars['k0'].value - pars['k_de'].value * ac
    d_l_ac_dt = pars['k1'].value - pars['k_de'].value * l_ac

    d_nonac_dt = k_d1 * (na_ac12 + na_ac13) - (k_a1) * (non_ac * (ac + l_ac))

    d_na_ac12_dt = k_a1 * non_ac * ac - k_d1 * na_ac12

    d_na_ac13_dt = k_a1 * non_ac * l_ac - k_d1 * na_ac13

    return [d_ac_dt, d_l_ac_dt, d_nonac_dt, d_na_ac12_dt, d_na_ac13_dt]


def model_acetylation_2sites_uncorr(t, y, pars, ac=0.0, l_ac=1.0):
    '''
    A simple model for the labelling of Acetyl-CoA and histone acetylation

    Parameters:
    - t: Time
    - y: List containing the concentrations of non-acetylated histone peptides and labeled peptides
    - pars: Dictionary containing the parameters of the model
    - ac: Concentration of Acetyl-CoA
    - l_ac: Concentration of histone peptide with lysine acetylation

    Returns:
    - List containing the rate of change of concentrations
    '''
    non_ac, na_ac12, ac12_na, na_ac13, ac13_na = y

    k_d1 = pars['k_d1'].value
    k_d2 = pars['k_d2'].value
    k_a1 = pars['k_a1'].value
    k_a2 = pars['k_a2'].value

    d_nonac_dt = k_d1 * (na_ac12 + na_ac13) + k_d2 * (ac12_na + + ac13_na) - \
        (k_a1 + k_a2) * non_ac * (ac + l_ac)

    d_na_ac12_dt = k_a1 * non_ac * ac - k_d1 * na_ac12

    d_ac12_na_dt = k_a2 * non_ac * ac - k_d2 * ac12_na

    d_na_ac13_dt = k_a1 * non_ac * l_ac - k_d1 * na_ac13

    d_ac13_na_dt = k_a2 * non_ac * l_ac - k_d2 * ac13_na

    return [d_nonac_dt, d_na_ac12_dt, d_ac12_na_dt, d_na_ac13_dt, d_ac13_na_dt]


def model_acetylation_2sites(t, y, pars):
    '''
    A simple model for the labelling of Acetyl-CoA and histone acetylation

    Parameters:
    - t: Time
    - y: List containing the concentrations of non-acetylated histone peptides and labeled peptides
    - pars: Dictionary containing the parameters of the model

    Returns:
    - List containing the rate of change of concentrations
    '''
    ac, l_ac, non_ac, na_ac12, ac12_na, na_ac13, ac13_na = y

    k_d1 = pars['k_d1'].value
    k_d2 = pars['k_d2'].value
    k_a1 = pars['k_a1'].value
    k_a2 = pars['k_a2'].value

    d_ac_dt = pars['k0'].value - pars['k_de'].value * ac
    d_l_ac_dt = pars['k1'].value - pars['k_de'].value * l_ac

    d_nonac_dt = k_d1 * (na_ac12 + na_ac13) + k_d2 * (ac12_na + + ac13_na) - \
        (k_a1 + k_a2) * non_ac * (ac + l_ac)

    d_na_ac12_dt = k_a1 * non_ac * ac - k_d1 * na_ac12

    d_ac12_na_dt = k_a2 * non_ac * ac - k_d2 * ac12_na

    d_na_ac13_dt = k_a1 * non_ac * l_ac - k_d1 * na_ac13

    d_ac13_na_dt = k_a2 * non_ac * l_ac - k_d2 * ac13_na

    return [d_ac_dt, d_l_ac_dt, d_nonac_dt, d_na_ac12_dt, d_ac12_na_dt, d_na_ac13_dt, d_ac13_na_dt]
