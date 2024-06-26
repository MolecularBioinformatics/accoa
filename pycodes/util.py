from scipy.integrate import solve_ivp
import numpy as np
import pandas as pd
from lmfit import minimize, report_fit


def residual(params, fun_ode_model, t, y0, data):
    """
    Calculate the residual between the model and the data

    Args:
        params (lmfit.parameter): model parameters to be optimized
        fun_ode_model (function): function that describes the ODE model
        t (list): time points
        y0 (list): initial conditions
        data (list): measured data

    Returns:
        array: residual between the model and the data
    """
    model = solve_ivp(fun_ode_model, t_span=[t[0], t[-1]],
                      y0=y0, args=(params,), t_eval=t)

    return (model.y.T - data).ravel()


def get_measured_data(d1, d2, carrier, cells):
    """
    Get measured data from the dataset

    Args:
        d1 (pd.DataFrame): dataset 1
        d2 (pd.DataFrame): dataset 2
        carrier (str): column name of the carrier
        cells (str): column name of the cells

    Returns:
        _type_: _description_
    """
    m_label = d1[(d1.carrier == carrier) & (d1.cells == cells) &
                 (d1.state == 'ac13_CoA')].groupby('time').mean(numeric_only=True).relative_label
    m_nolabel = d1[(d1.carrier == carrier) & (d1.cells == cells)
                   & (d1.state == 'ac12_CoA')].groupby('time').mean(numeric_only=True).relative_label
    measured_d1 = pd.concat((m_nolabel, m_label), axis=1)
    measured_d1.columns = ['ac12_CoA', 'ac13_CoA']
    measured_d1.index = m_label.index/60.
    measured_d2 = d2[(d2.carrier == carrier) & (d2.cells == cells) & (
        d2.rep == 1)].pivot(index='time', values='rel_area', columns='combination')
    measured = pd.concat((measured_d1, measured_d2), axis=1)
    return measured


def parameter_estimation(fun_ode_model, residual, params, d1, d2, carrier='DMSO', cells='TSCctrl',
                         method='least_squares', display_statistics=True, uncorr=False,
                         columns=['nolabel', 'label', 'non_ac', 'ac*_ac12', 'ac12_ac*',
                                  'ac*_ac13', 'ac13_ac*']):
    """
    Estimate the parameters of the ODE model

    Args:
        fun_ode_model (function): function that describes the ODE model
        residual (function): function that calculates the residual between the model and the data
        params (lmfit.parameter): model parameters to be optimized
        d1 (pd.DataFrame): dataset 1
        d2 (pd.DataFrame): dataset 2
        carrier (str, optional): _description_. Defaults to 'DMSO'.
        cells (str, optional): _description_. Defaults to 'TSCctrl'.
        method (str, optional): _description_. Defaults to 'least_squares'.
        display_statistics (bool, optional): _description_. Defaults to True.
        uncorr (bool, optional): _description_. Defaults to False.
        columns (list, optional): _description_. Defaults to ['nolabel', 'label', 'non_ac', 'ac*_ac12', 'ac12_ac*', 'ac*_ac13', 'ac13_ac*'].

    Returns:
        tuple: optimized parameters and statistics
    """

    # measured data
    measured = get_measured_data(d1=d1, d2=d2, carrier=carrier, cells=cells)
    if uncorr:
        t_measured = measured[columns[2:]].index
        x2_measured = measured[columns[2:]].values
        y0 = measured[columns[2:]].iloc[0]
        fit_colnames = columns[2:]
    else:
        t_measured = measured[columns].index
        x2_measured = measured[columns].values
        y0 = measured[columns].iloc[0]
        fit_colnames = columns

    # fit model
    result = minimize(residual, params, args=(fun_ode_model, t_measured, y0, x2_measured),
                      method=method)  # leastsq nelder

    # check results of the fit
    sol = solve_ivp(fun=fun_ode_model, t_span=[t_measured[0], t_measured[-1]], y0=y0,
                    t_eval=np.linspace(t_measured[0], t_measured[-1], 100),  args=(result.params,))
    df = pd.DataFrame(sol.y.T, columns=fit_colnames)
    df['time'] = sol.t
    df = df.set_index('time')

    if display_statistics:
        return report_fit(result)
    else:
        return df, result
