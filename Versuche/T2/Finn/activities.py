import numpy as np
from uncertainties import ufloat
from uncertainties import umath

def calculate_activity(o_act, o_date, half_life):
    """
    Calculate the activity of a radioactive substance at a given date.

    Parameters:
    o_act (ufloat): The original activity of the substance with uncertainty, in kBq.
    o_date (str): The original date in the format 'YYYY-MM-DD'.
    half_life (ufloat): The half-life of the substance with uncertainty, in days.

    Returns:
    ufloat: The calculated activity at the current date with uncertainty.
    """
    from datetime import datetime

    # Convert string dates to datetime objects
    o_date = datetime.strptime(o_date, '%Y-%m-%d')
    current_date = datetime.now()

    # Calculate the time difference in days
    time_diff = (current_date - o_date).days 
    time_diff_unc = 12/np.sqrt(12)/24
    # Calculate the decay constant
    decay_constant = np.log(2) / half_life

    # Calculate the activity at the current date
    activity = o_act * umath.exp(-decay_constant * time_diff)

    return activity





names = [
    'Cs-137; 4/1B',
    'Cs-137; OI 797',
    'Co-60, LP 213',
    'Cs-137; MH 851',
    'Eu-152, MH 850',
    'Na-22, MH 852'
]
original_activities = [
    ufloat(44400, 0),
    ufloat(24800, 0),
    ufloat(37, 0),
    ufloat(37, 0),
    ufloat(37, 0),
    ufloat(37, 0)
]
original_dates = [
    '1964-04-28',
    '2008-10-01',
    '2003-04-15',
    '2004-03-10',
    '2004-03-10',
    '2004-03-10'
    ]
half_lives = [
    ufloat(11000, 90),
    ufloat(11000, 90),
    ufloat(1925.3, 0.4),
    ufloat(11000, 90),
    ufloat(4943, 5),
    ufloat(950.5, 0.4)
]
activities = []
for o_act, o_date, half_life in zip(original_activities, original_dates, half_lives):
    activity = calculate_activity(o_act, o_date, half_life)
    activities.append(activity)
print("Activity of the substances at the current date:")
for i, activity in enumerate(activities):
    print(f"{names[i]}: {activity}")