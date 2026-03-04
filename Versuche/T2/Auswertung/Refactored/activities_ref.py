from dataclasses import dataclass
from datetime import datetime
import numpy as np
from uncertainties import ufloat, UFloat
from uncertainties import umath


@dataclass
class Source:
    name: str
    a_ref: UFloat          # reference activity (kBq)
    ref_date: str          # 'YYYY-MM-DD'
    half_life_days: UFloat

    def activity_on(self, date=None):
        """Return activity at a given date (default: today)."""
        if date is None:
            date = datetime.now()

        ref_date = datetime.strptime(self.ref_date, "%Y-%m-%d")

        dt_days = (date - ref_date).total_seconds() / (24 * 3600)
        decay_const = np.log(2) / self.half_life_days

        return self.a_ref * umath.exp(-decay_const * dt_days)


# --- Define all sources here ---

SOURCES = {
    "Cs-137_1": Source(
        name="Cs-137; 4/1B",
        a_ref=ufloat(44400, 0),
        ref_date="1964-04-28",
        half_life_days=ufloat(11000, 90),
    ),

    "Cs-137_2": Source(
        name="Cs-137; OI 797",
        a_ref=ufloat(24800, 0),
        ref_date="2008-10-01",
        half_life_days=ufloat(11000, 90),
    ),

    "Co-60": Source(
        name="Co-60, LP 213",
        a_ref=ufloat(37, 0),
        ref_date="2003-04-15",
        half_life_days=ufloat(1925.3, 0.4),
    ),

    "Cs-137_3": Source(
        name="Cs-137; MH 851",
        a_ref=ufloat(37, 0),
        ref_date="2004-03-10",
        half_life_days=ufloat(11000, 90),
    ),

    "Eu-152": Source(
        name="Eu-152, MH 850",
        a_ref=ufloat(37, 0),
        ref_date="2004-03-10",
        half_life_days=ufloat(4943, 5),
    ),

    "Na-22": Source(
        name="Na-22, MH 852",
        a_ref=ufloat(37, 0),
        ref_date="2004-03-10",
        half_life_days=ufloat(950.5, 0.4),
    ),
}