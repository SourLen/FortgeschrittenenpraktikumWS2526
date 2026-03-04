# plotting.py
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Callable, Sequence, Optional

@dataclass(frozen=True)
class PullPlotResult:
    chi2: float
    ndof: int
    chi2_red: float

def plot_data_fit_and_pulls(
    x: np.ndarray,
    y: np.ndarray,
    yerr: np.ndarray,
    model: Callable[[np.ndarray, Sequence[float]], np.ndarray],
    popt: Sequence[float],
    fit_mask: np.ndarray,
    title: str = "",
    xlabel: str = "channel",
    ylabel: str = "counts",
    context_pad: int | float = 0,
    vlines: Optional[list[tuple[float, str]]] = None,
    show: bool = True,
) -> PullPlotResult:
    """
    Top: data+fit (with optional context outside fit range)
    Bottom: pulls only for fit_mask
    """
    x = np.asarray(x)
    y = np.asarray(y)
    yerr = np.asarray(yerr)

    if vlines is None:
        vlines = []

    # context range
    if context_pad != 0:
        x_min = x[fit_mask].min() - context_pad
        x_max = x[fit_mask].max() + context_pad
        context_mask = (x >= x_min) & (x <= x_max)
    else:
        context_mask = fit_mask

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=False)

    # Top panel
    ax1.errorbar(x[context_mask], y[context_mask], yerr=yerr[context_mask], color="blue",fmt=".", label="Data")
    xfit = np.linspace(x[fit_mask].min(), x[fit_mask].max(), 1200)
    ax1.plot(xfit, model(xfit, popt), color="red", label="Fit")
    for xv, lab in vlines:
        ax1.axvline(xv, linestyle="--", label=lab, color="purple")

    ax1.set_title(title)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.legend()

    # Pulls (fit points only)
    resid = y[fit_mask] - model(x[fit_mask], popt)
    pulls = resid / yerr[fit_mask]
    chi2 = float(np.sum(pulls**2))
    ndof = int(np.sum(fit_mask) - len(popt))
    chi2_red = chi2 / ndof if ndof > 0 else float("nan")

    ax2.errorbar(x[fit_mask], pulls, yerr=np.ones_like(pulls), color="blue",fmt=".")
    ax2.axhline(0, color="red", linestyle="--")
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(r"$(y_i-f_i)/\sigma_i$")
    ax2.set_title(f"Pulls (χ²/ndf = {chi2_red:.2f})")

    plt.tight_layout()
    if show:
        plt.show()
    else:
        plt.close(fig)

    return PullPlotResult(chi2=chi2, ndof=ndof, chi2_red=chi2_red)