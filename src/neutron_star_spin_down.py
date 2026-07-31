#!/usr/bin/env python3
"""Generate validated neutron-star spin-down figures."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Final

import matplotlib
import numpy as np
import numpy.typing as npt

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


# Physical constants and model parameters in CGS units.
SECONDS_PER_YEAR: Final[float] = 365.25 * 24.0 * 3600.0
SPIN_DOWN_CONSTANT: Final[float] = 2.44e-40  # s G^-2
WEAK_MAGNETIC_FIELD: Final[float] = 1.0e13  # G
STRONG_MAGNETIC_FIELD: Final[float] = 5.0e14  # G
DECAY_TIMESCALE_YEARS: Final[float] = 1.0e4  # yr
MAXIMUM_TIME_YEARS: Final[float] = 1.0e6  # yr
INITIAL_PERIODS: Final[tuple[float, float]] = (1.0e-3, 3.0e-3)  # s

STELLAR_RADIUS: Final[float] = 1.0e6  # cm
MOMENT_OF_INERTIA: Final[float] = 1.0e45  # g cm^2
SOLAR_MASS: Final[float] = 2.0e33  # g
STELLAR_MASS: Final[float] = 1.25 * SOLAR_MASS  # g
GRAVITATIONAL_CONSTANT: Final[float] = 6.67259e-8  # cm^3 g^-1 s^-2
PDOT_CUTOFF: Final[float] = 1.0e-14  # s/s

ArrayLike = npt.ArrayLike
FloatArray = npt.NDArray[np.float64]


@dataclass(frozen=True)
class EvolutionTrack:
    """Computed quantities for one magnetic-field prescription."""

    period: FloatArray
    period_derivative: FloatArray
    relative_variation: FloatArray
    normalized_derivative: FloatArray
    asymptotic_period: float | None = None


def _time_array(time_seconds: ArrayLike) -> FloatArray:
    """Return a finite, non-negative time array in seconds."""

    values = np.asarray(time_seconds, dtype=np.float64)
    if np.any(~np.isfinite(values)):
        raise ValueError("Time values must all be finite.")
    if np.any(values < 0.0):
        raise ValueError("Time values must be non-negative.")
    return values


def _positive_scalar(value: float, name: str) -> float:
    """Validate and return a finite positive scalar."""

    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive.")
    return float(value)


def constant_field_period(
    time_seconds: ArrayLike,
    initial_period: float,
    magnetic_field: float,
    spin_down_constant: float = SPIN_DOWN_CONSTANT,
) -> FloatArray:
    """Return the period for a constant dipolar magnetic field."""

    time = _time_array(time_seconds)
    period_zero = _positive_scalar(initial_period, "Initial period")
    field = _positive_scalar(magnetic_field, "Magnetic field")
    coefficient = _positive_scalar(spin_down_constant, "Spin-down constant")
    return np.sqrt(period_zero**2 + 2.0 * coefficient * field**2 * time)


def constant_field_period_derivative(
    time_seconds: ArrayLike,
    initial_period: float,
    magnetic_field: float,
    spin_down_constant: float = SPIN_DOWN_CONSTANT,
) -> FloatArray:
    """Return Pdot for a constant dipolar magnetic field."""

    period = constant_field_period(
        time_seconds, initial_period, magnetic_field, spin_down_constant
    )
    field = _positive_scalar(magnetic_field, "Magnetic field")
    coefficient = _positive_scalar(spin_down_constant, "Spin-down constant")
    return coefficient * field**2 / period


def decaying_field_period(
    time_seconds: ArrayLike,
    initial_period: float,
    initial_magnetic_field: float,
    decay_timescale_seconds: float,
    spin_down_constant: float = SPIN_DOWN_CONSTANT,
) -> FloatArray:
    """Return the period for B(t) = B0 exp(-t/tau)."""

    time = _time_array(time_seconds)
    period_zero = _positive_scalar(initial_period, "Initial period")
    field_zero = _positive_scalar(initial_magnetic_field, "Initial magnetic field")
    timescale = _positive_scalar(decay_timescale_seconds, "Decay timescale")
    coefficient = _positive_scalar(spin_down_constant, "Spin-down constant")
    integrated_torque = coefficient * field_zero**2 * timescale
    return np.sqrt(
        period_zero**2
        + integrated_torque * (1.0 - np.exp(-2.0 * time / timescale))
    )


def decaying_field_period_derivative(
    time_seconds: ArrayLike,
    initial_period: float,
    initial_magnetic_field: float,
    decay_timescale_seconds: float,
    spin_down_constant: float = SPIN_DOWN_CONSTANT,
) -> FloatArray:
    """Return Pdot for B(t) = B0 exp(-t/tau)."""

    time = _time_array(time_seconds)
    period = decaying_field_period(
        time,
        initial_period,
        initial_magnetic_field,
        decay_timescale_seconds,
        spin_down_constant,
    )
    field_zero = _positive_scalar(initial_magnetic_field, "Initial magnetic field")
    timescale = _positive_scalar(decay_timescale_seconds, "Decay timescale")
    coefficient = _positive_scalar(spin_down_constant, "Spin-down constant")
    return coefficient * field_zero**2 * np.exp(-2.0 * time / timescale) / period


def angular_velocity(period_seconds: ArrayLike) -> FloatArray:
    """Return angular velocity 2 pi/P in radians per second."""

    period = np.asarray(period_seconds, dtype=np.float64)
    if np.any(~np.isfinite(period)) or np.any(period <= 0.0):
        raise ValueError("Periods must be finite and positive.")
    return 2.0 * np.pi / period


def relative_percentage_variation(
    period_seconds: ArrayLike, initial_period: float
) -> FloatArray:
    """Return 100 [P(t) - P0]/P0."""

    period = np.asarray(period_seconds, dtype=np.float64)
    period_zero = _positive_scalar(initial_period, "Initial period")
    return 100.0 * (period - period_zero) / period_zero


def normalized_period_derivative(period_derivative: ArrayLike) -> FloatArray:
    """Normalize Pdot by its initial value."""

    derivative = np.asarray(period_derivative, dtype=np.float64)
    if derivative.size == 0:
        raise ValueError("The period-derivative array must not be empty.")
    if np.any(~np.isfinite(derivative)) or np.any(derivative < 0.0):
        raise ValueError("Period derivatives must be finite and non-negative.")
    initial_derivative = _positive_scalar(
        float(derivative.flat[0]), "Initial period derivative"
    )
    return derivative / initial_derivative


def asymptotic_decaying_field_period(
    initial_period: float,
    initial_magnetic_field: float,
    decay_timescale_seconds: float,
    spin_down_constant: float = SPIN_DOWN_CONSTANT,
) -> float:
    """Return the t -> infinity period for an exponential field decay."""

    period_zero = _positive_scalar(initial_period, "Initial period")
    field_zero = _positive_scalar(initial_magnetic_field, "Initial magnetic field")
    timescale = _positive_scalar(decay_timescale_seconds, "Decay timescale")
    coefficient = _positive_scalar(spin_down_constant, "Spin-down constant")
    return float(
        np.sqrt(period_zero**2 + coefficient * field_zero**2 * timescale)
    )


def breakup_angular_velocity(
    gravitational_constant: float = GRAVITATIONAL_CONSTANT,
    stellar_mass: float = STELLAR_MASS,
    stellar_radius: float = STELLAR_RADIUS,
) -> float:
    """Estimate the centrifugal breakup angular velocity."""

    gravity = _positive_scalar(gravitational_constant, "Gravitational constant")
    mass = _positive_scalar(stellar_mass, "Stellar mass")
    radius = _positive_scalar(stellar_radius, "Stellar radius")
    return float(np.sqrt(gravity * mass / radius**3))


def breakup_period(
    gravitational_constant: float = GRAVITATIONAL_CONSTANT,
    stellar_mass: float = STELLAR_MASS,
    stellar_radius: float = STELLAR_RADIUS,
) -> float:
    """Estimate the centrifugal breakup period in seconds."""

    return float(
        2.0
        * np.pi
        / breakup_angular_velocity(gravitational_constant, stellar_mass, stellar_radius)
    )


def _time_grid() -> tuple[FloatArray, FloatArray]:
    """Return zero plus a dense logarithmic grid in years and seconds."""

    positive_years = np.geomspace(1.0e-6, MAXIMUM_TIME_YEARS, 3000)
    time_years = np.concatenate((np.array([0.0]), positive_years))
    return time_years, time_years * SECONDS_PER_YEAR


def _validate_track(
    label: str,
    time_seconds: FloatArray,
    initial_period: float,
    track: EvolutionTrack,
) -> None:
    """Raise a clear error if a computed evolution is inconsistent."""

    arrays = (
        track.period,
        track.period_derivative,
        track.relative_variation,
        track.normalized_derivative,
    )
    if any(np.any(~np.isfinite(values)) for values in arrays):
        raise ValueError(f"{label}: a generated array contains NaN or infinity.")
    if np.any(track.period <= 0.0):
        raise ValueError(f"{label}: every period must be positive.")
    if np.any(track.period_derivative < 0.0):
        raise ValueError(f"{label}: period derivatives must not be negative.")
    monotonic_tolerance = 64.0 * np.finfo(np.float64).eps * np.max(track.period)
    if np.any(np.diff(track.period) < -monotonic_tolerance):
        raise ValueError(f"{label}: the period decreases with time.")
    if not np.isclose(track.period[0], initial_period, rtol=0.0, atol=1.0e-15):
        raise ValueError(f"{label}: the initial period does not equal the requested P0.")
    if track.period.shape != time_seconds.shape:
        raise ValueError(f"{label}: the period and time grids have different shapes.")
    if track.asymptotic_period is not None:
        asymptote = track.asymptotic_period
        if np.any(track.period > asymptote * (1.0 + 1.0e-12)):
            raise ValueError(f"{label}: the decaying-field solution exceeds its asymptote.")
        if not np.isclose(track.period[-1], asymptote, rtol=1.0e-10, atol=1.0e-12):
            raise ValueError(
                f"{label}: the late-time solution disagrees with its analytical asymptote."
            )


def _compute_tracks(time_seconds: FloatArray, initial_period: float) -> dict[str, EvolutionTrack]:
    """Compute and validate all three magnetic-field prescriptions."""

    decay_timescale_seconds = DECAY_TIMESCALE_YEARS * SECONDS_PER_YEAR
    specifications = {
        "Constant weak field": ("constant", WEAK_MAGNETIC_FIELD),
        "Constant strong field": ("constant", STRONG_MAGNETIC_FIELD),
        "Exponentially decaying strong field": ("decaying", STRONG_MAGNETIC_FIELD),
    }
    tracks: dict[str, EvolutionTrack] = {}
    for label, (model, field) in specifications.items():
        if model == "constant":
            period = constant_field_period(time_seconds, initial_period, field)
            derivative = constant_field_period_derivative(
                time_seconds, initial_period, field
            )
            asymptote = None
        else:
            period = decaying_field_period(
                time_seconds, initial_period, field, decay_timescale_seconds
            )
            derivative = decaying_field_period_derivative(
                time_seconds, initial_period, field, decay_timescale_seconds
            )
            asymptote = asymptotic_decaying_field_period(
                initial_period, field, decay_timescale_seconds
            )
        track = EvolutionTrack(
            period=period,
            period_derivative=derivative,
            relative_variation=relative_percentage_variation(period, initial_period),
            normalized_derivative=normalized_period_derivative(derivative),
            asymptotic_period=asymptote,
        )
        _validate_track(label, time_seconds, initial_period, track)
        tracks[label] = track
    return tracks


FIELD_STYLES: Final[dict[str, dict[str, object]]] = {
    "Constant weak field": {"color": "#31688e", "linestyle": "-"},
    "Constant strong field": {"color": "#b63679", "linestyle": "--"},
    "Exponentially decaying strong field": {
        "color": "#35b779",
        "linestyle": "-.",
    },
}


def _configure_time_axis(axis: plt.Axes) -> None:
    """Use a valid symmetric-logarithmic time axis that includes zero."""

    axis.set_xscale("symlog", linthresh=1.0e-3, linscale=1.0)
    axis.set_xlim(0.0, MAXIMUM_TIME_YEARS)
    axis.set_xlabel("Time [yr]")
    axis.grid(True, which="both", alpha=0.25, linewidth=0.6)


def _save_figure(figure: plt.Figure, path: Path) -> None:
    """Save one publication-quality PNG and close its figure."""

    figure.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def _spin_evolution_figure(
    time_years: FloatArray,
    initial_period: float,
    tracks: dict[str, EvolutionTrack],
    output_path: Path,
) -> None:
    """Create the four-panel spin-evolution summary."""

    figure, axes = plt.subplots(2, 2, figsize=(11.5, 8.2), constrained_layout=True)
    quantities = (
        ("relative_variation", "Relative period variation [%]", False),
        ("period", "Period [s]", True),
        ("period_derivative", "Period derivative [s/s]", True),
        ("normalized_derivative", r"Normalized derivative $\dot{P}(t)/\dot{P}(0)$", True),
    )
    for axis, (attribute, y_label, logarithmic_y) in zip(axes.flat, quantities):
        for label, track in tracks.items():
            axis.plot(
                time_years,
                getattr(track, attribute),
                linewidth=2.0,
                label=label,
                **FIELD_STYLES[label],
            )
        _configure_time_axis(axis)
        axis.set_ylabel(y_label)
        if logarithmic_y:
            axis.set_yscale("log")
        axis.legend(frameon=False, fontsize=8)
    figure.suptitle(
        f"Neutron-star spin-down evolution ($P_0={initial_period * 1.0e3:.0f}$ ms)",
        fontsize=15,
    )
    _save_figure(figure, output_path)


def _initial_period_comparison_figure(
    time_years: FloatArray,
    all_tracks: dict[float, dict[str, EvolutionTrack]],
    output_path: Path,
) -> None:
    """Compare relative variation for two initial periods at fixed field."""

    figure, axis = plt.subplots(figsize=(8.2, 5.4), constrained_layout=True)
    colors = {1.0e-3: "#3b528b", 3.0e-3: "#e76f51"}
    for initial_period in INITIAL_PERIODS:
        track = all_tracks[initial_period]["Constant strong field"]
        axis.plot(
            time_years,
            track.relative_variation,
            color=colors[initial_period],
            linewidth=2.2,
            label=f"$P_0={initial_period * 1.0e3:.0f}$ ms",
        )
    _configure_time_axis(axis)
    axis.set_ylabel("Relative period variation [%]")
    axis.set_title(r"Initial-period comparison for constant $B=5\times10^{14}$ G")
    axis.legend(frameon=False)
    # A smaller P0 yields a larger normalized change, while the absolute
    # constant-field periods become less sensitive to P0 at late times.
    axis.text(
        0.02,
        0.97,
        "Smaller $P_0$ gives a larger relative variation;\n"
        "the absolute late-time periods progressively converge.",
        transform=axis.transAxes,
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.8, "edgecolor": "0.8"},
    )
    _save_figure(figure, output_path)


def _p_pdot_figure(
    all_tracks: dict[float, dict[str, EvolutionTrack]], output_path: Path
) -> None:
    """Create the physical P-Pdot diagram with an explicit display cutoff."""

    figure, axis = plt.subplots(figsize=(8.4, 6.0), constrained_layout=True)
    period_colors = {1.0e-3: "#2a6fbb", 3.0e-3: "#d95f02"}
    field_line_styles = {
        "Constant weak field": "-",
        "Constant strong field": "--",
        "Exponentially decaying strong field": "-.",
    }
    displayed_periods: list[FloatArray] = []
    for initial_period in INITIAL_PERIODS:
        for field_label, track in all_tracks[initial_period].items():
            visible = track.period_derivative >= PDOT_CUTOFF
            if not np.any(visible):
                raise ValueError(
                    f"No P-Pdot samples remain above the cutoff for {field_label}."
                )
            displayed_periods.append(track.period[visible])
            axis.plot(
                track.period[visible],
                track.period_derivative[visible],
                color=period_colors[initial_period],
                linestyle=field_line_styles[field_label],
                linewidth=2.0,
                label=(
                    f"$P_0={initial_period * 1.0e3:.0f}$ ms — "
                    f"{field_label.lower()}"
                ),
            )
    minimum_period = min(float(np.min(values)) for values in displayed_periods)
    maximum_period = max(float(np.max(values)) for values in displayed_periods)
    axis.hlines(
        PDOT_CUTOFF,
        minimum_period,
        maximum_period,
        colors="0.35",
        linestyles=":",
        linewidth=1.4,
        label=r"Adopted spin-down cutoff: $\dot{P}=10^{-14}$ s/s",
    )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("Period $P$ [s]")
    axis.set_ylabel(r"Period derivative $\dot{P}$ [s/s]")
    axis.set_title(r"Neutron-star evolution in the $P$–$\dot{P}$ plane")
    axis.grid(True, which="both", alpha=0.25, linewidth=0.6)
    axis.legend(frameon=False, fontsize=8, ncol=2)
    _save_figure(figure, output_path)


def _print_summary(
    all_tracks: dict[float, dict[str, EvolutionTrack]], output_paths: list[Path]
) -> None:
    """Print a concise numerical and output summary."""

    print("Neutron-star spin-down summary")
    print(f"  seconds per year: {SECONDS_PER_YEAR:.0f} s")
    print(f"  K: {SPIN_DOWN_CONSTANT:.3e} s G^-2")
    print(f"  weak magnetic field: {WEAK_MAGNETIC_FIELD:.3e} G")
    print(f"  strong magnetic field: {STRONG_MAGNETIC_FIELD:.3e} G")
    print(f"  decay timescale: {DECAY_TIMESCALE_YEARS:.3e} yr")
    print(f"  maximum evolution time: {MAXIMUM_TIME_YEARS:.3e} yr")
    print(f"  stellar radius: {STELLAR_RADIUS:.3e} cm")
    print(f"  moment of inertia: {MOMENT_OF_INERTIA:.3e} g cm^2")
    print(f"  stellar mass: {STELLAR_MASS:.3e} g")
    print(f"  P-Pdot display cutoff: {PDOT_CUTOFF:.3e} s/s")
    for initial_period in INITIAL_PERIODS:
        print(f"  initial period: {initial_period:.6f} s ({initial_period * 1.0e3:.0f} ms)")
        for label, track in all_tracks[initial_period].items():
            print(f"    {label}, P(1 Myr): {track.period[-1]:.9e} s")
        asymptote = all_tracks[initial_period][
            "Exponentially decaying strong field"
        ].asymptotic_period
        if asymptote is None:
            raise RuntimeError("The decaying-field asymptote was not computed.")
        print(f"    decaying-field asymptotic period: {asymptote:.9e} s")
    angular_breakup = breakup_angular_velocity()
    period_breakup = breakup_period()
    print(f"  breakup angular velocity: {angular_breakup:.9e} rad/s")
    print(f"  breakup period: {period_breakup:.9e} s ({period_breakup * 1.0e3:.6f} ms)")
    print("Generated figures:")
    for path in output_paths:
        print(f"  {path}")


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Generate validated neutron-star spin-down figures."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory, relative to the repository root by default.",
    )
    return parser.parse_args()


def main() -> None:
    """Validate the model, generate all figures, and print a summary."""

    arguments = parse_arguments()
    repository_root = Path(__file__).resolve().parents[1]
    if arguments.output_dir is None:
        output_directory = repository_root / "plots"
    elif arguments.output_dir.is_absolute():
        output_directory = arguments.output_dir
    else:
        output_directory = repository_root / arguments.output_dir
    output_directory.mkdir(parents=True, exist_ok=True)

    time_years, time_seconds = _time_grid()
    if np.any(~np.isfinite(time_years)) or np.any(np.diff(time_years) <= 0.0):
        raise ValueError("The time grid must be finite and strictly increasing.")
    all_tracks = {
        initial_period: _compute_tracks(time_seconds, initial_period)
        for initial_period in INITIAL_PERIODS
    }

    expected_paths = [
        output_directory / "spin_evolution_1ms.png",
        output_directory / "spin_evolution_3ms.png",
        output_directory / "initial_period_comparison.png",
        output_directory / "p_pdot_diagram.png",
    ]
    _spin_evolution_figure(
        time_years, INITIAL_PERIODS[0], all_tracks[INITIAL_PERIODS[0]], expected_paths[0]
    )
    _spin_evolution_figure(
        time_years, INITIAL_PERIODS[1], all_tracks[INITIAL_PERIODS[1]], expected_paths[1]
    )
    _initial_period_comparison_figure(time_years, all_tracks, expected_paths[2])
    _p_pdot_figure(all_tracks, expected_paths[3])
    _print_summary(all_tracks, expected_paths)


if __name__ == "__main__":
    main()
