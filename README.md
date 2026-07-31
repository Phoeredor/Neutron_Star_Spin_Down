# Neutron-Star Spin-Down Evolution

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Wolfram Language](https://img.shields.io/badge/Wolfram%20Language-Scientific%20Computing-DD1100.svg)](https://www.wolfram.com/language/)
[![Physics](https://img.shields.io/badge/Physics-Neutron%20Stars-5C2D91.svg)](https://en.wikipedia.org/wiki/Neutron_star)

🌀 A compact theoretical-astrophysics project studying the rotational evolution of neutron stars under magnetic-dipole spin-down, including constant and exponentially decaying dipolar magnetic fields.

## Project report

The original project report is available here: **[Neutron-Star Spin-Down Report](report/neutron_star_spin_down_report.pdf)**.

## Physical model

The evolution is described by

$$
P\dot{P} = K B_{\mathrm d}^{2}(t),
$$

where $P$ is the rotation period, $\dot{P}$ is its time derivative, $B_{\mathrm d}(t)$ is the dipolar magnetic-field strength, and $K$ collects the model constants. The equation relates rotational slowing directly to the squared dipolar field.

The project considers three magnetic-field prescriptions:

- $B_{\mathrm d} = 10^{13}\,\mathrm{G}$;
- $B_{\mathrm d} = 5 \times 10^{14}\,\mathrm{G}$;
- $B_{\mathrm d}(t) = 5 \times 10^{14}\exp(-t/\tau)\,\mathrm{G}$, with $\tau = 10^4\,\mathrm{yr}$.

Each case is studied for initial periods of $1\,\mathrm{ms}$ and $3\,\mathrm{ms}$.

## Repository structure

The intended public repository layout is:

```text
.
├── README.md
├── LICENSE
├── .gitignore
├── report/
│   └── neutron_star_spin_down_report.pdf
├── notebooks/
│   └── neutron_star_spin_down.nb
├── src/
│   └── neutron_star_spin_down.wl
└── plots/
    └── generated figures
```

- `report/`: complete original project report in PDF format;
- `notebooks/`: curated Mathematica notebook;
- `src/`: cleaned and documented Wolfram Language implementation;
- `plots/`: final figures generated from the validated model.

## Project status

🔬 The complete PDF report is already included. The computational files and validated figures are being prepared for the public release, so reproducibility is not yet claimed.

## License

This repository is released under the [MIT License](LICENSE). © 2026 Federico Malara.
