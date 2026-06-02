"""
navier_stokes_branching.py
--------------------------
Finite-difference Navier-Stokes-inspired complex branching renderer.

The original renderer sampled a closed-form toy residual directly, which
made the deeper recursion levels contract into a nearly uniform field.
This version builds a periodic 2D vorticity field on an auxiliary grid,
evolves it with a finite-difference advection-diffusion stepper, solves
streamfunction / pressure Poisson problems numerically, and then evaluates
the branching kernel by midpoint quadrature along short streamlines.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parent.parent
OUTPUTS_DIR = ROOT_DIR / "outputs"
FLOW_EXTENT = (-3.0, 3.0)
_FLOW_CACHE: dict[tuple[float, ...], dict[str, np.ndarray | float]] = {}


def branch_pow(z, a, beta, cut_angle=0.0):
    """Return (z - a)^beta with a rotatable branch cut."""
    w = z - a
    mag = np.abs(w) ** beta
    phase = np.angle(w)
    phase_wrapped = (phase - cut_angle + np.pi) % (2 * np.pi) - np.pi + cut_angle
    return mag * np.exp(1j * beta * phase_wrapped)


def finite_difference_gradient(field, dx, dy):
    """Second-order periodic central differences."""
    grad_x = (np.roll(field, -1, axis=1) - np.roll(field, 1, axis=1)) / (2.0 * dx)
    grad_y = (np.roll(field, -1, axis=0) - np.roll(field, 1, axis=0)) / (2.0 * dy)
    return grad_x, grad_y


def finite_difference_laplacian(field, dx, dy):
    """Periodic five-point Laplacian."""
    return (
        (np.roll(field, -1, axis=1) - 2.0 * field + np.roll(field, 1, axis=1)) / dx**2
        + (np.roll(field, -1, axis=0) - 2.0 * field + np.roll(field, 1, axis=0)) / dy**2
    )


def bilinear_sample_periodic(field, x_query, y_query, x_min, y_min, dx, dy):
    """Sample a periodic field with bilinear interpolation."""
    ny, nx = field.shape
    x_idx = np.mod((x_query - x_min) / dx, nx)
    y_idx = np.mod((y_query - y_min) / dy, ny)

    x0 = np.floor(x_idx).astype(int)
    y0 = np.floor(y_idx).astype(int)
    x1 = (x0 + 1) % nx
    y1 = (y0 + 1) % ny

    tx = x_idx - x0
    ty = y_idx - y0

    f00 = field[y0, x0]
    f10 = field[y0, x1]
    f01 = field[y1, x0]
    f11 = field[y1, x1]

    return (
        (1.0 - tx) * (1.0 - ty) * f00
        + tx * (1.0 - ty) * f10
        + (1.0 - tx) * ty * f01
        + tx * ty * f11
    )


def poisson_solve_periodic(rhs, dx, dy, iterations):
    """Solve Δu = rhs on a periodic grid with Jacobi relaxation."""
    rhs = rhs - np.mean(rhs)
    solution = np.zeros_like(rhs, dtype=float)
    dx2 = dx * dx
    dy2 = dy * dy
    denom = 2.0 * (dx2 + dy2)

    for _ in range(iterations):
        solution = (
            (np.roll(solution, -1, axis=1) + np.roll(solution, 1, axis=1)) * dy2
            + (np.roll(solution, -1, axis=0) + np.roll(solution, 1, axis=0)) * dx2
            - rhs * dx2 * dy2
        ) / denom
        solution -= np.mean(solution)

    return solution


def initial_vorticity_field(X, Y, frequency, time):
    """Seed a smooth, asymmetric vorticity field for the numerical solve."""
    radius = np.hypot(X, Y)
    theta = np.angle(X + 1j * Y)

    omega = (
        np.sin(frequency * X + time) * np.cos((frequency + 0.35) * Y - time)
        + 0.55 * np.cos(2.0 * frequency * X - 0.5 * time)
        * np.sin((frequency - 0.2) * Y + 0.75 * time)
        + 0.25 * np.sin(frequency * radius - 1.5 * theta + 0.4 * time)
    )

    width_a = 0.35 + 0.12 / (frequency + 0.5)
    width_b = 0.45 + 0.10 / (frequency + 0.5)
    blob_a = np.exp(-((X - 1.15) ** 2 + (Y + 0.75) ** 2) / width_a**2)
    blob_b = np.exp(-((X + 1.35) ** 2 + (Y - 0.95) ** 2) / width_b**2)
    omega += 0.65 * (blob_a - blob_b)

    omega -= np.mean(omega)
    omega /= np.std(omega) + 1e-8
    return omega


def forcing_field(X, Y, frequency, time):
    """Mild forcing keeps the auxiliary flow structured at deeper recursion."""
    return (
        0.35 * np.sin(1.6 * frequency * X - 0.3 * time) * np.sin((frequency + 0.4) * Y + 0.2 * time)
        + 0.15 * np.cos(frequency * (X + Y) + 0.5 * time)
    )


def finite_difference_flow(viscosity=0.08, frequency=1.5, time=0.0,
                           fd_grid=144, solver_steps=28, poisson_iters=80):
    """Build and cache a periodic finite-difference flow state."""
    cache_key = (
        round(viscosity, 8),
        round(frequency, 8),
        round(time, 8),
        int(fd_grid),
        int(solver_steps),
        int(poisson_iters),
    )
    if cache_key in _FLOW_CACHE:
        return _FLOW_CACHE[cache_key]

    x_min, x_max = FLOW_EXTENT
    y_min, y_max = FLOW_EXTENT
    x = np.linspace(x_min, x_max, fd_grid, endpoint=False)
    y = np.linspace(y_min, y_max, fd_grid, endpoint=False)
    dx = (x_max - x_min) / fd_grid
    dy = (y_max - y_min) / fd_grid
    X, Y = np.meshgrid(x, y)

    omega = initial_vorticity_field(X, Y, frequency, time)
    forcing = forcing_field(X, Y, frequency, time)

    advect_dt = 0.42 * min(dx, dy) / (1.0 + 0.65 * frequency)
    diffusion_scale = min(viscosity * advect_dt, 0.18 * min(dx, dy) ** 2)

    for _ in range(solver_steps):
        psi = poisson_solve_periodic(-omega, dx, dy, poisson_iters)
        dpsi_dx, dpsi_dy = finite_difference_gradient(psi, dx, dy)
        vel_x = dpsi_dy
        vel_y = -dpsi_dx

        departure_x = X - advect_dt * vel_x
        departure_y = Y - advect_dt * vel_y
        omega_adv = bilinear_sample_periodic(omega, departure_x, departure_y, x_min, y_min, dx, dy)
        omega = omega_adv + diffusion_scale * finite_difference_laplacian(omega_adv, dx, dy) + advect_dt * forcing
        omega -= np.mean(omega)
        omega /= np.std(omega) + 1e-8

    psi = poisson_solve_periodic(-omega, dx, dy, poisson_iters)
    dpsi_dx, dpsi_dy = finite_difference_gradient(psi, dx, dy)
    vel_x = dpsi_dy
    vel_y = -dpsi_dx

    dvelx_dx, dvelx_dy = finite_difference_gradient(vel_x, dx, dy)
    dvely_dx, dvely_dy = finite_difference_gradient(vel_y, dx, dy)
    pressure_rhs = -(dvelx_dx**2 + 2.0 * dvelx_dy * dvely_dx + dvely_dy**2)
    pressure = poisson_solve_periodic(pressure_rhs, dx, dy, poisson_iters)

    speed = np.hypot(vel_x, vel_y)
    omega_norm = omega / (np.std(omega) + 1e-8)
    pressure_norm = pressure / (np.std(pressure) + 1e-8)
    speed_norm = (speed - np.mean(speed)) / (np.std(speed) + 1e-8)

    scalar = 0.65 * np.tanh(0.95 * omega_norm) + 0.35j * np.tanh(0.85 * pressure_norm + 0.25 * speed_norm)
    state = {
        "x_min": x_min,
        "y_min": y_min,
        "dx": dx,
        "dy": dy,
        "vel_x": vel_x,
        "vel_y": vel_y,
        "scalar": scalar,
    }
    _FLOW_CACHE[cache_key] = state
    return state


def navier_stokes_scalar(z, viscosity=0.08, frequency=1.5, time=0.0,
                         fd_grid=144, solver_steps=28, poisson_iters=80):
    """Sample the finite-difference scalar field at complex coordinates."""
    flow_state = finite_difference_flow(
        viscosity=viscosity,
        frequency=frequency,
        time=time,
        fd_grid=fd_grid,
        solver_steps=solver_steps,
        poisson_iters=poisson_iters,
    )
    return bilinear_sample_periodic(
        flow_state["scalar"],
        np.real(z),
        np.imag(z),
        flow_state["x_min"],
        flow_state["y_min"],
        flow_state["dx"],
        flow_state["dy"],
    )


def phase_integral(z, viscosity=0.08, frequency=1.5, time=0.0,
                   fd_grid=144, solver_steps=28, poisson_iters=80,
                   quad_steps=7, path_scale=0.8):
    """Midpoint quadrature of the complex phase along short streamlines."""
    flow_state = finite_difference_flow(
        viscosity=viscosity,
        frequency=frequency,
        time=time,
        fd_grid=fd_grid,
        solver_steps=solver_steps,
        poisson_iters=poisson_iters,
    )
    x_min = flow_state["x_min"]
    y_min = flow_state["y_min"]
    dx = flow_state["dx"]
    dy = flow_state["dy"]
    vel_x = flow_state["vel_x"]
    vel_y = flow_state["vel_y"]
    scalar = flow_state["scalar"]

    x = np.real(z).copy()
    y = np.imag(z).copy()
    integral = np.zeros_like(z, dtype=complex)
    segment = 1.0 / quad_steps

    for _ in range(quad_steps):
        vel_x0 = bilinear_sample_periodic(vel_x, x, y, x_min, y_min, dx, dy)
        vel_y0 = bilinear_sample_periodic(vel_y, x, y, x_min, y_min, dx, dy)
        mid_x = x + 0.5 * segment * path_scale * vel_x0
        mid_y = y + 0.5 * segment * path_scale * vel_y0

        scalar_mid = bilinear_sample_periodic(scalar, mid_x, mid_y, x_min, y_min, dx, dy)
        integral += np.exp(1j * 2.0 * np.pi * scalar_mid)

        vel_x_mid = bilinear_sample_periodic(vel_x, mid_x, mid_y, x_min, y_min, dx, dy)
        vel_y_mid = bilinear_sample_periodic(vel_y, mid_x, mid_y, x_min, y_min, dx, dy)
        x += segment * path_scale * vel_x_mid
        y += segment * path_scale * vel_y_mid

    return segment * integral


def ns_branch_kernel(z, a=0.0, beta=0.5, cut_angle=0.0, viscosity=0.08,
                     frequency=1.5, time=0.0, fd_grid=144, solver_steps=28,
                     poisson_iters=80, quad_steps=7, path_scale=0.8):
    """One branching step driven by the Navier-Stokes-inspired scalar."""
    branched = branch_pow(z, a, beta, cut_angle)
    return phase_integral(
        branched,
        viscosity=viscosity,
        frequency=frequency,
        time=time,
        fd_grid=fd_grid,
        solver_steps=solver_steps,
        poisson_iters=poisson_iters,
        quad_steps=quad_steps,
        path_scale=path_scale,
    )


def ns_branch_recursion(z, depth, **kwargs):
    """Apply the kernel repeatedly to build a recursive complex field."""
    for _ in range(depth):
        z = ns_branch_kernel(z, **kwargs)
    return z


def domain_color(z, title="", ax=None):
    """Render a complex array with perceptually-uniform domain coloring.

    Phase → twilight_shifted (cyclic, muted, not rainbow).
    Magnitude → brightness with subtle isochromatic contour rings.
    """
    angle = np.angle(z)
    log_mag = np.log10(np.abs(z) + 1e-300)

    phase_norm = np.nan_to_num((angle + np.pi) / (2 * np.pi), nan=0.0)
    rgb = plt.get_cmap("twilight_shifted")(phase_norm)[..., :3]

    clipped = np.clip(log_mag, -2, 2)
    brightness = np.nan_to_num(0.30 + 0.70 * (clipped + 2) / 4, nan=0.30)
    frac = log_mag - np.floor(log_mag)
    ring = 0.80 + 0.20 * np.cos(2 * np.pi * frac)
    modulation = np.clip(ring * brightness, 0, 1)[..., np.newaxis]
    rgb = np.clip(rgb * modulation, 0, 1)

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(rgb, origin="lower", extent=[-3, 3, -3, 3], interpolation="bilinear")
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("Re")
    ax.set_ylabel("Im")
    ax.set_aspect("equal")
    return ax


def make_grid(resolution):
    x = np.linspace(-3, 3, resolution)
    y = np.linspace(-3, 3, resolution)
    X, Y = np.meshgrid(x, y)
    return X + 1j * Y


def render_sweep_panel(grid, depths, kernel_kwargs, title, output_path):
    fig, axes = plt.subplots(2, 3, figsize=(15, 10), facecolor="white")
    fig.suptitle(title, fontsize=12, color="#222222")
    for ax, depth in zip(axes.flat, depths):
        print(f"  depth={depth} ...", flush=True)
        field = ns_branch_recursion(grid.copy(), depth, **kernel_kwargs)
        domain_color(field, title=f"depth = {depth}", ax=ax)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(output_path, dpi=130, facecolor="white")
    print(f"Saved {output_path}")


def render_parameter_panel(grid, values, value_name, fixed_kwargs, output_path):
    fig, axes = plt.subplots(2, 3, figsize=(15, 10), facecolor="white")
    fig.suptitle(f"{value_name} sweep", fontsize=12, color="#222222")
    for ax, value in zip(axes.flat, values):
        print(f"  {value_name}={value:.3f} ...", flush=True)
        field = ns_branch_recursion(grid.copy(), 4, **{**fixed_kwargs, value_name: value})
        domain_color(field, title=f"{value_name} = {value:.3f}", ax=ax)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(output_path, dpi=130, facecolor="white")
    print(f"Saved {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Render Navier-Stokes-inspired branching sweeps.")
    parser.add_argument("--resolution", type=int, default=220, help="Grid resolution for sweeps")
    parser.add_argument("--hero-resolution", type=int, default=380, help="Grid resolution for the hero shot")
    parser.add_argument("--depth", type=int, default=5, help="Depth for the hero shot")
    parser.add_argument("--beta", type=float, default=0.5, help="Branch exponent")
    parser.add_argument("--cut-angle", type=float, default=0.0, help="Branch-cut rotation in radians")
    parser.add_argument("--viscosity", type=float, default=0.08, help="Viscosity for the toy flow field")
    parser.add_argument("--frequency", type=float, default=1.5, help="Spatial frequency for the toy flow field")
    parser.add_argument("--time", type=float, default=0.0, help="Phase shift for the toy flow field")
    parser.add_argument("--fd-grid", type=int, default=144, help="Auxiliary grid size for the finite-difference flow solve")
    parser.add_argument("--solver-steps", type=int, default=28, help="Advection-diffusion steps for the auxiliary flow solve")
    parser.add_argument("--poisson-iters", type=int, default=80, help="Jacobi iterations for periodic Poisson solves")
    parser.add_argument("--quad-steps", type=int, default=7, help="Midpoint quadrature segments for the branching kernel")
    parser.add_argument("--path-scale", type=float, default=0.8, help="Streamline length used inside the numerical phase integral")
    parser.add_argument("--skip-show", action="store_true", help="Do not open a plot window")
    args = parser.parse_args()

    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)

    sweep_grid = make_grid(args.resolution)
    hero_grid = make_grid(args.hero_resolution)

    base_kwargs = dict(
        a=0.0,
        beta=args.beta,
        cut_angle=args.cut_angle,
        viscosity=args.viscosity,
        frequency=args.frequency,
        time=args.time,
        fd_grid=args.fd_grid,
        solver_steps=args.solver_steps,
        poisson_iters=args.poisson_iters,
        quad_steps=args.quad_steps,
        path_scale=args.path_scale,
    )

    depths = [1, 2, 3, 4, 5, 6]
    render_sweep_panel(
        sweep_grid,
        depths,
        base_kwargs,
        title=f"Navier-Stokes branching depth evolution (finite-difference, beta={args.beta}, nu={args.viscosity})",
        output_path=OUTPUTS_DIR / "navier_stokes_depth_evolution.png",
    )

    viscosities = [0.01, 0.03, 0.08, 0.12, 0.20, 0.35]
    render_parameter_panel(
        sweep_grid,
        viscosities,
        "viscosity",
        dict(
            a=0.0,
            beta=args.beta,
            cut_angle=args.cut_angle,
            frequency=args.frequency,
            time=args.time,
            fd_grid=args.fd_grid,
            solver_steps=args.solver_steps,
            poisson_iters=args.poisson_iters,
            quad_steps=args.quad_steps,
            path_scale=args.path_scale,
        ),
        OUTPUTS_DIR / "navier_stokes_viscosity_sweep.png",
    )

    frequencies = [0.75, 1.0, 1.5, 2.0, 2.5, 3.0]
    render_parameter_panel(
        sweep_grid,
        frequencies,
        "frequency",
        dict(
            a=0.0,
            beta=args.beta,
            cut_angle=args.cut_angle,
            viscosity=args.viscosity,
            time=args.time,
            fd_grid=args.fd_grid,
            solver_steps=args.solver_steps,
            poisson_iters=args.poisson_iters,
            quad_steps=args.quad_steps,
            path_scale=args.path_scale,
        ),
        OUTPUTS_DIR / "navier_stokes_frequency_sweep.png",
    )

    print(f"  hero shot depth={args.depth} ...", flush=True)
    hero_field = ns_branch_recursion(hero_grid.copy(), args.depth, **base_kwargs)
    fig, ax = plt.subplots(figsize=(9, 8))
    domain_color(hero_field, title=f"Navier-Stokes branch depth={args.depth}, beta={args.beta}", ax=ax)
    plt.tight_layout()
    hero_path = OUTPUTS_DIR / "navier_stokes_hero.png"
    fig.savefig(hero_path, dpi=150)
    print(f"Saved {hero_path}")

    print("\nAll plots saved.")
    if not args.skip_show:
        plt.show()


if __name__ == "__main__":
    main()
