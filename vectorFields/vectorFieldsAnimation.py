#!/usr/bin/env python3
"""
Animated visualization of PAIRS of 2D vector fields (one pair per frame).

- Each frame shows exactly TWO vector fields (left/right subplots).
- You supply a provider function: vf_provider(frame_idx) -> (VF1, VF2, titles?)
- Works with SymPy expressions or strings in variables (x, y).
- Uses Matplotlib (quiver or streamplot).

Requires: numpy, sympy, matplotlib
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, Iterable, Optional, Tuple, Union, List
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from basicClasses.CommutatorSearchSymbolic import *

# ---------------------------
# Core data structures
# ---------------------------

SymLike = Union[str, sp.Expr]

class ConstantFunction:
    """Wrap a constant expression so it broadcasts over meshgrids cleanly."""
    def __init__(self, const_callable):
        self.const = const_callable
    def __call__(self, *args):
        if not args:
            raise AssertionError("ConstantFunction expects at least one array argument.")
        return self.const(*args) * np.ones_like(args[0], dtype=float)

@dataclass
class VectorField2D:
    """
    Vector field in the plane:
      F(x, y) = (P(x,y), Q(x,y))
    P, Q may be SymPy expressions or strings using variables x, y (and optional params).
    """
    P: SymLike
    Q: SymLike
    params: Optional[dict] = None  # e.g. {"a": 2.0}

    def _sympify(self) -> Tuple[sp.Expr, sp.Expr, sp.Symbol, sp.Symbol]:
        x, y = sp.symbols("x_0 x_1", real=True)
        local = {"x": x, "y": y}
        if self.params:
            for k in self.params:
                local[k] = sp.Symbol(k, real=True)
        P_expr = sp.sympify(self.P, locals=local)
        Q_expr = sp.sympify(self.Q, locals=local)
        return P_expr, Q_expr, x, y

    def _lambdas(self):
        """Return two NumPy-callables for evaluating on arrays."""
        P_expr, Q_expr, x, y = self._sympify()
        subs_params = {}
        if self.params:
            subs_params = {sp.Symbol(k): float(v) for k, v in self.params.items()}
            P_expr = P_expr.subs(subs_params)
            Q_expr = Q_expr.subs(subs_params)

        # If expression is constant w.r.t x,y, wrap it to broadcast correctly
        if P_expr.is_constant(x, y):
            P_fn = ConstantFunction(sp.lambdify((x, y), P_expr, modules="numpy"))
        else:
            P_fn = sp.lambdify((x, y), P_expr, modules="numpy")

        if Q_expr.is_constant(x, y):
            Q_fn = ConstantFunction(sp.lambdify((x, y), Q_expr, modules="numpy"))
        else:
            Q_fn = sp.lambdify((x, y), Q_expr, modules="numpy")
        return P_fn, Q_fn

    def evaluate_on_grid(self, X: np.ndarray, Y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        P_fn, Q_fn = self._lambdas()
        U = np.array(P_fn(X, Y), dtype=float)
        V = np.array(Q_fn(X, Y), dtype=float)
        U = np.nan_to_num(U, nan=0.0, posinf=0.0, neginf=0.0)
        V = np.nan_to_num(V, nan=0.0, posinf=0.0, neginf=0.0)
        return U, V


# ---------------------------
# Animator
# ---------------------------

class TwoFieldsAnimator:
    """
    Animate a sequence of frames; each frame shows TWO vector fields in 1×2 subplots.

    Parameters
    ----------
    vf_provider : callable
        vf_provider(frame_idx) -> (VF1, VF2, titles)
        - VF1, VF2 must be VectorField2D
        - titles is optional; if provided, should be a (title_left, title_right) tuple
    frames : int or iterable
        Number of frames or explicit sequence of frame indices.
    bounds : (xmin, xmax, ymin, ymax)
    density : int
        Number of points per axis (grid resolution).
    mode : {"stream", "quiver"}
        Plot style.
    normalize : bool
        If True (quiver mode), normalize vectors to unit length.
    interval : int
        Delay between frames (ms).
    save_path : Optional[str]
        If set (e.g., "out.mp4" or "out.gif"), animation is saved.
    writer : Optional[str]
        Matplotlib writer name (e.g., "ffmpeg", "imagemagick", "pillow").
    """

    def __init__(
        self,
        vf_provider: Callable[[int], Tuple[VectorField2D, VectorField2D, Optional[Tuple[str, str]]]],
        frames: Union[int, Iterable[int]],
        bounds: Tuple[float, float, float, float] = (-5, 5, -5, 5),
        density: int = 25,
        mode: str = "stream",
        normalize: bool = False,
        interval: int = 800,
        save_path: Optional[str] = None,
        writer: Optional[str] = None,
    ):
        self.vf_provider = vf_provider
        self.frames = range(frames) if isinstance(frames, int) else list(frames)
        self.bounds = bounds
        self.density = int(density)
        self.mode = mode
        self.normalize = normalize
        self.interval = int(interval)
        self.save_path = save_path
        self.writer = writer

        # Grid
        xmin, xmax, ymin, ymax = self.bounds
        x = np.linspace(xmin, xmax, self.density)
        y = np.linspace(ymin, ymax, self.density)
        self.X, self.Y = np.meshgrid(x, y)

        # Figure: 1 row, 2 columns
        self.fig, self.axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)
        for ax in self.axes:
            self._style_axis(ax)

        # First draw
        self._artists_cache = None

    @staticmethod
    def _style_axis(ax):
        ax.spines["bottom"].set_position("zero")
        ax.spines["left"].set_position("zero")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, lw=0.3, alpha=0.4)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x")
        ax.set_ylabel("y")

    def _draw_field(self, ax, U, V, title: Optional[str]):
        ax.clear()
        self._style_axis(ax)

        if self.mode == "stream":
            ax.streamplot(self.X[0], self.Y[:, 0], U, V, density=1.1, linewidth=0.8, arrowsize=1.0)
        elif self.mode == "quiver":
            if self.normalize:
                M = np.hypot(U, V)
                M[M == 0] = 1.0
                U = U / M
                V = V / M
            ax.quiver(self.X, self.Y, U, V, angles="xy", scale_units="xy", pivot="mid")
        else:
            raise ValueError("mode must be 'stream' or 'quiver'")

        if title:
            ax.set_title(title)

    def _frame(self, idx):
        vf1, vf2, titles = self.vf_provider(idx)
        U1, V1 = vf1.evaluate_on_grid(self.X, self.Y)
        U2, V2 = vf2.evaluate_on_grid(self.X, self.Y)

        t1, t2 = (titles if titles is not None else (f"Frame {idx} — Field 1", f"Frame {idx} — Field 2"))

        self._draw_field(self.axes[0], U1, V1, t1)
        self._draw_field(self.axes[1], U2, V2, t2)
        self.fig.suptitle(f"Vector fields — frame {idx}")
        print(f'frame {idx + 1} is completed')
        return self.axes

    def animate(self):
        anim = FuncAnimation(self.fig, self._frame, frames=self.frames, interval=self.interval, blit=False, repeat=True)
        if self.save_path:
            anim.save(self.save_path, writer=(self.writer or "pillow"))
            print(f"Animation saved to: {self.save_path}")
        plt.show()


# ---------------------------
# Example provider
# ---------------------------

def my_provider(frame_idx: int):
    """
    Returns two fields that vary with frame index.
    You can replace this with your own function that returns arbitrary pairs.
    """
    # Left: rotated linear field
    N = frame_idx

    k = np.random.randint(1,10)
    n = 0
    l = 0
    m = np.random.randint(1,10)

    alpha = 1
    beta = 1

    powers1 = [k, n]
    powers2 = [l, m]

    x, y = symbols("x"), symbols("y")
    variables = [x, y]

    monomial1 = alpha * x ** k * y ** n
    monomial2 = beta * x ** l * y ** m

    polynomial1 = Polynomial(poly_symbols=monomial1, variables=variables)
    polynomial2 = Polynomial(poly_symbols=monomial2, variables=variables)

    der = Derivation([polynomial1, polynomial2], variables)
    K = 2
    commutator = Commutator(der, K)
    res, isProportional = commutator.get_commutator()

    t = frame_idx
    VF1 = VectorField2D(P=monomial1, Q=monomial2)
    # Right: a nonlinear field with time-varying coefficient

    VF2 = VectorField2D(P=res.polynomials[0].polynomial_symbolic, Q=res.polynomials[1].polynomial_symbolic)
    titles = (f"F₁(x,y; t={t:.1f})", f"F₂(x,y; t={t:.1f})")

    print(f"Variables: {polynomial1.variables_polynom}")
    print("========Given derivation=======")
    for i in range(len(der.polynomials)):
        print(f'poly {i}: {der.polynomials[i].polynomial_symbolic}')

    print("==========Unknown derivation=======")

    print(f"proportional: {isProportional}")
    print(f"is Solution correct: {commutator.isSolution(derivation1=der, derivation2=res)}")
    print("COMMUTATOR")
    for i in range(len(res.polynomials)):
        print(f'poly {i}: {simplify(res.polynomials[i].polynomial_symbolic)}')
    print("=" * 100)
    return VF1, VF2, titles


# ---------------------------
# Run directly
# ---------------------------

if __name__ == "__main__":
    animator = TwoFieldsAnimator(
        vf_provider=my_provider,    # <— swap in your own function here
        frames=10,                       # number of frames
        bounds=(-5, 5, -5, 5),
        density=30,
        mode="stream",                   # "stream" or "quiver"
        normalize=False,                 # only used in "quiver" mode
        interval=1000,                    # ms
        save_path="animations/animation1.mp4",                  # e.g. "two_fields.gif" or "two_fields.mp4"
        writer="ffmpeg",                     # e.g. "ffmpeg", "imagemagick", "pillow"
    )
    animator.animate()
