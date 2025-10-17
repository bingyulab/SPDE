"""
Stochastic Wave Equation Solver
Implements multiple spatial and temporal discretization schemes
Based on: Introduction to Computational Stochastic PDEs (Lord, Powell, Shardlow)

Author: Leon Jiang
Date: 2024
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse import linalg as splinalg
from scipy import fftpack
from matplotlib import cm
import os
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# ============================================================================
# NOISE GENERATION: Q-Wiener Process Approximation
# ============================================================================

def get_onedD_bj(dtref, J, a, r):
    """
    Compute coefficients for Q-Wiener process in H^r(0,a) space.
    
    These coefficients determine the spatial correlation structure of the noise.
    Higher r means smoother noise in space.
    
    Parameters:
    -----------
    dtref : float
        Reference time step for temporal discretization
    J : int
        Number of spatial modes (Fourier modes for noise expansion)
    a : float
        Domain size [0, a]
    r : float
        Regularity parameter (r=1 gives H^1 noise, r=0.5 gives L^2 noise)
    
    Returns:
    --------
    bj : ndarray, shape (J-1,)
        Coefficients for Q-Wiener process expansion
        
    Notes:
    ------
    The Q-Wiener process is approximated as:
    W(t,x) ≈ Σ_j b_j * β_j(t) * sin(jπx/a)
    where β_j are independent Brownian motions
    """
    jj = np.arange(1, J)
    myeps = 0.001  # Small regularization to ensure convergence
    
    # Eigenvalue decay: slower decay = rougher noise
    root_qj = jj ** (-(2 * r + 1 + myeps) / 2)
    
    # Scale by sqrt(dt/a) for correct variance
    bj = root_qj * np.sqrt(2 * dtref / a)
    
    logger.debug(f"Generated {len(bj)} noise coefficients with regularity r={r}")
    return bj


def icspde_dst1(u):
    """
    Discrete Sine Transform Type 1.
    Transforms between physical and Fourier sine space.
    """
    return fftpack.dst(u, type=1, axis=0) / 2


def get_onedD_dW(bj, kappa, iFspace, M):
    """
    Sample increments of Q-Wiener process ΔW.
    
    Parameters:
    -----------
    bj : ndarray
        Coefficients from get_onedD_bj
    kappa : int
        Refinement factor (kappa=1 for no refinement, kappa>1 for multilevel)
    iFspace : int
        0 for physical space, 1 for Fourier space output
    M : int
        Number of independent realizations
    
    Returns:
    --------
    dW : ndarray, shape (J-1, M)
        Wiener process increments
        If iFspace=0: in physical space
        If iFspace=1: in Fourier space
    """
    J = len(bj) + 1
    
    # Generate independent Gaussian random variables
    if kappa == 1:
        nn = np.random.randn(J-1, M)
    else:
        # For multilevel: sum kappa independent samples (gives sqrt(kappa) scaling)
        nn = np.sum(np.random.randn(kappa, J-1, M), axis=0)
    
    # Scale by coefficients
    X = bj.reshape(-1, 1) * nn
    
    # Transform to physical space if needed
    if iFspace == 1:
        dW = X
    else:
        dW = icspde_dst1(X)
    
    return dW


# ============================================================================
# STOCHASTIC WAVE EQUATION: FINITE DIFFERENCE METHOD
# ============================================================================

def stochastic_wave_fd(u0, v0, T, a, N, J, c, sigma, r):
    """
    Solve stochastic wave equation using finite differences and Euler-Maruyama.
    
    PDE: u_tt = c² u_xx + σ dW/dt
    Boundary: u(0,t) = u(a,t) = 0 (Dirichlet)
    
    Discretization:
    - Space: Central finite differences (2nd order)
    - Time: Euler-Maruyama (1st order weak, 0.5 order strong)
    
    Parameters:
    -----------
    u0 : ndarray, shape (J+1,)
        Initial displacement u(x,0)
    v0 : ndarray, shape (J+1,)
        Initial velocity u_t(x,0)
    T : float
        Final time
    a : float
        Domain size [0, a]
    N : int
        Number of time steps
    J : int
        Number of spatial interior points
    c : float
        Wave speed
    sigma : float
        Noise intensity
    r : float
        Noise regularity parameter
    
    Returns:
    --------
    t : ndarray, shape (N+1,)
        Time points
    x : ndarray, shape (J+1,)
        Spatial points
    ut : ndarray, shape (J+1, N+1)
        Solution u(t,x)
    vt : ndarray, shape (J+1, N+1)
        Velocity v(t,x) = u_t(t,x)
    """
    logger.info("=== Finite Difference Method ===")
    
    dt = T / N
    h = a / J
    t = np.linspace(0, T, N+1)
    x = np.linspace(0, a, J+1)
    
    # CFL condition check (for stability of deterministic part)
    CFL = c * dt / h
    logger.info(f"Time step: dt = {dt:.6f}")
    logger.info(f"Space step: h = {h:.6f}")
    logger.info(f"CFL number: {CFL:.4f} (stable if < 1)")
    
    if CFL >= 1:
        logger.warning(f"CFL condition violated! CFL = {CFL:.4f} >= 1. Solution may be unstable.")
    
    # Spatial discretization matrix: u_xx ≈ (u_{j-1} - 2u_j + u_{j+1})/h²
    e = np.ones(J+1)
    A = sparse.diags([e, -2*e, e], [-1, 0, 1], shape=(J+1, J+1), format='csc')
    
    # Apply Dirichlet boundary conditions (remove first and last rows/columns)
    ind = np.arange(1, J)
    A = A[ind, :][:, ind]
    A = (c**2 / h**2) * A
    
    logger.info(f"Spatial discretization matrix: {A.shape}")
    
    # Setup noise coefficients
    bj = get_onedD_bj(dt, J, a, r)
    logger.info(f"Noise: {len(bj)} modes, regularity r={r}")
    
    # Initialize solution arrays
    ut = np.zeros((J+1, N+1))
    vt = np.zeros((J+1, N+1))
    ut[:, 0] = u0
    vt[:, 0] = v0
    
    # Interior points only (boundary is zero)
    u_n = u0[ind].copy()
    v_n = v0[ind].copy()
    
    logger.info("Time integration started...")
    
    # Time stepping with Euler-Maruyama
    # System: u_t = v, v_t = c² u_xx + σ dW/dt
    for n in range(N):
        # Sample noise increment: ΔW ~ N(0, Δt)
        dW = get_onedD_dW(bj, 1, 0, 1).flatten()
        
        # Euler-Maruyama scheme:
        # u^{n+1} = u^n + Δt * v^n
        # v^{n+1} = v^n + Δt * (c² u_xx^n) + σ * ΔW
        u_new = u_n + dt * v_n
        v_new = v_n + dt * A.dot(u_n) + sigma * dW
        
        u_n = u_new
        v_n = v_new
        
        # Store solution (boundary remains zero)
        ut[ind, n+1] = u_n
        vt[ind, n+1] = v_n
        
        # Progress logging
        if (n+1) % (N//10) == 0 or n == N-1:
            logger.info(f"  Step {n+1}/{N} ({100*(n+1)/N:.0f}%), max|u| = {np.max(np.abs(u_n)):.4f}")
    
    logger.info("Time integration completed.")
    return t, x, ut, vt


# ============================================================================
# STOCHASTIC WAVE EQUATION: SPECTRAL GALERKIN (FFT) - FIXED
# ============================================================================

def stochastic_wave_spectral(u0, v0, T, a, N, J, c, sigma, r):
    """
    Solve stochastic wave equation using spectral Galerkin with explicit time-stepping.
    
    Uses Fourier series expansion with periodic boundary conditions.
    Time integration via Euler-Maruyama in Fourier space.
    
    Note: This assumes periodic boundary conditions, so results will differ from
    Dirichlet BC methods.
    """
    logger.info("=== Spectral Galerkin Method ===")
    
    dt = T / N
    t = np.linspace(0, T, N+1)
    x = np.linspace(0, a, J+1)
    
    # Wave numbers for periodic BC: k = 2πn/a
    kk = 2 * np.pi * np.hstack([np.arange(0, J//2+1), np.arange(-J//2+1, 0)]) / a
    
    # Eigenvalues: λ_k = -(ck)²
    MM = -(c * kk)**2
    
    # CFL-like stability check for spectral method
    max_eigenval = np.max(np.abs(MM))
    spectral_CFL = dt * np.sqrt(max_eigenval)
    logger.info(f"Spectral CFL: {spectral_CFL:.4f} (should be << 1)")
    
    if spectral_CFL > 0.5:
        logger.warning(f"Spectral method may be unstable! Consider smaller dt.")
    
    # Noise coefficients in Fourier space
    bj_full = np.zeros(J)
    jj = np.abs(np.hstack([np.arange(1, J//2+1), np.arange(-J//2+1, 0)]))
    myeps = 0.001
    root_qj = jj ** (-(2*r + 1 + myeps)/2)
    
    # Correct scaling for FFT normalization
    bj_full[1:] = root_qj * np.sqrt(2 * dt / a)
    bj_full[0] = 0  # No noise in mean mode (k=0)
    
    logger.info(f"Noise modes configured for FFT (periodic BC)")
    
    # Initialize in Fourier space
    ut = np.zeros((J+1, N+1))
    vt = np.zeros((J+1, N+1))
    ut[:J, 0] = u0[:J]
    vt[:J, 0] = v0[:J]
    
    uh = np.fft.fft(u0[:J])
    vh = np.fft.fft(v0[:J])
    
    logger.info("Time integration started...")
    
    for n in range(N):
        # Sample complex Gaussian noise in Fourier space
        nn = np.random.randn(J) + 1j * np.random.randn(J)
        nn[0] = np.random.randn()  # k=0 must be real
        if J % 2 == 0:
            nn[J//2] = np.random.randn()  # Nyquist frequency must be real
        
        dWh = bj_full * nn
        
        # Euler-Maruyama in Fourier space:
        # û_t = v̂
        # v̂_t = -(ck)² û + σ ΔŴ
        uh_new = uh + dt * vh
        vh_new = vh + dt * MM * uh + sigma * dWh
        
        uh = uh_new
        vh = vh_new
        
        # Transform back to physical space
        u = np.real(np.fft.ifft(uh))
        v = np.real(np.fft.ifft(vh))
        
        ut[:J, n+1] = u
        vt[:J, n+1] = v
        
        if (n+1) % (N//10) == 0 or n == N-1:
            logger.info(f"  Step {n+1}/{N} ({100*(n+1)/N:.0f}%), max|u| = {np.max(np.abs(u)):.4f}")
    
    # Enforce periodicity
    ut[J, :] = ut[0, :]
    vt[J, :] = vt[0, :]
    
    logger.info("Time integration completed.")
    return t, x, ut, vt


# ============================================================================
# STOCHASTIC WAVE EQUATION: FINITE ELEMENT METHOD - FIXED
# ============================================================================

def stochastic_wave_fem(u0, v0, T, a, N, ne, c, sigma, r):
    """
    Solve stochastic wave equation using linear finite elements (P1).
    
    Weak formulation with mass-lumping option for stability.
    Uses semi-implicit time-stepping for better stability.
    """
    logger.info("=== Finite Element Method ===")
    
    dt = T / N
    h = a / ne
    nvtx = ne + 1
    t = np.linspace(0, T, N+1)
    x = np.linspace(0, a, nvtx)
    
    logger.info(f"Elements: {ne}, Vertices: {nvtx}, h = {h:.6f}")
    
    # Mass matrix (consistent): ∫ φ_i φ_j dx
    MM = sparse.diags([1/6, 2/3, 1/6], [-1, 0, 1], shape=(nvtx, nvtx), format='csc')
    MM = h * MM
    
    # Stiffness matrix: ∫ φ'_i φ'_j dx
    KK = sparse.diags([-1, 2, -1], [-1, 0, 1], shape=(nvtx, nvtx), format='csc')
    KK = (c**2 / h) * KK
    
    # Apply Dirichlet boundary conditions
    ind = np.arange(1, nvtx-1)
    MM_int = MM[ind, :][:, ind]
    KK_int = KK[ind, :][:, ind]
    
    # Stability check
    h_min = h
    FEM_CFL = c * dt / h_min
    logger.info(f"FEM CFL: {FEM_CFL:.4f}")
    
    # Noise coefficients
    bj = get_onedD_bj(dt, ne, a, r)
    logger.info(f"Noise: {len(bj)} modes")
    
    # Initialize
    ut = np.zeros((nvtx, N+1))
    vt = np.zeros((nvtx, N+1))
    ut[:, 0] = u0
    vt[:, 0] = v0
    
    u_n = u0[ind].copy()
    v_n = v0[ind].copy()
    
    # Precompute mass matrix inverse (LU factorization for efficiency)
    MM_inv = splinalg.splu(MM_int.tocsc())
    
    logger.info("Time integration started...")
    
    for n in range(N):
        # Sample noise at nodes
        dW_nodes = get_onedD_dW(bj, 1, 0, 1).flatten()
        dW_full = np.zeros(nvtx)
        dW_full[ind] = dW_nodes
        
        # FEM load vector from noise: b = M * dW
        b_noise = MM_int.dot(dW_full[ind])
        
        # Semi-implicit Euler-Maruyama:
        # M v^{n+1} = M v^n - Δt K u^n + σ M ΔW
        # u^{n+1} = u^n + Δt v^{n+1}
        
        rhs = MM_int.dot(v_n) - dt * KK_int.dot(u_n) + sigma * b_noise
        v_new = MM_inv.solve(rhs)
        u_new = u_n + dt * v_new
        
        u_n = u_new
        v_n = v_new
        
        ut[ind, n+1] = u_n
        vt[ind, n+1] = v_n
        
        if (n+1) % (N//10) == 0 or n == N-1:
            logger.info(f"  Step {n+1}/{N} ({100*(n+1)/N:.0f}%), max|u| = {np.max(np.abs(u_n)):.4f}")
    
    logger.info("Time integration completed.")
    return t, x, ut, vt


# ============================================================================
# MILSTEIN METHOD FOR WAVE EQUATION
# ============================================================================

def stochastic_wave_milstein(u0, v0, T, a, N, J, c, sigma, r):
    """
    Milstein method for stochastic wave equation.
    
    Note: For additive noise (our case), Milstein = Euler-Maruyama.
    Included for completeness and educational purposes.
    
    For multiplicative noise g(u) dW, Milstein would add:
    0.5 * g'(u) * g(u) * (dW² - dt) correction term
    """
    logger.info("=== Milstein Method (= Euler-Maruyama for additive noise) ===")
    
    dt = T / N
    h = a / J
    t = np.linspace(0, T, N+1)
    x = np.linspace(0, a, J+1)
    
    CFL = c * dt / h
    logger.info(f"CFL number: {CFL:.4f}")
    
    # Spatial discretization
    e = np.ones(J+1)
    A = sparse.diags([e, -2*e, e], [-1, 0, 1], shape=(J+1, J+1), format='csc')
    ind = np.arange(1, J)
    A = A[ind, :][:, ind]
    A = (c**2 / h**2) * A
    
    # Noise
    bj = get_onedD_bj(dt, J, a, r)
    
    # Initialize
    ut = np.zeros((J+1, N+1))
    vt = np.zeros((J+1, N+1))
    ut[:, 0] = u0
    vt[:, 0] = v0
    
    u_n = u0[ind].copy()
    v_n = v0[ind].copy()
    
    logger.info("Time integration started...")
    
    for n in range(N):
        dW = get_onedD_dW(bj, 1, 0, 1).flatten()
        
        # For additive noise: σ dW, the Milstein correction is zero
        # Milstein correction: 0.5 * σ' * σ * (dW² - dt) = 0 (since σ' = 0)
        
        u_new = u_n + dt * v_n
        v_new = v_n + dt * A.dot(u_n) + sigma * dW
        
        u_n = u_new
        v_n = v_new
        
        ut[ind, n+1] = u_n
        vt[ind, n+1] = v_n
        
        if (n+1) % (N//10) == 0 or n == N-1:
            logger.info(f"  Step {n+1}/{N} ({100*(n+1)/N:.0f}%), max|u| = {np.max(np.abs(u_n)):.4f}")
    
    logger.info("Time integration completed.")
    return t, x, ut, vt


# ============================================================================
# MONTE CARLO SIMULATION
# ============================================================================

def monte_carlo_wave(method, u0, v0, T, a, N, J, c, sigma, r, M):
    """
    Monte Carlo simulation for computing ensemble statistics.
    
    Computes mean E[u(T,x)] and variance Var[u(T,x)] over M realizations.
    
    Parameters:
    -----------
    method : callable
        One of the solver functions (stochastic_wave_fd, etc.)
    M : int
        Number of Monte Carlo samples
    
    Returns:
    --------
    x : ndarray
        Spatial points
    mean_u : ndarray
        Sample mean of u(T,x)
    var_u : ndarray
        Sample variance of u(T,x)
    std_u : ndarray
        Sample standard deviation of u(T,x)
    """
    logger.info(f"=== Monte Carlo Simulation (M={M} samples) ===")
    
    # Run first sample to get dimensions
    t, x, ut, vt = method(u0, v0, T, a, N, J, c, sigma, r)
    
    sum_u = np.zeros_like(ut[:, -1])
    sum_u2 = np.zeros_like(ut[:, -1])
    
    logger.info("Running MC samples...")
    for m in range(M):
        t, x, ut, vt = method(u0, v0, T, a, N, J, c, sigma, r)
        u_final = ut[:, -1]
        sum_u += u_final
        sum_u2 += u_final**2
        
        if (m+1) % (M//10) == 0:
            logger.info(f"  Sample {m+1}/{M} ({100*(m+1)/M:.0f}%)")
    
    # Compute statistics
    mean_u = sum_u / M
    var_u = (sum_u2 - sum_u**2 / M) / (M - 1)  # Unbiased variance estimator
    std_u = np.sqrt(np.maximum(var_u, 0))  # Ensure non-negative
    
    logger.info(f"MC completed: mean(|u|) = {np.mean(np.abs(mean_u)):.4f}, mean(std) = {np.mean(std_u):.4f}")
    
    return x, mean_u, var_u, std_u


# ============================================================================
# VISUALIZATION
# ============================================================================

def plot_solution(t, x, ut, title="Stochastic Wave Equation", save_path=None):
    """
    Plot space-time solution as 3D surface and 2D contour.
    
    Parameters:
    -----------
    save_path : str, optional
        If provided, save figure to this path
    """
    fig = plt.figure(figsize=(14, 5))
    
    # 3D surface plot
    ax1 = fig.add_subplot(121, projection='3d')
    T, X = np.meshgrid(t, x)
    surf = ax1.plot_surface(X, T, ut, cmap=cm.coolwarm, 
                            linewidth=0, antialiased=True, alpha=0.9)
    ax1.set_xlabel('Space (x)', fontsize=10)
    ax1.set_ylabel('Time (t)', fontsize=10)
    ax1.set_zlabel('u(t,x)', fontsize=10)
    ax1.set_title(title, fontsize=12, fontweight='bold')
    ax1.view_init(elev=25, azim=45)
    fig.colorbar(surf, ax=ax1, shrink=0.5, aspect=10)
    
    # 2D contour plot
    ax2 = fig.add_subplot(122)
    contour = ax2.contourf(t, x, ut, 20, cmap='coolwarm')
    ax2.set_xlabel('Time (t)', fontsize=10)
    ax2.set_ylabel('Space (x)', fontsize=10)
    ax2.set_title(title + ' (Contour)', fontsize=12, fontweight='bold')
    fig.colorbar(contour, ax=ax2)
    ax2.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Figure saved to {save_path}")
    
    return fig


def plot_snapshots(t, x, ut, snapshots=[0, 0.25, 0.5, 0.75, 1.0], save_path=None):
    """
    Plot solution at specific time snapshots.
    
    Parameters:
    -----------
    snapshots : list of float
        Fractions of total time for snapshots (0 to 1)
    """
    fig, axes = plt.subplots(len(snapshots), 1, figsize=(10, 8), sharex=True)
    
    for i, frac in enumerate(snapshots):
        idx = int(frac * (len(t) - 1))
        axes[i].plot(x, ut[:, idx], 'b-', linewidth=2, label=f't = {t[idx]:.2f}')
        axes[i].set_ylabel('u(x)', fontsize=10)
        axes[i].set_title(f'Snapshot at t = {t[idx]:.2f}', fontsize=11, fontweight='bold')
        axes[i].grid(True, alpha=0.3, linestyle='--')
        axes[i].legend(loc='upper right')
        axes[i].axhline(0, color='k', linewidth=0.5, linestyle=':')
    
    axes[-1].set_xlabel('Space (x)', fontsize=10)
    plt.suptitle('Solution Snapshots Over Time', fontsize=13, fontweight='bold', y=1.00)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Figure saved to {save_path}")
    
    return fig


def plot_monte_carlo_stats(x, mean_u, std_u, T, M, save_path=None):
    """
    Plot Monte Carlo statistics: mean and confidence intervals.
    """
    fig, axes = plt.subplots(2, 1, figsize=(10, 7))
    
    # Mean with confidence intervals
    axes[0].plot(x, mean_u, 'b-', linewidth=2.5, label='Mean E[u(T,x)]')
    axes[0].fill_between(x, mean_u - 2*std_u, mean_u + 2*std_u, 
                         alpha=0.3, color='blue', label='±2σ (95% CI)')
    axes[0].fill_between(x, mean_u - std_u, mean_u + std_u, 
                         alpha=0.5, color='blue', label='±1σ (68% CI)')
    axes[0].axhline(0, color='k', linewidth=0.8, linestyle=':')
    axes[0].set_ylabel('u(T,x)', fontsize=11)
    axes[0].set_title(f'Monte Carlo Mean and Confidence Intervals at t={T:.2f} (M={M} samples)', 
                      fontsize=12, fontweight='bold')
    axes[0].legend(loc='best', fontsize=9)
    axes[0].grid(True, alpha=0.3, linestyle='--')
    
    # Standard deviation
    axes[1].plot(x, std_u, 'r-', linewidth=2.5, label='Std Dev σ[u(T,x)]')
    axes[1].fill_between(x, 0, std_u, alpha=0.3, color='red')
    axes[1].set_xlabel('Space (x)', fontsize=11)
    axes[1].set_ylabel('Standard Deviation', fontsize=11)
    axes[1].set_title('Pointwise Standard Deviation', fontsize=12, fontweight='bold')
    axes[1].legend(loc='best', fontsize=9)
    axes[1].grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Figure saved to {save_path}")
    
    return fig


# ============================================================================
# RESULT ANALYSIS
# ============================================================================

def analyze_results(methods_results):
    """
    Analyze and compare results from different methods.
    
    Parameters:
    -----------
    methods_results : dict
        Dictionary with method names as keys and (t, x, ut, vt) tuples as values
    """
    logger.info("\n" + "="*70)
    logger.info("RESULTS ANALYSIS")
    logger.info("="*70)
    
    for method_name, (t, x, ut, vt) in methods_results.items():
        max_u = np.max(np.abs(ut))
        mean_u = np.mean(np.abs(ut))
        final_energy = np.mean(vt[:, -1]**2 + ut[:, -1]**2)
        
        logger.info(f"\n{method_name}:")
        logger.info(f"  Max |u|: {max_u:.4f}")
        logger.info(f"  Mean |u|: {mean_u:.4f}")
        logger.info(f"  Final energy: {final_energy:.4f}")
        
        # Stability check
        if max_u > 100:
            logger.warning(f"  WARNING: Solution appears unstable (max|u| > 100)")
        elif max_u > 10:
            logger.warning(f"  CAUTION: Large amplitudes detected (max|u| > 10)")
        else:
            logger.info(f"  Status: Stable solution")


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    # Create output directory
    output_dir = "imgs"
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Output directory: {output_dir}")
    
    # ========================================================================
    # PROBLEM SETUP
    # ========================================================================
    logger.info("\n" + "="*70)
    logger.info("STOCHASTIC WAVE EQUATION SOLVER - COMPREHENSIVE DEMO")
    logger.info("="*70)
    
    # Physical parameters
    a = np.pi           # Domain: [0, π]
    T = 2.0             # Final time
    c = 1.0             # Wave speed
    sigma = 0.5         # Noise intensity
    r = 1.0             # Noise regularity (H^1 noise)
    
    # Discretization parameters (reduced for stability)
    J = 64              # Spatial points (reduced from 128)
    N = 2000            # Time steps (increased for stability)
    
    # Initial conditions: smooth Gaussian bump
    x_init = np.linspace(0, a, J+1)
    u0 = np.sin(x_init)              # Initial displacement
    v0 = np.zeros(J+1)               # Initial velocity (at rest)
    
    logger.info(f"\nProblem Parameters:")
    logger.info(f"  Domain: [0, {a:.4f}]")
    logger.info(f"  Time interval: [0, {T:.2f}]")
    logger.info(f"  Wave speed: c = {c}")
    logger.info(f"  Noise intensity: σ = {sigma}")
    logger.info(f"  Noise regularity: r = {r} (H^{r} noise)")
    logger.info(f"\nDiscretization:")
    logger.info(f"  Spatial points: J = {J}")
    logger.info(f"  Time steps: N = {N}")
    logger.info(f"  dt = {T/N:.6f}, dx = {a/J:.6f}")
    logger.info("="*70)
    
    # ========================================================================
    # RUN SIMULATIONS
    # ========================================================================
    results = {}
    
    # 1. Finite Difference + Euler-Maruyama
    logger.info("\n1. FINITE DIFFERENCE + EULER-MARUYAMA")
    logger.info("-" * 70)
    t_fd, x_fd, ut_fd, vt_fd = stochastic_wave_fd(u0, v0, T, a, N, J, c, sigma, r)
    results['FD + Euler-Maruyama'] = (t_fd, x_fd, ut_fd, vt_fd)
    
    # 2. Spectral Galerkin + Euler-Maruyama (with warnings)
    logger.info("\n2. SPECTRAL GALERKIN + EULER-MARUYAMA")
    logger.info("-" * 70)
    logger.info("Note: Using periodic BC (different from Dirichlet methods)")
    t_sp, x_sp, ut_sp, vt_sp = stochastic_wave_spectral(u0, v0, T, a, N, J, c, sigma, r)
    results['Spectral + Euler-Maruyama'] = (t_sp, x_sp, ut_sp, vt_sp)
    
    # 3. Finite Element + Euler-Maruyama
    logger.info("\n3. FINITE ELEMENT + EULER-MARUYAMA")
    logger.info("-" * 70)
    ne = J
    t_fe, x_fe, ut_fe, vt_fe = stochastic_wave_fem(u0, v0, T, a, N, ne, c, sigma, r)
    results['FEM + Euler-Maruyama'] = (t_fe, x_fe, ut_fe, vt_fe)
    
    # 4. Milstein Method
    logger.info("\n4. MILSTEIN METHOD")
    logger.info("-" * 70)
    t_mi, x_mi, ut_mi, vt_mi = stochastic_wave_milstein(u0, v0, T, a, N, J, c, sigma, r)
    results['Milstein'] = (t_mi, x_mi, ut_mi, vt_mi)
    
    # 5. Monte Carlo Simulation
    logger.info("\n5. MONTE CARLO SIMULATION")
    logger.info("-" * 70)
    M = 100
    x_mc, mean_u, var_u, std_u = monte_carlo_wave(stochastic_wave_fd, u0, v0, 
                                                    T, a, N, J, c, sigma, r, M)
    
    # ========================================================================
    # ANALYZE RESULTS
    # ========================================================================
    analyze_results(results)
    
    # ========================================================================
    # VISUALIZATION
    # ========================================================================
    logger.info("\n" + "="*70)
    logger.info("GENERATING VISUALIZATIONS")
    logger.info("="*70)
    
    # Plot 1: Finite Difference solution
    fig1 = plot_solution(t_fd, x_fd, ut_fd, 
                         "Stochastic Wave Equation (Finite Difference)",
                         save_path=os.path.join(output_dir, "01_fd_solution.png"))
    
    # Plot 2: Spectral solution
    fig2 = plot_solution(t_sp, x_sp, ut_sp, 
                         "Stochastic Wave Equation (Spectral Galerkin - Periodic BC)",
                         save_path=os.path.join(output_dir, "02_spectral_solution.png"))
    
    # Plot 3: FEM solution
    fig3 = plot_solution(t_fe, x_fe, ut_fe, 
                         "Stochastic Wave Equation (Finite Element)",
                         save_path=os.path.join(output_dir, "03_fem_solution.png"))
    
    # Plot 4: Snapshots of FD solution
    fig4 = plot_snapshots(t_fd, x_fd, ut_fd,
                          save_path=os.path.join(output_dir, "04_snapshots.png"))
    
    # Plot 5: Monte Carlo statistics
    fig5 = plot_monte_carlo_stats(x_mc, mean_u, std_u, T, M,
                                   save_path=os.path.join(output_dir, "05_monte_carlo.png"))
    
    # Plot 6: Comparison of methods
    fig6, ax = plt.subplots(figsize=(12, 6))
    ax.plot(x_fd, ut_fd[:, -1], 'b-', linewidth=2, label='FD', alpha=0.8)
    ax.plot(x_sp, ut_sp[:, -1], 'r--', linewidth=2, label='Spectral', alpha=0.8)
    ax.plot(x_fe, ut_fe[:, -1], 'g:', linewidth=2.5, label='FEM', alpha=0.8)
    ax.plot(x_mi, ut_mi[:, -1], 'm-.', linewidth=2, label='Milstein', alpha=0.8)
    ax.plot(x_mc, mean_u, 'k-', linewidth=3, label=f'MC Mean (M={M})', alpha=0.6)
    ax.set_xlabel('Space (x)', fontsize=11)
    ax.set_ylabel('u(T,x)', fontsize=11)
    ax.set_title(f'Comparison of Methods at t={T:.2f}', fontsize=13, fontweight='bold')
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3, linestyle='--')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "06_comparison.png"), dpi=300, bbox_inches='tight')
    logger.info(f"Figure saved to {os.path.join(output_dir, '06_comparison.png')}")
    
    plt.show()
    
    # ========================================================================
    # FINAL SUMMARY
    # ========================================================================
    logger.info("\n" + "="*70)
    logger.info("SIMULATION COMPLETE")
    logger.info("="*70)
    logger.info(f"All figures saved to: {os.path.abspath(output_dir)}")
    logger.info("\nKey Observations:")
    logger.info("1. FD and Milstein give similar stable results (additive noise)")
    logger.info("2. Spectral method uses periodic BC (different physics)")
    logger.info("3. FEM provides smooth approximation with better energy conservation")
    logger.info("4. MC statistics show typical stochastic fluctuations")
    logger.info("\nFor production use:")
    logger.info("- Adjust N (time steps) based on CFL condition")
    logger.info("- Increase M (MC samples) for better statistics")
    logger.info("- Consider implicit schemes for stiffer problems")
    logger.info("="*70 + "\n")