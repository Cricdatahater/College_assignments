"""Run: python solve_q1_q2.py [--image your_photo.png]
Dependencies: numpy scipy matplotlib pillow
All output is written beside this script in results/.
"""
from pathlib import Path
import argparse
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded
from scipy.integrate import solve_bvp
from PIL import Image

OUT = Path(__file__).resolve().parent / 'results'

def save_csv(name, header, rows):
    with (OUT / name).open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)

def exact(x):
    return x * np.sin(2 * np.pi * x)

def manufactured_f(x):
    w = 2 * np.pi
    return ((w*w*x*(1+x*x)-2*x+x*np.exp(x))*np.sin(w*x)
            - 2*w*(1+2*x*x)*np.cos(w*x))

def fd(N, f):
    """N intervals, N-1 unknowns; zero Dirichlet boundary values."""
    x = np.linspace(0, 1, N+1)
    h = 1 / N
    xi = x[1:-1]
    am = 1 + (xi-h/2)**2
    ap = 1 + (xi+h/2)**2
    diag = (am+ap)/h**2 + np.exp(xi)
    lower, upper = -am[1:]/h**2, -ap[:-1]/h**2
    ab = np.zeros((3, N-1))
    ab[0, 1:], ab[1], ab[2, :-1] = upper, diag, lower
    b = np.asarray(f(xi))
    v = solve_banded((1, 1), ab, b)
    residual = diag*v-b
    residual[1:] += lower*v[:-1]
    residual[:-1] += upper*v[1:]
    return x, np.r_[0., v, 0.], np.max(np.abs(residual))

def question1():
    Ns = [20, 40, 80, 160, 320, 640]
    rows = []
    old = None
    for N in Ns:
        x, v, res = fd(N, manufactured_f)
        e = v-exact(x)
        einf, e2 = np.max(np.abs(e)), np.sqrt(np.sum(e*e)/N)
        p = np.log2(old/einf) if old is not None else float('nan')
        rows.append([N, 1/N, einf, e2, p, res])
        old = einf
    save_csv('q1_manufactured.csv', ['N', 'h', 'max_error', 'discrete_L2_error', 'order', 'residual_inf'], rows)
    assert 1.95 < rows[-1][4] < 2.05
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    x, v, _ = fd(40, manufactured_f)
    xx = np.linspace(0, 1, 1001)
    ax[0].plot(xx, exact(xx), label='Exact')
    ax[0].plot(x, v, '.', label='FD, N=40')
    ax[0].set(xlabel='x', ylabel='u', title='Manufactured solution')
    ax[0].legend()
    hs, errors = np.array(rows)[:, 1], np.array(rows)[:, 2]
    ax[1].loglog(hs, errors, 'o-', label='Maximum error')
    ax[1].loglog(hs, errors[-1]*(hs/hs[-1])**2, '--', label='Slope 2 reference')
    ax[1].set(xlabel='h', ylabel='Error', title='Second-order convergence')
    ax[1].legend()
    fig.tight_layout(); fig.savefig(OUT/'q1_manufactured.png', dpi=180); plt.close(fig)

    # Independent collocation check: y=(u, q), q=(1+x^2)u'.
    def ode(x, y):
        return np.vstack((y[1]/(1+x*x), np.exp(x)*y[0]-1))
    ref = solve_bvp(ode, lambda ya, yb: np.array([ya[0], yb[0]]),
                    np.linspace(0, 1, 51), np.zeros((2, 51)), tol=1e-10, max_nodes=20000)
    assert ref.success, ref.message
    ones = lambda x: np.ones_like(x)
    rows2 = []
    for N in Ns:
        x, v, res = fd(N, ones)
        xfine, vf, _ = fd(2*N, ones)
        delta = np.max(np.abs(v-vf[::2]))
        rows2.append([N, v.max(), x[v.argmax()], delta, np.max(np.abs(v-ref.sol(x)[0])), res])
        assert np.all(v[1:-1] > 0)
        assert np.all(v <= x*(1-x)*0.75 + 1e-12)
    save_csv('q1_f_equals_1.csv', ['N', 'max_u', 'x_at_grid_max', 'coarse_fine_difference', 'collocation_difference', 'residual_inf'], rows2)
    x, v, _ = fd(640, ones)
    save_csv('q1_f_equals_1_solution.csv', ['x', 'u'], zip(x, v))
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(x, v, label='FD solution, N=640')
    ax.plot(xx, xx*(1-xx)*0.75, '--', label='Upper bound 3x(1-x)/4')
    ax.set(xlabel='x', ylabel='u', title='Question 1(iv): f(x)=1')
    ax.legend(); fig.tight_layout(); fig.savefig(OUT/'q1_f_equals_1.png', dpi=180); plt.close(fig)
    print('\nQ1 manufactured: N, h, max error, L2 error, order, residual')
    for row in rows: print(row)
    print('\nQ1 f=1: N, max u, location, mesh difference, reference difference, residual')
    for row in rows2: print(row)
    return {'manufactured': rows, 'constant_rhs': rows2}

def question2(image_path=None):
    # Matplotlib ships this photograph; no external download is needed.
    source = Path(image_path) if image_path else Path(matplotlib.get_data_path())/'sample_data'/'grace_hopper.jpg'
    im = Image.open(source).convert('L')
    im.thumbnail((256, 256))
    if im.width == im.height:
        im = im.crop((0, 0, im.width-1, im.height))
    im.save(OUT/'input_grayscale.png')
    A = np.asarray(im, dtype=np.float64)/255.0
    m, n = A.shape
    U, s, Vt = np.linalg.svd(A, full_matrices=False)
    p = len(s)
    tails = np.sqrt(np.r_[np.cumsum((s*s)[::-1])[::-1], 0.])
    Ak = np.zeros_like(A)
    errors = []
    ranks = sorted(set(min(k, p) for k in [1, 5, 10, 20, 50, 100, p]))
    snapshots = {}
    rows = []
    for k in range(1, p+1):
        Ak += s[k-1]*np.outer(U[:, k-1], Vt[k-1])
        direct = np.linalg.norm(A-Ak, 'fro')
        errors.append(direct)
        rows.append([k, direct, tails[k], direct/np.linalg.norm(A, 'fro'), k*(m+n+1), m*n/(k*(m+n+1))])
        if k in ranks: snapshots[k] = Ak.copy()
    errors = np.array(errors)
    assert np.allclose(errors, tails[1:], rtol=1e-9, atol=1e-11)
    assert np.all(np.diff(errors) <= 1e-11)
    save_csv('q2_errors.csv', ['k', 'direct_Frobenius_error', 'singular_value_tail', 'relative_error', 'stored_scalars', 'scalar_compression_ratio'], rows)
    fig, axes = plt.subplots(2, 4, figsize=(12, 8))
    axes.flat[0].imshow(A, cmap='gray', vmin=0, vmax=1)
    axes.flat[0].set_title(f'Original ({m} x {n})')
    for ax, k in zip(list(axes.flat)[1:], ranks):
        # Clipping is for display only, never for evaluating the theorem.
        ax.imshow(np.clip(snapshots[k], 0, 1), cmap='gray', vmin=0, vmax=1)
        ax.set_title(f'k={k}; rel. error={errors[k-1]/np.linalg.norm(A):.3f}')
    for ax in axes.flat: ax.axis('off')
    fig.tight_layout(); fig.savefig(OUT/'q2_reconstructions.png', dpi=180); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4))
    ks = np.arange(1, p+1)
    ax.plot(ks, errors, label='Direct Frobenius error')
    ax.plot(ks[::6], tails[1:][::6], 'o', ms=3, label='Singular-value tail formula')
    ax.set(xlabel='Rank k', ylabel='E_k', title='SVD error identity')
    ax.legend(); fig.tight_layout(); fig.savefig(OUT/'q2_errors.png', dpi=180); plt.close(fig)
    summary = {'shape': [m, n], 'source': str(source), 'selected_results': [rows[k-1] for k in ranks],
               'max_identity_discrepancy': float(np.max(np.abs(errors-tails[1:])))}
    print('\nQ2:', json.dumps(summary, indent=2))
    return summary

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--image', type=Path, help='Optional replacement image')
    args = parser.parse_args()
    OUT.mkdir(exist_ok=True)
    summary = {'q1': question1(), 'q2': question2(args.image)}
    # Tables have NaN only for the undefined first convergence order.
    (OUT/'summary.json').write_text(json.dumps(summary, indent=2))
    print(f'\nResults saved in {OUT}')


if __name__ == '__main__':
    main()
