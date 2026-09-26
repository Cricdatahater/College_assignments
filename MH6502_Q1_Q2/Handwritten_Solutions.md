# MH6502 — Questions 1 and 2

These are worked derivations to copy and explain by hand, together with results from the accompanying Python simulation. “First two questions” means all of Question 1(i)–(iv) and Question 2. The PDF supplies the mathematical problems; its administrative instructions are not additional requests from you.

## Question 1(i): Taylor derivation

Write $g(x)=a(x)u'(x)$. Assume $a,u$ have sufficiently many continuous derivatives for the following Taylor expansions and uniform remainders.

Taylor expansion at $x$ gives

$$g(x\pm h/2)=g(x)\pm\frac h2g'(x)+\frac{h^2}{8}g''(x)\pm\frac{h^3}{48}g'''(x)+O(h^4).$$

Subtract the two expansions and divide by $h$ (using one further Taylor term if writing the remainder as $O(h^4)$):

$$\frac{g(x+h/2)-g(x-h/2)}h=g'(x)+\frac{h^2}{24}g'''(x)+O(h^4).$$

Hence

$$P(x)=\frac{a(x+h/2)u'(x+h/2)-a(x-h/2)u'(x-h/2)}h+O(h^2).$$

Next expand $u$ about the midpoint $y=x+h/2$:

$$\frac{u(x+h)-u(x)}h=u'(x+h/2)+\frac{h^2}{24}u'''(x+h/2)+O(h^4).$$

Similarly,

$$\frac{u(x)-u(x-h)}h=u'(x-h/2)+\frac{h^2}{24}u'''(x-h/2)+O(h^4).$$

**Why the outer division by $h$ preserves second order:** the leading error in each approximate flux is $h^2a(y)u'''(y)/24$. The difference of these leading errors is $O(h^3)$ because the two midpoint locations are distance $h$ apart and $au'''$ is smooth. Dividing their difference by $h$ therefore gives $O(h^2)$. Merely substituting two unrelated $O(h^2)$ errors would not establish this conclusion.

Combining and collecting coefficients yields

$$\boxed{P(x)=\frac{a(x+h/2)u(x+h)-[a(x+h/2)+a(x-h/2)]u(x)+a(x-h/2)u(x-h)}{h^2}+O(h^2).}$$

As a check, if $D_hu$ denotes the displayed difference operator, direct expansion gives

$$D_hu=(au')'+h^2\left(\frac{au^{(4)}}{12}+\frac{a'u'''}6+\frac{a''u''}8+\frac{a'''u'}{24}\right)+O(h^4).$$

All coefficients and derivatives here are evaluated at $x$.

## Question 1(ii): Finite-difference scheme and linear system

Divide $[0,1]$ into $N$ equal intervals. Set $h=1/N$, $x_j=jh$, and $v_j\approx u(x_j)$. There are $N-1$ interior unknowns, with $v_0=v_N=0$.

Define $a_{j\pm1/2}=a(x_j\pm h/2)$, $c_j=c(x_j)$ and $f_j=f(x_j)$. Substitution in $-(au')'+cu=f$ gives, for $j=1,\ldots,N-1$,

$$\boxed{-\frac{a_{j-1/2}}{h^2}v_{j-1}+\left(\frac{a_{j-1/2}+a_{j+1/2}}{h^2}+c_j\right)v_j-\frac{a_{j+1/2}}{h^2}v_{j+1}=f_j.}$$

Let $\mathbf v=(v_1,\ldots,v_{N-1})^T$ and $\mathbf f=(f_1,\ldots,f_{N-1})^T$. Then $B\mathbf v=\mathbf f$, where

$$B_{jj}=\frac{a_{j-1/2}+a_{j+1/2}}{h^2}+c_j,\qquad B_{j,j+1}=B_{j+1,j}=-\frac{a_{j+1/2}}{h^2}.$$

All other entries vanish. Zero boundary values add nothing to the right-hand side. For nonzero boundary values $\alpha,\beta$, add $a_{1/2}\alpha/h^2$ to the first entry and $a_{N-1/2}\beta/h^2$ to the last.

**Uniqueness and stability.** For any interior vector $z$, extended by $z_0=z_N=0$,

$$z^TBz=\frac1{h^2}\sum_{j=0}^{N-1}a_{j+1/2}(z_{j+1}-z_j)^2+\sum_{j=1}^{N-1}c_jz_j^2.$$

Since $a\ge a_0>0$ and $c\ge0$, this expression is strictly positive for $z\ne0$. Therefore $B$ is symmetric positive definite and the scheme has a unique solution.

For completeness, its smallest eigenvalue is bounded below by

$$\lambda_{\min}(B)\ge\frac{4a_0}{h^2}\sin^2\!\left(\frac{\pi}{2N}\right)\ge4a_0.$$

The first inequality follows by comparison with the constant-coefficient Dirichlet difference matrix; the second uses $\sin t\ge2t/\pi$ for $0\le t\le\pi/2$. If $\tau=B\mathbf u-\mathbf f$ is the $O(h^2)$ truncation error, then $B(\mathbf v-\mathbf u)=-\tau$. Hence in the mesh norm $\|z\|_{2,h}=\sqrt{h\sum_jz_j^2}$,

$$\|\mathbf v-\mathbf u\|_{2,h}\le\frac1{4a_0}\|\tau\|_{2,h}=O(h^2).$$

The computations below also demonstrate second order in the maximum norm.

**Algorithm to write by hand:** choose $N$; compute grid and midpoint coefficients; assemble the three diagonals and the load vector; solve the tridiagonal system; append the two boundary zeros. A tridiagonal solver uses $O(N)$ operations and storage. The Python code uses `solve_banded`, without constructing a dense inverse.

## Question 1(iii): Manufactured solution

Take $u(x)=x\sin(2\pi x)$, $a(x)=1+x^2$, $c(x)=e^x$. The boundary values are zero because $\sin0=\sin2\pi=0$.

Differentiation gives

$$u'=\sin(2\pi x)+2\pi x\cos(2\pi x),$$

$$u''=4\pi\cos(2\pi x)-4\pi^2x\sin(2\pi x).$$

Therefore

$$f=-a'u'-au''+cu=-2xu'-(1+x^2)u''+e^xu,$$

so the required right-hand side is

$$\boxed{f(x)=\left[4\pi^2x(1+x^2)-2x+xe^x\right]\sin(2\pi x)-4\pi(1+2x^2)\cos(2\pi x).}$$

Use this expression in the scheme, and compare against the exact solution at the same grid points. Define

$$E_\infty(h)=\max_j|v_j-u(x_j)|,\qquad p=\log_2\frac{E_\infty(h)}{E_\infty(h/2)}.$$

Actual Python results:

| Intervals $N$ | $h$ | Maximum error | Observed order |
|---:|---:|---:|---:|
|20|0.05|7.113454e-3|—|
|40|0.025|1.768342e-3|2.00815|
|80|0.0125|4.421295e-4|1.99986|
|160|0.00625|1.104927e-4|2.00052|
|320|0.003125|2.762358e-5|1.99998|
|640|0.0015625|6.905740e-6|2.00003|

Halving $h$ reduces the error by approximately four, confirming second-order convergence.

![Manufactured-solution test](results/q1_manufactured.png)

## Question 1(iv): Solve with $f(x)=1$

Keep exactly the same coefficient matrix, but replace the load vector by $(1,\ldots,1)^T$. The manufactured solution from (iii) is no longer the exact solution.

For a small example that can be worked entirely by hand, choose $N=4$, $h=1/4$. The midpoint values of $a$ are $65/64,73/64,89/64,113/64$. The system is

$$\begin{pmatrix}
34.5+e^{1/4}&-18.25&0\\
-18.25&40.5+e^{1/2}&-22.25\\
0&-22.25&50.5+e^{3/4}
\end{pmatrix}\begin{pmatrix}v_1\\v_2\\v_3\end{pmatrix}=\begin{pmatrix}1\\1\\1\end{pmatrix}.$$

Let $d_1=34.5+e^{1/4}$, $d_2=40.5+e^{1/2}$, $d_3=50.5+e^{3/4}$. Row elimination gives $\tilde d_2=d_2-18.25^2/d_1$ and $\tilde b_2=1+18.25/d_1$, then $\tilde d_3=d_3-22.25^2/\tilde d_2$ and $\tilde b_3=1+22.25\tilde b_2/\tilde d_2$. Back-substitute $v_3=\tilde b_3/\tilde d_3$, $v_2=(\tilde b_2+22.25v_3)/\tilde d_2$, and $v_1=(1+18.25v_2)/d_1$. Numerically, the interior values are $0.07001429$, $0.08248729$, $0.05388643$. This is the same tridiagonal algorithm used on fine grids.

**Justification without assuming an elementary exact solution:**

1. The energy identity above proves uniqueness. For the continuous problem, subtract two solutions, multiply by their difference $w$, and integrate by parts: $\int_0^1[a(w')^2+cw^2]dx=0$. Thus $w=0$.
2. The solution is strictly positive in $(0,1)$. A nonpositive interior minimum would have $u'=0$, $u''\ge0$, giving $-(1+x^2)u''+e^xu\le0$, contradicting the right-hand side $1$. The discrete system has the corresponding maximum principle, and the computed interior values are positive.
3. A useful upper bound is $b(x)=3x(1-x)/4$. Indeed,

$$Lb=\frac32(1-2x+3x^2)+\frac34e^xx(1-x)\ge1,$$

because $1-2x+3x^2=3(x-1/3)^2+2/3$. Both $b$ and $u$ vanish at the endpoints, so comparison gives $0\le u(x)\le b(x)\le3/16$.
4. Check mesh refinement: $D_N=\max_j|v^{(N)}_j-v^{(2N)}_{2j}|$ should decrease by four. This does not require knowing the exact solution.
5. Independently check against a high-accuracy collocation solver. Set $q=(1+x^2)u'$, so $u'=q/(1+x^2)$ and $q'=e^xu-1$, with $u(0)=u(1)=0$. The code uses `solve_bvp` with tolerance $10^{-10}$. This is a numerical reference, not a symbolic exact answer or a rigorous error bound.

| $N$ | Grid maximum of $v$ | $D_N$ | Max difference from collocation |
|---:|---:|---:|---:|
|20|0.0838191762|7.37234e-6|9.83020e-6|
|40|0.0838175261|1.84337e-6|2.45786e-6|
|80|0.0838740398|4.60861e-7|6.14483e-7|
|160|0.0838739104|1.15216e-7|1.53622e-7|
|320|0.0838738780|2.88042e-8|3.84055e-8|
|640|0.0838738699|7.20126e-9|9.60144e-9|

At $N=640$, the grid maximum is about $0.08387387$ at $x=0.4375$. This is the maximum over grid nodes, not an exact location of the continuous maximum. The curve need not be symmetric because $a$ and $c$ are not symmetric about $x=1/2$. Small nonmonotonic changes in grid maxima arise because different meshes sample the peak differently. The residual $\|Bv-\mathbf1\|_\infty$ is $2.18\times10^{-11}$; residual alone would not establish discretization accuracy, which is why the refinement check is also included.

![Solution for unit forcing](results/q1_f_equals_1.png)

## Question 2: SVD image approximation — handwritten proof

Read an 8-bit grayscale image into a nonsquare matrix $G\in\mathbb R^{m\times n}$, then normalize $A=G/255$. Thus $0\le A_{ij}\le1$. Let $r=\operatorname{rank}(A)$ and write

$$A=\sum_{i=1}^r\sigma_i u_iv_i^T,\qquad\sigma_1\ge\cdots\ge\sigma_r>0,$$

where the $u_i$ and $v_i$ are orthonormal. For $1\le k\le r$, define

$$A_k=\sum_{i=1}^k\sigma_i u_iv_i^T.$$

The residual is $A-A_k=\sum_{i=k+1}^r\sigma_i u_iv_i^T$. Use the Frobenius inner product $\langle X,Y\rangle_F=\operatorname{tr}(X^TY)$:

$$\langle u_iv_i^T,u_jv_j^T\rangle_F=(u_i^Tu_j)(v_i^Tv_j)=\delta_{ij}.$$

Therefore all cross terms vanish when the norm is squared:

$$\begin{aligned}E_k^2&=\left\|\sum_{i=k+1}^r\sigma_iu_iv_i^T\right\|_F^2\\&=\sum_{i=k+1}^r\sum_{j=k+1}^r\sigma_i\sigma_j\delta_{ij}\\&=\sum_{i=k+1}^r\sigma_i^2.\end{aligned}$$

Taking square roots proves

$$\boxed{E_k=\|A-A_k\|_F=\sqrt{\sigma_{k+1}^2+\cdots+\sigma_r^2}.}$$

This is an **equality**, stronger than an upper-bound estimate. Consequently $E_k^2-E_{k+1}^2=\sigma_{k+1}^2\ge0$, so error decreases with $k$; $E_r=0$ in exact arithmetic. Also $\|A\|_F^2=\sum_i\sigma_i^2$, hence

$$\frac{E_k}{\|A\|_F}=\sqrt{1-\frac{\sum_{i=1}^k\sigma_i^2}{\sum_{i=1}^r\sigma_i^2}}.$$

The truncated SVD is the best approximation of rank at most $k$ in Frobenius norm (Eckart–Young theorem). Larger $k$ retains more image detail.

**Small example to calculate by hand.** Consider the already normalized nonsquare matrix

$$A=\begin{pmatrix}1&0&0\\0&1/2&0\end{pmatrix}.$$

Here $AA^T=\operatorname{diag}(1,1/4)$, so its singular values are $1$ and $1/2$. Choose $U=I_2$, $\Sigma=\operatorname{diag}(1,1/2)$, and $V^T=\begin{pmatrix}1&0&0\\0&1&0\end{pmatrix}$. Then

$$A_1=\begin{pmatrix}1&0&0\\0&0&0\end{pmatrix},\quad E_1=1/2=\sqrt{\sigma_2^2};\qquad A_2=A,\quad E_2=0.$$

This illustrates the same calculation used for a large image; a full image SVD is computed numerically.

## Question 2: Python experiment and interpretation

The reproducible example uses Matplotlib's bundled Grace Hopper photograph, converted to grayscale and resized to $256\times218$. The script then saves and reads its pixel data as an image matrix normalized by 255. The saved `results/input_grayscale.png` is the exact input used. Supply `--image path/to/image.png` to run the experiment with your own image. A square input is cropped by one column so the assignment's nonsquare condition holds.

Algorithm: read and normalize; compute `U, s, Vt = np.linalg.svd(A, full_matrices=False)`; form `(U[:, :k] * s[:k]) @ Vt[:k, :]`; independently evaluate `np.linalg.norm(A-Ak, 'fro')` and `np.sqrt(np.sum(s[k:]**2))`; plot errors for every $k=1,\ldots,\min(m,n)$; display selected reconstructions. The implementation accumulates the rank-one terms to check every $k$ efficiently. Tiny singular values are retained, avoiding an arbitrary numerical-rank threshold.

| $k$ | Direct $E_k$ | Tail formula | Relative error | Scalar-count compression factor |
|---:|---:|---:|---:|---:|
|1|48.407298|48.407298|50.9091%|117.491|
|5|31.413778|31.413778|33.0373%|23.498|
|10|23.213898|23.213898|24.4137%|11.749|
|20|15.735348|15.735348|16.5486%|5.875|
|50|7.071610|7.071610|7.4371%|2.350|
|100|2.498058|2.498058|2.6272%|1.175|
|218|3.43e-13|0|3.61e-13%|0.539|

The largest absolute discrepancy between the two error calculations over every $k$ is $3.44\times10^{-13}$, consistent with floating-point roundoff. At full rank, the formula's tail is empty and exactly zero, while a floating-point reconstruction has a tiny nonzero residual.

The original matrix stores $mn$ scalars; factored $A_k$ stores $mk+k+nk=k(m+n+1)$ scalars. The table reports $mn/[k(m+n+1)]$. This is a scalar-count comparison, not a measured file compression ratio: float storage, image encoding, and quantization change byte counts. Actual scalar savings occur only when $k<mn/(m+n+1)$, here approximately $117.49$.

The low-rank approximation can have entries outside $[0,1]$. Clip only for display; calculate errors using the unmodified $A_k$, or the proved identity will no longer apply.

![SVD reconstructions](results/q2_reconstructions.png)

![SVD error verification](results/q2_errors.png)

## Running the supplied files

In a terminal opened in this folder:

```bash
python -m pip install numpy scipy matplotlib pillow
python solve_q1_q2.py
# Optional: use your own image
python solve_q1_q2.py --image "my_image.png"
```

The script writes figures, complete CSV tables, and a JSON summary into `results/`. It verifies the manufactured solution's observed order, positivity and an upper bound for the unit-forcing solution, successful independent collocation, monotonic SVD errors, and agreement with the singular-value formula. `Q1_Q2_Notebook.ipynb` contains these derivations and the same executable code for Jupyter.
