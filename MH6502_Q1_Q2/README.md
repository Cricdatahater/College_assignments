# MH6502: Questions 1 and 2

Finite-difference boundary-value problems and SVD image approximation.

- `Q1_Q2_Notebook.ipynb`: current notebook, preserved unchanged.
- `Handwritten_Solutions.md`: derivations and original sample-image results.
- `solve_q1_q2.py`: runnable Python script.
- `results/`: latest plots and tables, using the supplied replacement photograph for Question 2.

## Run

From this folder:

```sh
python -m pip install -r requirements.txt
python solve_q1_q2.py --image "results/jesus death img.jpg"
```

To use the same image in the notebook, change `q2_results = question2()` to `q2_results = question2("results/jesus death img.jpg")` before running its code cell. Open Jupyter with this folder as the working directory. The notebook defaults to Matplotlib's Grace Hopper sample; its written SVD table describes that original sample. The current CSV files and plots describe the replacement image (rank 50 relative error approximately 2.98%). Running either example replaces generated files in `results/`.

The error formula is checked against direct Frobenius errors at every rank. Storage comparisons count scalars, not encoded image-file bytes.
