# MH6502: Questions 1 and 2

Finite-difference boundary-value problems and SVD image approximation.

- `Q1_Q2_Notebook.ipynb`: current notebook, preserved unchanged.
- `Handwritten_Solutions.md`: step-by-step derivations and original sample-image results.
- `solve_q1_q2.py`: runnable Python script.
- `results/*.csv`: numerical tables; the SVD table uses the user's replacement image.

## Run

Open a terminal in this folder:

```sh
python -m pip install -r requirements.txt
python solve_q1_q2.py
jupyter notebook Q1_Q2_Notebook.ipynb
```

The default example uses Matplotlib's bundled Grace Hopper photograph. Running generates the figures referenced by the notebook and notes in `results/`. Image files are not included in this repository.

To use your own local image:

```sh
python solve_q1_q2.py --image "path/to/image.jpg"
```

In the notebook, change `q2_results = question2()` to `q2_results = question2("path/to/image.jpg")`.

The written SVD table describes the original sample. The committed SVD CSV describes the replacement image (rank 50 relative error approximately 2.98%). Running the example replaces the generated tables and figures. Open Jupyter with this assignment folder as the working directory.

Storage comparisons count scalars, not encoded file bytes.
