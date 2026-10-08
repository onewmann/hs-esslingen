# fftvel

Python implementation of the single-image velocity estimate (2-D FFT + Radon
transform of a pulsed-exposure image). It gives the same results as the
MATLAB/Octave functions in `../matlab` (checked by `tests/parity`).

```bash
pip install -e ".[plot]"
python -m fftvel simulate --vx 5 --vy -3 --plot fig.png
```

```python
from fftvel import estimate_velocity
est = estimate_velocity(bp, L=4)        # bp: mean of 2L+1 equally spaced exposures
print(est.angle_deg, est.speed, est.valid)
```

See the main README for the method, results and the live-camera example.
