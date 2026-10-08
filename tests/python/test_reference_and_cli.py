import json

import numpy as np
import pytest
from scipy.io import savemat

from fftvel import phase_correlation, simulate_pulsed_image
from fftvel.__main__ import main


def test_phase_correlation():
    _, fr = simulate_pulsed_image(128, (3, -2), 1, rng=5, return_frames=True)
    dx, dy, _ = phase_correlation(fr[:, :, 0], fr[:, :, 1])
    assert abs(dx - 3) < 0.05 and abs(dy + 2) < 0.05
    _, fr = simulate_pulsed_image(128, (-1.5, 0.75), 1, rng=6, return_frames=True)
    dx, dy, _ = phase_correlation(fr[:, :, 0], fr[:, :, 2])
    assert abs(dx + 3) < 0.15 and abs(dy - 1.5) < 0.15


def test_cli_simulate(capsys, tmp_path):
    png = tmp_path / "fig.png"
    assert main(["simulate", "--vx", "4", "--vy", "1", "--plot", str(png)]) == 0
    out = capsys.readouterr().out
    result = json.loads(out[out.index("{"):out.rindex("}") + 1])
    assert abs(result["speed"] - 4.1231) < 0.05
    assert png.stat().st_size > 10_000


def test_cli_estimate_mat(repo, capsys):
    # bp01 of the parity cases is a pulsed image with v = (3, -2)
    path = repo / "tests" / "parity" / "cases.mat"
    assert main(["estimate", str(path), "--var", "bp01"]) == 0
    out = capsys.readouterr().out
    assert '"valid": true' in out


def test_cli_stack(tmp_path, capsys):
    _, fr = simulate_pulsed_image(256, (3.0, -2.0), 2, rng=3, return_frames=True)   # 5 frames
    savemat(tmp_path / "stack.mat", {"frames": fr})
    with open(tmp_path / "STACK.NPY", "wb") as f:     # np.save would append .npy
        np.save(f, fr)
    for name in ("stack.mat", "STACK.NPY"):
        assert main(["estimate", str(tmp_path / name), "--stack"]) == 0
        out = capsys.readouterr().out
        r = json.loads(out[out.index("{"):out.rindex("}") + 1])
        assert r["L"] == 2 and abs(r["speed"] / np.hypot(3, 2) - 1) < 0.01
    with pytest.raises(SystemExit):
        main(["estimate", str(tmp_path / "stack.mat")])               # 3-D without --stack
    np.save(tmp_path / "even.npy", fr[:, :, :4])
    with pytest.raises(SystemExit):
        main(["estimate", str(tmp_path / "even.npy"), "--stack"])     # even depth needs --L
