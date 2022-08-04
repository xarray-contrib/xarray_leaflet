from pathlib import Path

import matplotlib.pyplot as plt

RESULTS = Path(__file__).parent / "results"
RESULTS.mkdir(exist_ok=True)
EXPECTED = Path(__file__).parent / "expected"
EXPECTED.mkdir(exist_ok=True)


def save_fig(da, name):
    result_path = RESULTS / name
    expected_path = EXPECTED / name
    da.plot()
    plt.savefig(result_path)
    plt.close()
    return result_path, expected_path
