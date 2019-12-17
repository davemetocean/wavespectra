import numpy as np
from wavespectra.construct.helpers import (
    spread,
    check_coordinates,
    arrange_inputs,
    make_dataset,
)


def contruct_mlm(
    frequency,
    df,
    a1,
    b1,
    a2,
    b2,
    varianceDensity,
    direction,
    directionalSpread,
    coordinates=[],
    sumpart=True,
):
    """Constructs spectra using the Maximum Likelyhodd Method (MLM, e.g. Longuet-Higgins et al., 1963)."""
    # check_coordinates(tp, coordinates)

    # Arrange inputs
    frequency_m, df_m, a1_m, b1_m, a2_m, b2_m, varianceDensity_m, direction_m, directionalSpread_m = arrange_inputs(
        frequency, df, a1, b1, a2, b2, varianceDensity, direction, directionalSpread
    )
    adir = np.array(np.radians(direction))
    spec = np.zeros([len(direction), len(frequency)])
    for ii, f in enumerate(frequency):
        spec[ii, :] = varianceDensity[ii] * (
            (a1 * np.cos(adir))
            + (b1 * np.sin(adir))
            + (a2 * np.cos(2 * adir))
            + (b2 * np.sin(2 * adir))
        )

    ret = make_dataset(
        spec, frequency, direction, coordinates=coordinates
    )
    return ret


if __name__ == "__main__":
    from wavespectra.input.spotter import read_spotter

    filename = "../../tests/sample_files/spotter_20180214.json"
    dset = read_spotter(filename)
    contruct_mlm(
        dset.frequency,
        dset.df,
        dset.a1,
        dset.b1,
        dset.a2,
        dset.b2,
        dset.varianceDensity,
        dsets.direction,
        dset.directionalSpread,
    )
