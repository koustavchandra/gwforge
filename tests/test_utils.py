import bilby
from GWForge import utils


def test_split_odd_even():
    odd, even = utils.split_odd_even(["a", "b", "c", "d", "e"])
    assert odd == ["b", "d"]
    assert even == ["a", "c", "e"]


def test_split_odd_even_single():
    # A single-element list must not produce an empty "odd" group that the
    # workflow generator would turn into an invalid "PARENT  CHILD x" DAG line.
    odd, even = utils.split_odd_even(["only"])
    assert odd == []
    assert even == ["only"]


def test_find_frame_files_one_sided_ranges(tmp_path):
    # Files span [1000, 1100] and [1100, 1200].
    (tmp_path / "H1-1000-100.gwf").touch()
    (tmp_path / "H1-1100-100.gwf").touch()

    # Only start_time given: previously raised TypeError on `x <= None`.
    _, paths = utils.find_frame_files(
        str(tmp_path), filePattern="*gwf", start_time=1050, end_time=None
    )
    assert len(paths) == 2

    # Only end_time given.
    _, paths = utils.find_frame_files(
        str(tmp_path), filePattern="*gwf", start_time=None, end_time=1050
    )
    assert len(paths) == 1

    # Both None: no filtering, return everything.
    _, paths = utils.find_frame_files(str(tmp_path), filePattern="*gwf")
    assert len(paths) == 2


def test_reference_priors_are_physical():
    # dec must be cosine-distributed on [-pi/2, pi/2]; theta_jn sine-distributed.
    assert isinstance(utils.reference_prior_dict["dec"], bilby.core.prior.Cosine)
    assert isinstance(utils.reference_prior_dict["theta_jn"], bilby.core.prior.Sine)
    assert utils.reference_prior_dict["dec"].minimum < 0
