"""Each module ships a --selftest of invariants and golden molecules; run them under pytest."""
from molcharge import charge, properties, validate


def test_charge_selftest():
    assert charge.selftest() == 0


def test_properties_selftest():
    assert properties.selftest() == 0


def test_validate_selftest():
    assert validate.selftest() == 0   # data-dependent gates skip when MOLCHARGE_DATA_DIR is absent
