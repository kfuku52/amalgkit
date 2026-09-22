import pandas

from amalgkit.text_utils import normalize_unique_text


def test_normalizes_mixed_values_and_preserves_first_seen_order():
    values = [' b ', '', None, 'a', pandas.NA, 'b', float('nan'), ' ', 1, 2.5, 'a']
    assert normalize_unique_text(values) == ['b', 'a', '1', '2.5']
