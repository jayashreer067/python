import pytest
from longest_sub import longest_contiguous, all_longest
from Bio.Seq import Seq


@pytest.mark.parametrize("seq, expected", [
    (Seq("a"), {"a": 1}),
    (Seq("aabbbbcc"), {"b": 4}),
    (Seq("aabb"), {"a": 2, "b": 2}),
    (Seq(""), {}),
])
def test_longest_contiguous(seq, expected):
    output = longest_contiguous(seq)

    assert output == expected


@pytest.mark.parametrize("records, expected", [
    ([{"a": 4}], {"a": 4}),
    ([{"a": 4}, {"b": 6, "a": 2}, {"b": 10}], {"a": 4, "b": 10}),
    ([], {}),
    ([{}], {}),
])
def test_all_longest(records, expected):
    output = all_longest(records)

    assert output == expected


def test_integration():
    output = all_longest(longest_contiguous(seq) for seq in [Seq("aaa"), Seq("aabbbbcc"), Seq("aabb")])

    assert output == {"a": 3, "b": 4}

