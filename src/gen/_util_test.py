from gen._util import write_file, xor_sorted


def test_write_file_preserves_bytes(tmp_path):
    content = b"\x00\x01\xffbinary-content"
    path = tmp_path / "payload.bin"

    write_file(path, content)

    assert path.read_bytes() == content


def test_xor_sorted():
    assert xor_sorted([]) == []
    assert xor_sorted([2]) == [2]
    assert xor_sorted([2, 3]) == [2, 3]
    assert xor_sorted([3, 2]) == [2, 3]
    assert xor_sorted([2, 2]) == []
    assert xor_sorted([2, 2, 2]) == [2]
    assert xor_sorted([2, 2, 2, 2]) == []
    assert xor_sorted([2, 2, 3]) == [3]
    assert xor_sorted([3, 2, 2]) == [3]
    assert xor_sorted([2, 3, 2]) == [3]
    assert xor_sorted([2, 3, 3]) == [2]
    assert xor_sorted([2, 3, 5, 7, 11, 13, 5]) == [2, 3, 7, 11, 13]
